#include "dixon_wrapper.h"

#include <windows.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define DIXON_PATH_CAP 4096
#define DIXON_OUTPUT_LIMIT (4 * 1024 * 1024)
#define DIXON_STARTUP_TIMEOUT_MS 30000
#define DIXON_WAIT_SLICE_MS 250

typedef struct {
    char *data;
    size_t len;
    size_t cap;
} string_builder_t;

static char *dup_string(const char *text)
{
    size_t len;
    char *copy;

    if (!text) {
        copy = (char *) malloc(1);
        if (copy) copy[0] = '\0';
        return copy;
    }

    len = strlen(text);
    copy = (char *) malloc(len + 1);
    if (!copy) return NULL;
    memcpy(copy, text, len + 1);
    return copy;
}

static void set_error(char *buffer, size_t size, const char *message)
{
    if (!buffer || size == 0) return;
    if (!message) message = "Unknown error";
    snprintf(buffer, size, "%s", message);
}

static void set_win32_error(char *buffer, size_t size, const char *prefix, DWORD error_code)
{
    char system_message[512];
    DWORD written = FormatMessageA(
        FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_IGNORE_INSERTS,
        NULL,
        error_code,
        MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
        system_message,
        (DWORD) sizeof(system_message),
        NULL);

    if (written == 0) {
        snprintf(system_message, sizeof(system_message), "Windows error %lu", (unsigned long) error_code);
    } else {
        while (written > 0 &&
               (system_message[written - 1] == '\r' || system_message[written - 1] == '\n')) {
            system_message[--written] = '\0';
        }
    }

    if (buffer && size > 0) {
        snprintf(buffer, size, "%s: %s", prefix, system_message);
    }
}

static int sb_ensure(string_builder_t *sb, size_t extra)
{
    size_t needed;
    size_t new_cap;
    char *new_data;

    needed = sb->len + extra + 1;
    if (needed <= sb->cap) return 1;

    new_cap = sb->cap ? sb->cap : 256;
    while (new_cap < needed) new_cap *= 2;

    new_data = (char *) realloc(sb->data, new_cap);
    if (!new_data) return 0;

    sb->data = new_data;
    sb->cap = new_cap;
    return 1;
}

static int sb_append_n(string_builder_t *sb, const char *text, size_t len)
{
    if (!sb_ensure(sb, len)) return 0;
    memcpy(sb->data + sb->len, text, len);
    sb->len += len;
    sb->data[sb->len] = '\0';
    return 1;
}

static int sb_append(string_builder_t *sb, const char *text)
{
    if (!text) text = "";
    return sb_append_n(sb, text, strlen(text));
}

static void sb_free(string_builder_t *sb)
{
    free(sb->data);
    sb->data = NULL;
    sb->len = 0;
    sb->cap = 0;
}

static int append_backslashes(string_builder_t *sb, size_t count)
{
    size_t i;

    for (i = 0; i < count; ++i) {
        if (!sb_append(sb, "\\")) return 0;
    }
    return 1;
}

/*
 * Quote a single Windows command-line argument using the rules expected by
 * CreateProcess / CommandLineToArgvW style parsing.
 */
static int append_quoted_arg(string_builder_t *sb, const char *arg)
{
    const char *p;
    int needs_quotes = 0;
    size_t backslashes = 0;

    if (!arg || arg[0] == '\0') {
        return sb_append(sb, "\"\"");
    }

    for (p = arg; *p; ++p) {
        if (*p == ' ' || *p == '\t' || *p == '"' || *p == '\n' || *p == '\r') {
            needs_quotes = 1;
            break;
        }
    }

    if (!needs_quotes) return sb_append(sb, arg);
    if (!sb_append(sb, "\"")) return 0;

    for (p = arg; *p; ++p) {
        if (*p == '\\') {
            backslashes++;
            continue;
        }

        if (*p == '"') {
            if (!append_backslashes(sb, backslashes * 2 + 1)) return 0;
            if (!sb_append_n(sb, "\"", 1)) return 0;
            backslashes = 0;
            continue;
        }

        if (!append_backslashes(sb, backslashes)) return 0;
        backslashes = 0;

        if (!sb_append_n(sb, p, 1)) return 0;
    }

    if (!append_backslashes(sb, backslashes * 2)) return 0;
    return sb_append(sb, "\"");
}

static int append_arg(string_builder_t *sb, const char *arg)
{
    if (sb->len > 0 && !sb_append(sb, " ")) return 0;
    return append_quoted_arg(sb, arg);
}

static const char *skip_first_arg(const char *command_line)
{
    const char *p = command_line;

    if (!p) return "";

    while (*p == ' ' || *p == '	') p++;

    if (*p == '"') {
        p++;
        while (*p) {
            if (*p == '"') {
                p++;
                break;
            }
            if (*p == '\\' && p[1] != '\0') p++;
            p++;
        }
    } else {
        while (*p && *p != ' ' && *p != '	') p++;
    }

    while (*p == ' ' || *p == '	') p++;
    return p;
}

static char *make_display_command_line(const char *actual_command_line)
{
    string_builder_t sb = {0};
    const char *tail = skip_first_arg(actual_command_line);

    if (!sb_append(&sb, "./dixon.exe")) {
        sb_free(&sb);
        return NULL;
    }

    if (*tail) {
        if (!sb_append(&sb, " ") || !sb_append(&sb, tail)) {
            sb_free(&sb);
            return NULL;
        }
    }

    return sb.data;
}

static char *compose_combined_output(const char *display_command_line,
                                     const char *stdout_text)
{
    string_builder_t sb = {0};
    const char *body = stdout_text ? stdout_text : "";

    if (display_command_line && display_command_line[0] != '\0') {
        if (!sb_append(&sb, "Command: ") ||
            !sb_append(&sb, display_command_line) ||
            !sb_append(&sb, "\r\n\r\n")) {
            sb_free(&sb);
            return NULL;
        }
    }

    if (!sb_append(&sb, body)) {
        sb_free(&sb);
        return NULL;
    }

    if (!sb.data) return dup_string("");
    return sb.data;
}

static char *trim_copy(const char *text)
{
    const char *start;
    const char *end;
    size_t len;
    char *copy;

    if (!text) return dup_string("");

    start = text;
    while (*start == ' ' || *start == '\t' || *start == '\r' || *start == '\n') start++;

    end = text + strlen(text);
    while (end > start &&
           (end[-1] == ' ' || end[-1] == '\t' || end[-1] == '\r' || end[-1] == '\n')) {
        end--;
    }

    len = (size_t) (end - start);
    copy = (char *) malloc(len + 1);
    if (!copy) return NULL;
    memcpy(copy, start, len);
    copy[len] = '\0';
    return copy;
}

static char *normalize_multiline_csv(const char *text)
{
    string_builder_t sb = {0};
    const char *line = text;
    const char *p = text;

    if (!text) return dup_string("");

    while (1) {
        if (*p == '\r' || *p == '\n' || *p == '\0') {
            size_t line_len = (size_t) (p - line);
            char *raw = (char *) malloc(line_len + 1);
            char *trimmed;

            if (!raw) {
                sb_free(&sb);
                return NULL;
            }

            memcpy(raw, line, line_len);
            raw[line_len] = '\0';
            trimmed = trim_copy(raw);
            free(raw);

            if (!trimmed) {
                sb_free(&sb);
                return NULL;
            }

            if (trimmed[0] != '\0') {
                if (sb.len > 0 && !sb_append(&sb, ", ")) {
                    free(trimmed);
                    sb_free(&sb);
                    return NULL;
                }
                if (!sb_append(&sb, trimmed)) {
                    free(trimmed);
                    sb_free(&sb);
                    return NULL;
                }
            }

            free(trimmed);

            if (*p == '\0') break;
            if (*p == '\r' && p[1] == '\n') p++;
            p++;
            line = p;
            continue;
        }

        p++;
    }

    if (!sb.data) return dup_string("");
    return sb.data;
}

void dixon_gui_response_init(dixon_gui_response_t *response)
{
    if (!response) return;
    response->exit_code = -1;
    response->resolved_dixon_path = NULL;
    response->command_line = NULL;
    response->stdout_text = NULL;
    response->stderr_text = NULL;
    response->combined_output = NULL;
}

void dixon_gui_response_clear(dixon_gui_response_t *response)
{
    if (!response) return;
    free(response->resolved_dixon_path);
    free(response->command_line);
    free(response->stdout_text);
    free(response->stderr_text);
    free(response->combined_output);
    dixon_gui_response_init(response);
}

int dixon_gui_find_default_cli(char *buffer, size_t buffer_size)
{
    DWORD len;
    char path[DIXON_PATH_CAP];
    char *slash;

    if (!buffer || buffer_size == 0) return 0;

    len = GetModuleFileNameA(NULL, path, (DWORD) sizeof(path));
    if (len == 0 || len >= sizeof(path)) return 0;

    slash = strrchr(path, '\\');
    if (!slash) slash = strrchr(path, '/');
    if (slash) {
        slash[1] = '\0';
    } else {
        path[0] = '\0';
    }

    if (snprintf(buffer, buffer_size, "%sdixon.exe", path) >= (int) buffer_size) {
        return 0;
    }

    return 1;
}

static int validate_request(const dixon_gui_request_t *request,
                            const char *resolved_path,
                            char *error_message,
                            size_t error_message_size)
{
    DWORD attrs;

    if (!request) {
        set_error(error_message, error_message_size, "Request is null.");
        return 0;
    }

    attrs = GetFileAttributesA(resolved_path);
    if (attrs == INVALID_FILE_ATTRIBUTES || (attrs & FILE_ATTRIBUTE_DIRECTORY)) {
        set_error(error_message, error_message_size,
                  "Cannot find dixon.exe. Put it next to the GUI or choose it manually.");
        return 0;
    }

    if (!request->field || request->field[0] == '\0') {
        set_error(error_message, error_message_size, "Field is required.");
        return 0;
    }

    switch (request->mode) {
        case DIXON_GUI_MODE_SOLVE:
            if (!request->polynomials || request->polynomials[0] == '\0') {
                set_error(error_message, error_message_size, "Polynomials are required.");
                return 0;
            }
            return 1;
        case DIXON_GUI_MODE_RESULTANT:
            if (!request->polynomials || request->polynomials[0] == '\0') {
                set_error(error_message, error_message_size, "Polynomials are required.");
                return 0;
            }
            if (!request->eliminate_vars || request->eliminate_vars[0] == '\0') {
                set_error(error_message, error_message_size,
                          "Elimination variables are required for this mode.");
                return 0;
            }
            return 1;
        case DIXON_GUI_MODE_COMPLEXITY:
            if (!request->polynomials || request->polynomials[0] == '\0') {
                set_error(error_message, error_message_size, "Polynomials are required.");
                return 0;
            }
            if (!request->eliminate_vars || request->eliminate_vars[0] == '\0') {
                set_error(error_message, error_message_size,
                          "Elimination variables are required for this mode.");
                return 0;
            }
            return 1;
        case DIXON_GUI_MODE_IDEAL:
            if (!request->polynomials || request->polynomials[0] == '\0') {
                set_error(error_message, error_message_size, "Polynomials are required.");
                return 0;
            }
            if (!request->eliminate_vars || request->eliminate_vars[0] == '\0') {
                set_error(error_message, error_message_size,
                          "Elimination variables are required for ideal reduction.");
                return 0;
            }
            if (!request->ideal_generators || request->ideal_generators[0] == '\0') {
                set_error(error_message, error_message_size,
                          "Ideal generators are required for ideal reduction.");
                return 0;
            }
            return 1;
        case DIXON_GUI_MODE_RANDOM:
            if (!request->random_degrees || request->random_degrees[0] == '\0') {
                set_error(error_message, error_message_size,
                          "Degree list is required for random mode.");
                return 0;
            }
            if (request->random_mode < DIXON_GUI_RANDOM_RESULTANT ||
                request->random_mode > DIXON_GUI_RANDOM_COMPLEXITY) {
                set_error(error_message, error_message_size,
                          "Random mode selection is invalid.");
                return 0;
            }
            return 1;
        default:
            set_error(error_message, error_message_size, "Unsupported mode.");
            return 0;
    }
}

static int build_temp_output_path(char *buffer,
                                  size_t buffer_size,
                                  char *error_message,
                                  size_t error_message_size)
{
    char temp_dir[MAX_PATH];
    char temp_file[MAX_PATH];
    DWORD dir_len;

    dir_len = GetTempPathA((DWORD) sizeof(temp_dir), temp_dir);
    if (dir_len == 0 || dir_len >= sizeof(temp_dir)) {
        set_win32_error(error_message, error_message_size,
                        "Failed to locate the temporary directory", GetLastError());
        return 0;
    }

    if (GetTempFileNameA(temp_dir, "dix", 0, temp_file) == 0) {
        set_win32_error(error_message, error_message_size,
                        "Failed to create a temporary output file", GetLastError());
        return 0;
    }

    if (snprintf(buffer, buffer_size, "%s", temp_file) >= (int) buffer_size) {
        DeleteFileA(temp_file);
        set_error(error_message, error_message_size, "Temporary output path is too long.");
        return 0;
    }

    return 1;
}

static int get_file_size_ull(const char *path, unsigned long long *size_out)
{
    WIN32_FILE_ATTRIBUTE_DATA info;

    if (!GetFileAttributesExA(path, GetFileExInfoStandard, &info)) {
        return 0;
    }

    *size_out = (((unsigned long long) info.nFileSizeHigh) << 32) |
                (unsigned long long) info.nFileSizeLow;
    return 1;
}

static int read_file_into_buffer_limited(const char *path,
                                         char **output,
                                         char *error_message,
                                         size_t error_message_size)
{
    string_builder_t sb = {0};
    char chunk[4096];
    size_t total = 0;
    int truncated = 0;
    FILE *fp;

    fp = fopen(path, "rb");
    if (!fp) {
        set_error(error_message, error_message_size,
                  "Failed to open the temporary GUI output file.");
        return 0;
    }

    while (!truncated) {
        size_t bytes_read = fread(chunk, 1, sizeof(chunk), fp);
        if (bytes_read == 0) break;

        if (total + bytes_read > DIXON_OUTPUT_LIMIT) {
            bytes_read = DIXON_OUTPUT_LIMIT - total;
            truncated = 1;
        }

        if (bytes_read > 0 && !sb_append_n(&sb, chunk, bytes_read)) {
            fclose(fp);
            sb_free(&sb);
            set_error(error_message, error_message_size,
                      "Out of memory while capturing GUI output.");
            return 0;
        }

        total += bytes_read;
        if (total >= DIXON_OUTPUT_LIMIT) truncated = 1;
    }

    if (ferror(fp)) {
        fclose(fp);
        sb_free(&sb);
        set_error(error_message, error_message_size,
                  "Failed while reading the temporary GUI output file.");
        return 0;
    }

    fclose(fp);

    if (truncated) {
        char note[128];
        snprintf(note, sizeof(note), "\n\n[GUI output truncated at %d bytes]\n", DIXON_OUTPUT_LIMIT);
        if (!sb_append(&sb, note)) {
            sb_free(&sb);
            set_error(error_message, error_message_size,
                      "Out of memory while finalizing GUI output.");
            return 0;
        }
    }

    if (!sb.data) sb.data = dup_string("");
    if (!sb.data) {
        set_error(error_message, error_message_size,
                  "Out of memory while finalizing GUI output.");
        return 0;
    }

    *output = sb.data;
    return 1;
}

int dixon_gui_run_request(const dixon_gui_request_t *request,
                          dixon_gui_response_t *response,
                          char *error_message,
                          size_t error_message_size)
{
    char default_path[DIXON_PATH_CAP];
    const char *path_to_use;
    char *trimmed_path = NULL;
    char *field = NULL;
    char *polynomials = NULL;
    char *vars = NULL;
    char *ideal = NULL;
    char *omega = NULL;
    char *random_degrees = NULL;
    string_builder_t command = {0};
    SECURITY_ATTRIBUTES sa = {0};
    HANDLE output_file = INVALID_HANDLE_VALUE;
    STARTUPINFOA si;
    PROCESS_INFORMATION pi;
    DWORD exit_code = 0;
    DWORD wait_status;
    DWORD startup_waited = 0;
    unsigned long long output_size = 0;
    char output_path[DIXON_PATH_CAP] = {0};
    char working_dir[DIXON_PATH_CAP];
    char *last_slash;
    char *actual_command_line = NULL;
    int saw_child_output = 0;
    int ok = 0;

    if (!response) {
        set_error(error_message, error_message_size, "Response object is required.");
        return 0;
    }

    dixon_gui_response_clear(response);
    ZeroMemory(&pi, sizeof(pi));

    if (request->dixon_path && request->dixon_path[0] != '\0') {
        trimmed_path = trim_copy(request->dixon_path);
        if (!trimmed_path) {
            set_error(error_message, error_message_size,
                      "Out of memory while processing the CLI path.");
            goto cleanup;
        }
        path_to_use = trimmed_path;
    } else {
        if (!dixon_gui_find_default_cli(default_path, sizeof(default_path))) {
            set_error(error_message, error_message_size,
                      "Failed to resolve the default dixon.exe path.");
            goto cleanup;
        }
        path_to_use = default_path;
    }

    if (!validate_request(request, path_to_use, error_message, error_message_size)) {
        goto cleanup;
    }

    field = trim_copy(request->field);
    if (!field) {
        set_error(error_message, error_message_size, "Out of memory while copying the field input.");
        goto cleanup;
    }

    polynomials = normalize_multiline_csv(request->polynomials);
    if (!polynomials) {
        set_error(error_message, error_message_size, "Out of memory while normalizing the polynomial input.");
        goto cleanup;
    }

    if (request->eliminate_vars) {
        vars = normalize_multiline_csv(request->eliminate_vars);
        if (!vars) {
            set_error(error_message, error_message_size, "Out of memory while normalizing the elimination variables.");
            goto cleanup;
        }
    }

    if (request->ideal_generators) {
        ideal = normalize_multiline_csv(request->ideal_generators);
        if (!ideal) {
            set_error(error_message, error_message_size, "Out of memory while normalizing the ideal generators.");
            goto cleanup;
        }
    }

    if (request->omega) {
        omega = trim_copy(request->omega);
        if (!omega) {
            set_error(error_message, error_message_size, "Out of memory while copying the omega parameter.");
            goto cleanup;
        }
    }

    if (request->random_degrees) {
        random_degrees = trim_copy(request->random_degrees);
        if (!random_degrees) {
            set_error(error_message, error_message_size, "Out of memory while copying the random degree list.");
            goto cleanup;
        }
    }

    if (!append_arg(&command, path_to_use)) {
        set_error(error_message, error_message_size,
                  "Out of memory while appending the dixon.exe path to the command line.");
        goto cleanup;
    }

    switch (request->mode) {
        case DIXON_GUI_MODE_SOLVE:
            if (!append_arg(&command, request->solve_verbose ? "--solve-verbose" : "--solve")) {
                set_error(error_message, error_message_size, "Out of memory while building the solve command line.");
                goto cleanup;
            }
            if (!append_arg(&command, polynomials)) {
                set_error(error_message, error_message_size, "Out of memory while appending polynomials to the command line.");
                goto cleanup;
            }
            if (!append_arg(&command, field)) {
                set_error(error_message, error_message_size, "Out of memory while appending the field to the command line.");
                goto cleanup;
            }
            break;
        case DIXON_GUI_MODE_RESULTANT:
            if (!append_arg(&command, polynomials)) {
                set_error(error_message, error_message_size, "Out of memory while appending polynomials to the command line.");
                goto cleanup;
            }
            if (!append_arg(&command, vars ? vars : "")) {
                set_error(error_message, error_message_size, "Out of memory while appending elimination variables to the command line.");
                goto cleanup;
            }
            if (!append_arg(&command, field)) {
                set_error(error_message, error_message_size, "Out of memory while appending the field to the command line.");
                goto cleanup;
            }
            break;
        case DIXON_GUI_MODE_COMPLEXITY:
            if (!append_arg(&command, "--comp")) {
                set_error(error_message, error_message_size, "Out of memory while building the complexity command line.");
                goto cleanup;
            }
            if (omega && omega[0] != '\0') {
                if (!append_arg(&command, "--omega")) {
                    set_error(error_message, error_message_size, "Out of memory while adding --omega to the command line.");
                    goto cleanup;
                }
                if (!append_arg(&command, omega)) {
                    set_error(error_message, error_message_size, "Out of memory while appending omega to the command line.");
                    goto cleanup;
                }
            }
            if (!append_arg(&command, polynomials)) {
                set_error(error_message, error_message_size, "Out of memory while appending polynomials to the command line.");
                goto cleanup;
            }
            if (!append_arg(&command, vars ? vars : "")) {
                set_error(error_message, error_message_size, "Out of memory while appending elimination variables to the command line.");
                goto cleanup;
            }
            if (!append_arg(&command, field)) {
                set_error(error_message, error_message_size, "Out of memory while appending the field to the command line.");
                goto cleanup;
            }
            break;
        case DIXON_GUI_MODE_IDEAL:
            if (!append_arg(&command, "--ideal")) {
                set_error(error_message, error_message_size, "Out of memory while building the ideal command line.");
                goto cleanup;
            }
            if (!append_arg(&command, ideal ? ideal : "")) {
                set_error(error_message, error_message_size, "Out of memory while appending ideal generators to the command line.");
                goto cleanup;
            }
            if (!append_arg(&command, polynomials)) {
                set_error(error_message, error_message_size, "Out of memory while appending polynomials to the command line.");
                goto cleanup;
            }
            if (!append_arg(&command, vars ? vars : "")) {
                set_error(error_message, error_message_size, "Out of memory while appending elimination variables to the command line.");
                goto cleanup;
            }
            if (!append_arg(&command, field)) {
                set_error(error_message, error_message_size, "Out of memory while appending the field to the command line.");
                goto cleanup;
            }
            break;
        case DIXON_GUI_MODE_RANDOM:
            if (request->random_mode == DIXON_GUI_RANDOM_SOLVE) {
                if (!append_arg(&command, "--solve")) {
                    set_error(error_message, error_message_size, "Out of memory while building the random solve command line.");
                    goto cleanup;
                }
            } else if (request->random_mode == DIXON_GUI_RANDOM_COMPLEXITY) {
                if (!append_arg(&command, "--comp")) {
                    set_error(error_message, error_message_size, "Out of memory while building the random complexity command line.");
                    goto cleanup;
                }
                if (omega && omega[0] != '\0') {
                    if (!append_arg(&command, "--omega")) {
                        set_error(error_message, error_message_size, "Out of memory while adding --omega to the random command line.");
                        goto cleanup;
                    }
                    if (!append_arg(&command, omega)) {
                        set_error(error_message, error_message_size, "Out of memory while appending omega to the random command line.");
                        goto cleanup;
                    }
                }
            }
            if (!append_arg(&command, "--random")) {
                set_error(error_message, error_message_size, "Out of memory while building the random command line.");
                goto cleanup;
            }
            if (!append_arg(&command, random_degrees ? random_degrees : "")) {
                set_error(error_message, error_message_size, "Out of memory while appending the random degree list to the command line.");
                goto cleanup;
            }
            if (!append_arg(&command, field)) {
                set_error(error_message, error_message_size, "Out of memory while appending the field to the random command line.");
                goto cleanup;
            }
            break;
        default:
            set_error(error_message, error_message_size, "Unsupported mode.");
            goto cleanup;
    }

    if (snprintf(working_dir, sizeof(working_dir), "%s", path_to_use) >= (int) sizeof(working_dir)) {
        set_error(error_message, error_message_size, "The dixon.exe path is too long.");
        goto cleanup;
    }
    last_slash = strrchr(working_dir, '\\');
    if (!last_slash) last_slash = strrchr(working_dir, '/');
    if (last_slash) {
        *last_slash = '\0';
    } else {
        snprintf(working_dir, sizeof(working_dir), ".");
    }

    if (!build_temp_output_path(output_path, sizeof(output_path), error_message, error_message_size)) {
        goto cleanup;
    }

    sa.nLength = sizeof(sa);
    sa.lpSecurityDescriptor = NULL;
    sa.bInheritHandle = TRUE;

    output_file = CreateFileA(output_path,
                              GENERIC_WRITE,
                              FILE_SHARE_READ | FILE_SHARE_WRITE,
                              &sa,
                              CREATE_ALWAYS,
                              FILE_ATTRIBUTE_TEMPORARY,
                              NULL);
    if (output_file == INVALID_HANDLE_VALUE) {
        set_win32_error(error_message, error_message_size,
                        "Failed to create the temporary GUI output file", GetLastError());
        goto cleanup;
    }

    ZeroMemory(&si, sizeof(si));
    si.cb = sizeof(si);
    si.dwFlags = STARTF_USESTDHANDLES;
    si.hStdInput = GetStdHandle(STD_INPUT_HANDLE);
    si.hStdOutput = output_file;
    si.hStdError = output_file;

    response->resolved_dixon_path = dup_string(path_to_use);
    actual_command_line = dup_string(command.data ? command.data : "");
    response->command_line = make_display_command_line(actual_command_line);
    if (!response->resolved_dixon_path || !actual_command_line || !response->command_line) {
        set_error(error_message, error_message_size,
                  "Out of memory while preparing the child process command line.");
        goto cleanup;
    }

    if (!CreateProcessA(
            response->resolved_dixon_path,
            actual_command_line,
            NULL,
            NULL,
            TRUE,
            CREATE_NO_WINDOW,
            NULL,
            working_dir,
            &si,
            &pi)) {
        set_win32_error(error_message, error_message_size,
                        "Failed to start dixon.exe", GetLastError());
        goto cleanup;
    }

    CloseHandle(output_file);
    output_file = INVALID_HANDLE_VALUE;

    for (;;) {
        wait_status = WaitForSingleObject(pi.hProcess, DIXON_WAIT_SLICE_MS);
        if (wait_status == WAIT_OBJECT_0) break;
        if (wait_status != WAIT_TIMEOUT) {
            set_win32_error(error_message, error_message_size,
                            "Failed while waiting for dixon.exe", GetLastError());
            goto cleanup;
        }

        if (!saw_child_output && get_file_size_ull(output_path, &output_size) && output_size > 0) {
            saw_child_output = 1;
        }

        if (!saw_child_output) {
            startup_waited += DIXON_WAIT_SLICE_MS;
            if (startup_waited >= DIXON_STARTUP_TIMEOUT_MS) {
                TerminateProcess(pi.hProcess, 125);
                WaitForSingleObject(pi.hProcess, 5000);
                set_error(error_message, error_message_size,
                          "dixon.exe did not produce any output within 30 seconds. If Windows security is holding the program, run dixon.exe manually once, allow it to start, then retry from the GUI.");
                goto cleanup;
            }
        }
    }

    if (!GetExitCodeProcess(pi.hProcess, &exit_code)) {
        set_win32_error(error_message, error_message_size,
                        "Failed to obtain dixon.exe exit code", GetLastError());
        goto cleanup;
    }
    response->exit_code = (int) exit_code;

    if (!read_file_into_buffer_limited(output_path, &response->stdout_text,
                                       error_message, error_message_size)) {
        goto cleanup;
    }

    response->stderr_text = dup_string("");
    response->combined_output = compose_combined_output(response->command_line,
                                                       response->stdout_text);
    if (!response->stderr_text || !response->combined_output) {
        set_error(error_message, error_message_size, "Out of memory while finalizing the GUI output text.");
        goto cleanup;
    }

    ok = 1;

cleanup:
    if (output_file != INVALID_HANDLE_VALUE) CloseHandle(output_file);
    if (pi.hThread) CloseHandle(pi.hThread);
    if (pi.hProcess) CloseHandle(pi.hProcess);
    if (output_path[0] != '\0') DeleteFileA(output_path);
    free(trimmed_path);
    free(field);
    free(polynomials);
    free(vars);
    free(ideal);
    free(omega);
    free(random_degrees);
    free(actual_command_line);
    sb_free(&command);

    if (!ok) {
        dixon_gui_response_clear(response);
    }

    return ok;
}
