#include <windows.h>
#include <commctrl.h>
#include <commdlg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <mmsystem.h>

#include "dixon_wrapper.h"

#define APP_TITLE "DRsolve Windows GUI"
#define TAB_COUNT 6

#define ID_TAB 100
#define ID_PATH_EDIT 102
#define ID_PATH_BROWSE 103
#define ID_RUN 104
#define ID_CLEAR 105
#define ID_SAVE 106
#define ID_STATUS 107
#define ID_OUTPUT 108
#define ID_OUTPUT_RESULT 109
#define ID_STEP1_METHOD 110
#define ID_STEP4_METHOD 111

#define ID_SOLVE_POLYS 200
#define ID_SOLVE_FIELD 201

#define ID_RES_POLYS 300
#define ID_RES_VARS 301
#define ID_RES_FIELD 302

#define ID_COMP_POLYS 400
#define ID_COMP_VARS 401
#define ID_COMP_FIELD 402
#define ID_COMP_OMEGA 403

#define ID_IDEAL_GENS 500
#define ID_IDEAL_POLYS 501
#define ID_IDEAL_VARS 502
#define ID_IDEAL_FIELD 503

#define ID_RAND_DEGREES 600
#define ID_RAND_FIELD 601
#define ID_RAND_MODE 602
#define ID_RAND_OMEGA 603

#define ID_MUSIC_PLAY 700
#define ID_MUSIC_PAUSE 701
#define ID_MUSIC_STOP 702
#define ID_MUSIC_PROGRESS 703
#define ID_MUSIC_LYRICS 704

#define WM_DIXON_COMPLETE (WM_APP + 1)
#define ID_TAB_TIMER 1
#define ID_MUSIC_TIMER 2

typedef struct {
    HWND panel;
    HWND label1;
    HWND edit1;
    HWND label2;
    HWND edit2;
    HWND label3;
    HWND edit3;
    HWND label4;
    HWND edit4;
} tab_controls_t;

typedef struct {
    dixon_gui_request_t request;
    dixon_gui_response_t response;
    char error[512];
    int ok;
} worker_result_t;

static HINSTANCE g_instance;
static HWND g_main_window;
static HWND g_path_label;
static HWND g_step1_label;
static HWND g_step1_combo;
static HWND g_step4_label;
static HWND g_step4_combo;
static HWND g_tab;
static HWND g_status;
static HWND g_output;
static HWND g_output_result;
static HWND g_path_edit;
static HWND g_run_button;
static HWND g_save_button;
static HWND g_input_group;
static HWND g_output_group;
static HFONT g_ui_font;
static HFONT g_mono_font;
static tab_controls_t g_tabs[TAB_COUNT];
static HANDLE g_worker_thread = NULL;
static int g_is_running = 0;
static WNDPROC g_tab_wndproc = NULL;
static int g_current_tab_index = -1;
static int g_requested_step1_method = -1;
static int g_requested_step4_method = -1;

static HWND g_music_play_button;
static HWND g_music_pause_button;
static HWND g_music_stop_button;
static HWND g_music_progress;
static HWND g_music_lyrics;
static HWND g_music_time_label;
static int g_music_playing = 0;
static WNDPROC g_original_panel_proc = NULL;

// 歌词结构
typedef struct {
    int time;
    char text[256];
} lyric_t;

static lyric_t g_lyrics[100];
static int g_lyric_count = 0;
static int g_current_lyric = -1;

static void show_active_tab(void);
static void init_music(void);
static void play_music(void);
static void pause_music(void);
static void stop_music(void);
static void parse_lyrics(void);
static void update_lyrics(void);
static int get_music_position(void);
static LRESULT CALLBACK panel_wnd_proc(HWND hwnd, UINT msg, WPARAM wparam, LPARAM lparam);

static void populate_method_combo(HWND combo)
{
    SendMessageA(combo, CB_ADDSTRING, 0, (LPARAM) "Default");
    SendMessageA(combo, CB_ADDSTRING, 0, (LPARAM) "0 - Recursive");
    SendMessageA(combo, CB_ADDSTRING, 0, (LPARAM) "1 - Kronecker+HNF");
    SendMessageA(combo, CB_ADDSTRING, 0, (LPARAM) "2 - Interpolation");
    SendMessageA(combo, CB_ADDSTRING, 0, (LPARAM) "3 - Sparse interpolation");
    SendMessageA(combo, CB_SETCURSEL, 0, 0);
}

static int get_method_combo_value(HWND combo)
{
    int index = (int) SendMessageA(combo, CB_GETCURSEL, 0, 0);
    if (index == CB_ERR || index <= 0) return -1;
    return index - 1;
}

#define DRSOLVE_PATH_CAP 4096
#define DRSOLVE_OUTPUT_LIMIT (4 * 1024 * 1024)
#define DRSOLVE_STARTUP_TIMEOUT_MS 30000
#define DRSOLVE_WAIT_SLICE_MS 250

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

    if (!sb_append(&sb, "./drsolve.exe")) {
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
    char path[DRSOLVE_PATH_CAP];
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

    if (snprintf(buffer, buffer_size, "%sdrsolve.exe", path) >= (int) buffer_size) {
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
                  "Cannot find drsolve.exe. Put it next to the GUI or choose it manually.");
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

    if (GetTempFileNameA(temp_dir, "drs", 0, temp_file) == 0) {
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

        if (total + bytes_read > DRSOLVE_OUTPUT_LIMIT) {
            bytes_read = DRSOLVE_OUTPUT_LIMIT - total;
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
        if (total >= DRSOLVE_OUTPUT_LIMIT) truncated = 1;
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
        snprintf(note, sizeof(note), "\n\n[GUI output truncated at %d bytes]\n", DRSOLVE_OUTPUT_LIMIT);
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
    char default_path[DRSOLVE_PATH_CAP];
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
    char output_path[DRSOLVE_PATH_CAP] = {0};
    char working_dir[DRSOLVE_PATH_CAP];
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
                      "Failed to resolve the default drsolve.exe path.");
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
                  "Out of memory while appending the drsolve.exe path to the command line.");
        goto cleanup;
    }

    if (g_requested_step1_method >= 0) {
        char step1_value[16];
        snprintf(step1_value, sizeof(step1_value), "%d", g_requested_step1_method);
        if (!append_arg(&command, "--step1") || !append_arg(&command, step1_value)) {
            set_error(error_message, error_message_size,
                      "Out of memory while appending the step 1 method to the command line.");
            goto cleanup;
        }
    }

    if (g_requested_step4_method >= 0) {
        char step4_value[16];
        snprintf(step4_value, sizeof(step4_value), "%d", g_requested_step4_method);
        if (!append_arg(&command, "--step4") || !append_arg(&command, step4_value)) {
            set_error(error_message, error_message_size,
                      "Out of memory while appending the step 4 method to the command line.");
            goto cleanup;
        }
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
        set_error(error_message, error_message_size, "The drsolve.exe path is too long.");
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
                        "Failed to start drsolve.exe", GetLastError());
        goto cleanup;
    }

    CloseHandle(output_file);
    output_file = INVALID_HANDLE_VALUE;

    for (;;) {
        wait_status = WaitForSingleObject(pi.hProcess, DRSOLVE_WAIT_SLICE_MS);
        if (wait_status == WAIT_OBJECT_0) break;
        if (wait_status != WAIT_TIMEOUT) {
            set_win32_error(error_message, error_message_size,
                            "Failed while waiting for drsolve.exe", GetLastError());
            goto cleanup;
        }

        if (!saw_child_output && get_file_size_ull(output_path, &output_size) && output_size > 0) {
            saw_child_output = 1;
        }

        if (!saw_child_output) {
            startup_waited += DRSOLVE_WAIT_SLICE_MS;
            if (startup_waited >= DRSOLVE_STARTUP_TIMEOUT_MS) {
                TerminateProcess(pi.hProcess, 125);
                WaitForSingleObject(pi.hProcess, 5000);
                set_error(error_message, error_message_size,
                          "drsolve.exe did not produce any output within 30 seconds. If Windows security is holding the program, run drsolve.exe manually once, allow it to start, then retry from the GUI.");
                goto cleanup;
            }
        }
    }

    if (!GetExitCodeProcess(pi.hProcess, &exit_code)) {
        set_win32_error(error_message, error_message_size,
                        "Failed to obtain drsolve.exe exit code", GetLastError());
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


static int max_int(int a, int b)
{
    return a > b ? a : b;
}

static void set_status(const char *text)
{
    SendMessageA(g_status, SB_SETTEXTA, 0, (LPARAM) (text ? text : ""));
}

static char *normalize_edit_text_newlines(const char *text)
{
    const char *src;
    char *buffer;
    char *dst;
    size_t extra = 0;
    size_t len;

    if (!text) return NULL;

    for (src = text; *src; ++src) {
        if (*src == '\n' && (src == text || src[-1] != '\r')) extra++;
    }

    len = strlen(text);
    buffer = (char *) malloc(len + extra + 1);
    if (!buffer) return NULL;

    dst = buffer;
    for (src = text; *src; ++src) {
        if (*src == '\n' && (src == text || src[-1] != '\r')) *dst++ = '\r';
        *dst++ = *src;
    }
    *dst = '\0';
    return buffer;
}

static void set_edit_text(HWND edit, const char *text)
{
    char *normalized;

    if (!text) {
        SetWindowTextA(edit, "");
        return;
    }

    normalized = normalize_edit_text_newlines(text);
    if (!normalized) {
        SetWindowTextA(edit, text);
        return;
    }

    SetWindowTextA(edit, normalized);
    free(normalized);
}

static int is_space_char(char ch)
{
    return ch == ' ' || ch == '\t' || ch == '\r' || ch == '\n';
}

static void trim_span(const char **start, const char **end)
{
    while (*start < *end && is_space_char(**start)) (*start)++;
    while (*end > *start && is_space_char((*end)[-1])) (*end)--;
}

static int extract_saved_path(const char *text, char *out, size_t out_size)
{
    const char *p = text;
    int found = 0;

    if (!out || out_size == 0) return 0;
    out[0] = '\0';
    if (!text) return 0;

    while (p && *p) {
        const char *line_start = p;
        const char *line_end = strchr(line_start, '\n');
        const char *line_trim_start = line_start;
        const char *line_trim_end = line_end ? line_end : line_start + strlen(line_start);
        const char *saved_to;
        const char *colon;
        const char *value_start;
        const char *value_end;
        size_t len;

        trim_span(&line_trim_start, &line_trim_end);
        saved_to = strstr(line_trim_start, "saved to");
        if (!saved_to || saved_to >= line_trim_end) {
            if (!line_end) break;
            p = line_end + 1;
            continue;
        }

        colon = strchr(saved_to, ':');
        if (!colon || colon >= line_trim_end) {
            if (!line_end) break;
            p = line_end + 1;
            continue;
        }

        value_start = colon + 1;
        value_end = line_trim_end;
        trim_span(&value_start, &value_end);
        len = (size_t) (value_end - value_start);
        if (len > 0) {
            if (len >= out_size) len = out_size - 1;
            memcpy(out, value_start, len);
            out[len] = '\0';
            found = 1;
        }

        if (!line_end) break;
        p = line_end + 1;
    }

    return found;
}

static int is_absolute_path(const char *path)
{
    if (!path || !path[0]) return 0;
    if ((path[0] == '\\' && path[1] == '\\') || path[0] == '/') return 1;
    return isalpha((unsigned char) path[0]) && path[1] == ':';
}

static int join_result_path(const char *dixon_path, const char *saved_path, char *out, size_t out_size)
{
    char base[MAX_PATH];
    char *last_slash;

    if (!saved_path || !saved_path[0] || !out || out_size == 0) return 0;
    if (is_absolute_path(saved_path)) {
        if (snprintf(out, out_size, "%s", saved_path) >= (int) out_size) return 0;
        return 1;
    }

    if (!dixon_path || !dixon_path[0]) return 0;
    if (snprintf(base, sizeof(base), "%s", dixon_path) >= (int) sizeof(base)) return 0;
    last_slash = strrchr(base, '\\');
    if (!last_slash) last_slash = strrchr(base, '/');
    if (last_slash) {
        *last_slash = '\0';
    } else {
        snprintf(base, sizeof(base), ".");
    }

    if (snprintf(out, out_size, "%s\\%s", base, saved_path) >= (int) out_size) return 0;
    return 1;
}

static char *read_text_file_limited(const char *path, size_t max_bytes)
{
    FILE *fp;
    long file_size;
    size_t read_size;
    char *buffer;
    size_t got;

    if (!path || !path[0]) return NULL;
    fp = fopen(path, "rb");
    if (!fp) return NULL;
    if (fseek(fp, 0, SEEK_END) != 0) {
        fclose(fp);
        return NULL;
    }
    file_size = ftell(fp);
    if (file_size < 0) {
        fclose(fp);
        return NULL;
    }
    if (fseek(fp, 0, SEEK_SET) != 0) {
        fclose(fp);
        return NULL;
    }

    read_size = (size_t) file_size;
    if (read_size > max_bytes) read_size = max_bytes;
    buffer = (char *) malloc(read_size + 1);
    if (!buffer) {
        fclose(fp);
        return NULL;
    }

    got = fread(buffer, 1, read_size, fp);
    buffer[got] = '\0';
    fclose(fp);
    return buffer;
}

static void update_result_output(worker_result_t *result)
{
    char saved_path[MAX_PATH];
    char full_path[MAX_PATH];
    char *content;

    if (!result || !result->ok) {
        set_edit_text(g_output_result, "");
        return;
    }

    if (!extract_saved_path(result->response.stdout_text, saved_path, sizeof(saved_path))) {
        set_edit_text(g_output_result, "No result/report file path found in output.");
        return;
    }

    if (!join_result_path(result->response.resolved_dixon_path ? result->response.resolved_dixon_path : result->request.dixon_path,
                          saved_path, full_path, sizeof(full_path))) {
        set_edit_text(g_output_result, "Failed to resolve result/report file path.");
        return;
    }

    content = read_text_file_limited(full_path, 2 * 1024 * 1024);
    if (!content) {
        char message[512];
        snprintf(message, sizeof(message),
                 "Detected result/report file, but failed to read:\n%s", full_path);
        set_edit_text(g_output_result, message);
        return;
    }

    set_edit_text(g_output_result, content);
    free(content);
}

static char *get_edit_text(HWND edit)
{
    int len = GetWindowTextLengthA(edit);
    char *buffer = (char *) malloc((size_t) len + 1);
    if (!buffer) return NULL;
    GetWindowTextA(edit, buffer, len + 1);
    return buffer;
}

static void apply_font(HWND hwnd, HFONT font)
{
    SendMessage(hwnd, WM_SETFONT, (WPARAM) font, TRUE);
}

static HWND create_label(HWND parent, const char *text, int x, int y, int w, int h)
{
    HWND hwnd = CreateWindowExA(
        0, "STATIC", text,
        WS_CHILD | WS_VISIBLE,
        x, y, w, h,
        parent, NULL, g_instance, NULL);
    apply_font(hwnd, g_ui_font);
    return hwnd;
}

static HWND create_groupbox(HWND parent, const char *text)
{
    HWND hwnd = CreateWindowExA(
        0, "BUTTON", text,
        WS_CHILD | WS_VISIBLE | BS_GROUPBOX,
        0, 0, 10, 10,
        parent, NULL, g_instance, NULL);
    apply_font(hwnd, g_ui_font);
    return hwnd;
}

static HWND create_edit(HWND parent, int id, DWORD style, int x, int y, int w, int h, const char *text)
{
    HWND hwnd = CreateWindowExA(
        WS_EX_CLIENTEDGE, "EDIT", text ? text : "",
        WS_CHILD | WS_VISIBLE | style,
        x, y, w, h,
        parent, (HMENU) (INT_PTR) id, g_instance, NULL);
    apply_font(hwnd, (style & ES_MULTILINE) ? g_mono_font : g_ui_font);
    return hwnd;
}

static HWND create_button(HWND parent, int id, const char *text, int x, int y, int w, int h)
{
    HWND hwnd = CreateWindowExA(
        0, "BUTTON", text,
        WS_CHILD | WS_VISIBLE | BS_PUSHBUTTON,
        x, y, w, h,
        parent, (HMENU) (INT_PTR) id, g_instance, NULL);
    apply_font(hwnd, g_ui_font);
    return hwnd;
}

static HWND create_combobox(HWND parent, int id, int x, int y, int w, int h)
{
    HWND hwnd = CreateWindowExA(
        0, "COMBOBOX", "",
        WS_CHILD | WS_VISIBLE | WS_TABSTOP | CBS_DROPDOWNLIST | WS_VSCROLL,
        x, y, w, h,
        parent, (HMENU) (INT_PTR) id, g_instance, NULL);
    apply_font(hwnd, g_ui_font);
    return hwnd;
}

static void update_random_mode_controls(void)
{
    int index = (int) SendMessageA(g_tabs[4].edit3, CB_GETCURSEL, 0, 0);
    EnableWindow(g_tabs[4].label4, index == DIXON_GUI_RANDOM_COMPLEXITY);
    EnableWindow(g_tabs[4].edit4, index == DIXON_GUI_RANDOM_COMPLEXITY);
}

static LRESULT CALLBACK tab_wnd_proc(HWND hwnd, UINT msg, WPARAM wparam, LPARAM lparam)
{
    LRESULT result;

    result = CallWindowProcA(g_tab_wndproc, hwnd, msg, wparam, lparam);

    switch (msg) {
        case WM_LBUTTONUP:
        case WM_KEYUP:
        case WM_MOUSEWHEEL:
            show_active_tab();
            break;
        default:
            break;
    }

    return result;
}

static void show_active_tab(void)
{
    int current = TabCtrl_GetCurSel(g_tab);
    int i;

    if (current < 0 || current >= TAB_COUNT) current = 0;
    g_current_tab_index = current;

    for (i = 0; i < TAB_COUNT; ++i) {
        ShowWindow(g_tabs[i].panel, i == current ? SW_SHOW : SW_HIDE);
    }

    if (current == 4) update_random_mode_controls();
}

static void layout_tab_pages(void)
{
    RECT rc;
    int page_w;
    int page_h;
    int margin = 12;
    int i;

    GetClientRect(g_tab, &rc);
    TabCtrl_AdjustRect(g_tab, FALSE, &rc);

    page_w = max_int(rc.right - rc.left, 240);
    page_h = max_int(rc.bottom - rc.top, 240);

    for (i = 0; i < TAB_COUNT; ++i) {
        MoveWindow(g_tabs[i].panel, rc.left, rc.top, page_w, page_h, TRUE);
    }

    {
        int label_y = 12;
        int poly_h = max_int(page_h - 118, 130);
        int row_y = 34 + poly_h + 14;
        int field_w = 220;

        MoveWindow(g_tabs[0].label1, margin, label_y, 120, 20, TRUE);
        MoveWindow(g_tabs[0].edit1, margin, 34, page_w - margin * 2, poly_h, TRUE);
        MoveWindow(g_tabs[0].label2, margin, row_y, 90, 20, TRUE);
        MoveWindow(g_tabs[0].edit2, margin, row_y + 22, field_w, 26, TRUE);
    }

    {
        int label_y = 12;
        int poly_h = max_int(page_h - 150, 104);
        int row1_y = 34 + poly_h + 12;
        int row2_y = row1_y + 58;

        MoveWindow(g_tabs[1].label1, margin, label_y, 120, 20, TRUE);
        MoveWindow(g_tabs[1].edit1, margin, 34, page_w - margin * 2, poly_h, TRUE);
        MoveWindow(g_tabs[1].label2, margin, row1_y, 140, 20, TRUE);
        MoveWindow(g_tabs[1].edit2, margin, row1_y + 22, page_w - margin * 2, 26, TRUE);
        MoveWindow(g_tabs[1].label3, margin, row2_y, 80, 20, TRUE);
        MoveWindow(g_tabs[1].edit3, margin, row2_y + 22, 220, 26, TRUE);
    }

    {
        int label_y = 12;
        int poly_h = max_int(page_h - 150, 96);
        int row1_y = 34 + poly_h + 12;
        int row2_y = row1_y + 58;

        MoveWindow(g_tabs[2].label1, margin, label_y, 120, 20, TRUE);
        MoveWindow(g_tabs[2].edit1, margin, 34, page_w - margin * 2, poly_h, TRUE);
        MoveWindow(g_tabs[2].label2, margin, row1_y, 140, 20, TRUE);
        MoveWindow(g_tabs[2].edit2, margin, row1_y + 22, page_w - margin * 2, 26, TRUE);
        MoveWindow(g_tabs[2].label3, margin, row2_y, 80, 20, TRUE);
        MoveWindow(g_tabs[2].edit3, margin, row2_y + 22, 220, 26, TRUE);
        MoveWindow(g_tabs[2].label4, 252, row2_y, 80, 20, TRUE);
        MoveWindow(g_tabs[2].edit4, 252, row2_y + 22, 140, 26, TRUE);
    }

    {
        int label_y = 12;
        int usable_h = page_h - 152;
        int gens_h = usable_h / 2;
        int polys_h;
        int row_y;

        if (gens_h < 72) gens_h = 72;
        polys_h = usable_h - gens_h;
        if (polys_h < 72) {
            polys_h = 72;
            gens_h = max_int(usable_h - polys_h, 72);
        }

        MoveWindow(g_tabs[3].label1, margin, label_y, 140, 20, TRUE);
        MoveWindow(g_tabs[3].edit1, margin, 34, page_w - margin * 2, gens_h, TRUE);
        MoveWindow(g_tabs[3].label2, margin, 46 + gens_h, 120, 20, TRUE);
        MoveWindow(g_tabs[3].edit2, margin, 68 + gens_h, page_w - margin * 2, polys_h, TRUE);

        row_y = 80 + gens_h + polys_h;
        MoveWindow(g_tabs[3].label3, margin, row_y, 140, 20, TRUE);
        MoveWindow(g_tabs[3].edit3, margin, row_y + 22, 220, 26, TRUE);
        MoveWindow(g_tabs[3].label4, 252, row_y, 80, 20, TRUE);
        MoveWindow(g_tabs[3].edit4, 252, row_y + 22, 180, 26, TRUE);
    }

    {
        int label_y = 18;
        int row1_y = 48;
        int row2_y = row1_y + 60;
        int row3_y = row2_y + 60;

        MoveWindow(g_tabs[4].label1, margin, label_y, 120, 20, TRUE);
        MoveWindow(g_tabs[4].edit1, margin, row1_y, page_w - margin * 2, 28, TRUE);
        MoveWindow(g_tabs[4].label2, margin, row2_y, 80, 20, TRUE);
        MoveWindow(g_tabs[4].edit2, margin, row2_y + 22, 220, 26, TRUE);
        MoveWindow(g_tabs[4].label3, 264, row2_y, 80, 20, TRUE);
        MoveWindow(g_tabs[4].edit3, 264, row2_y + 22, 220, 220, TRUE);
        MoveWindow(g_tabs[4].label4, margin, row3_y, 80, 20, TRUE);
        MoveWindow(g_tabs[4].edit4, margin, row3_y + 22, 140, 26, TRUE);
    }

    { // Music tab layout
        int button_y = 20;
        int button_w = 80;
        int button_h = 30;
        int button_gap = 10;
        int progress_y = button_y + button_h + 20;
        int lyrics_y = progress_y + 40;
        int lyrics_h = page_h - lyrics_y - 20;

        MoveWindow(g_music_play_button, margin, button_y, button_w, button_h, TRUE);
        MoveWindow(g_music_pause_button, margin + button_w + button_gap, button_y, button_w, button_h, TRUE);
        MoveWindow(g_music_stop_button, margin + (button_w + button_gap) * 2, button_y, button_w, button_h, TRUE);
        MoveWindow(g_music_progress, margin, progress_y, page_w - margin * 2 - 100, 25, TRUE);
        MoveWindow(g_music_time_label, page_w - margin - 100, progress_y, 100, 25, TRUE);
        MoveWindow(g_music_lyrics, margin, lyrics_y, page_w - margin * 2, lyrics_h, TRUE);
        
        // Force redraw of all music controls
        RedrawWindow(g_music_progress, NULL, NULL, RDW_FRAME | RDW_INVALIDATE | RDW_UPDATENOW);
        RedrawWindow(g_music_time_label, NULL, NULL, RDW_FRAME | RDW_INVALIDATE | RDW_UPDATENOW);
        RedrawWindow(g_music_lyrics, NULL, NULL, RDW_FRAME | RDW_INVALIDATE | RDW_UPDATENOW);
    }
}

static void create_tab_pages(void)
{
    TCITEMA item;

    item.mask = TCIF_TEXT;

    item.pszText = "Solve";
    TabCtrl_InsertItem(g_tab, 0, &item);
    item.pszText = "Dixon";
    TabCtrl_InsertItem(g_tab, 1, &item);
    item.pszText = "Complexity";
    TabCtrl_InsertItem(g_tab, 2, &item);
    item.pszText = "Ideal";
    TabCtrl_InsertItem(g_tab, 3, &item);
    item.pszText = "Random";
    TabCtrl_InsertItem(g_tab, 4, &item);
    item.pszText = "Music";
    TabCtrl_InsertItem(g_tab, 5, &item);

    g_tabs[0].panel = CreateWindowExA(0, "STATIC", "", WS_CHILD | WS_VISIBLE,
                                      0, 0, 10, 10, g_tab, NULL, g_instance, NULL);
    g_tabs[0].label1 = create_label(g_tabs[0].panel, "Polynomials", 0, 0, 10, 10);
    g_tabs[0].edit1 = create_edit(g_tabs[0].panel, ID_SOLVE_POLYS,
                                  ES_LEFT | ES_MULTILINE | ES_AUTOVSCROLL | WS_VSCROLL,
                                  0, 0, 10, 10,
                                  "x^2+y^2+z^2-6\r\nx+y+z-4\r\nx*y*z-x-1");
    g_tabs[0].label2 = create_label(g_tabs[0].panel, "Field", 0, 0, 10, 10);
    g_tabs[0].edit2 = create_edit(g_tabs[0].panel, ID_SOLVE_FIELD, ES_LEFT | WS_TABSTOP,
                                  0, 0, 10, 10, "257");

    g_tabs[1].panel = CreateWindowExA(0, "STATIC", "", WS_CHILD | WS_VISIBLE,
                                      0, 0, 10, 10, g_tab, NULL, g_instance, NULL);
    g_tabs[1].label1 = create_label(g_tabs[1].panel, "Polynomials", 0, 0, 10, 10);
    g_tabs[1].edit1 = create_edit(g_tabs[1].panel, ID_RES_POLYS,
                                  ES_LEFT | ES_MULTILINE | ES_AUTOVSCROLL | WS_VSCROLL,
                                  0, 0, 10, 10,
                                  "x+y+z\r\nx*y+y*z+z*x\r\nx*y*z+1");
    g_tabs[1].label2 = create_label(g_tabs[1].panel, "Eliminate Variables", 0, 0, 10, 10);
    g_tabs[1].edit2 = create_edit(g_tabs[1].panel, ID_RES_VARS, ES_LEFT | WS_TABSTOP,
                                  0, 0, 10, 10, "x,y");
    g_tabs[1].label3 = create_label(g_tabs[1].panel, "Field", 0, 0, 10, 10);
    g_tabs[1].edit3 = create_edit(g_tabs[1].panel, ID_RES_FIELD, ES_LEFT | WS_TABSTOP,
                                  0, 0, 10, 10, "257");

    g_tabs[2].panel = CreateWindowExA(0, "STATIC", "", WS_CHILD | WS_VISIBLE,
                                      0, 0, 10, 10, g_tab, NULL, g_instance, NULL);
    g_tabs[2].label1 = create_label(g_tabs[2].panel, "Polynomials", 0, 0, 10, 10);
    g_tabs[2].edit1 = create_edit(g_tabs[2].panel, ID_COMP_POLYS,
                                  ES_LEFT | ES_MULTILINE | ES_AUTOVSCROLL | WS_VSCROLL,
                                  0, 0, 10, 10,
                                  "x^2+y^2+1\r\nx*y+z\r\nx+y+z^2");
    g_tabs[2].label2 = create_label(g_tabs[2].panel, "Eliminate Variables", 0, 0, 10, 10);
    g_tabs[2].edit2 = create_edit(g_tabs[2].panel, ID_COMP_VARS, ES_LEFT | WS_TABSTOP,
                                  0, 0, 10, 10, "x,y");
    g_tabs[2].label3 = create_label(g_tabs[2].panel, "Field", 0, 0, 10, 10);
    g_tabs[2].edit3 = create_edit(g_tabs[2].panel, ID_COMP_FIELD, ES_LEFT | WS_TABSTOP,
                                  0, 0, 10, 10, "257");
    g_tabs[2].label4 = create_label(g_tabs[2].panel, "Omega", 0, 0, 10, 10);
    g_tabs[2].edit4 = create_edit(g_tabs[2].panel, ID_COMP_OMEGA, ES_LEFT | WS_TABSTOP,
                                  0, 0, 10, 10, "2.373");

    g_tabs[3].panel = CreateWindowExA(0, "STATIC", "", WS_CHILD | WS_VISIBLE,
                                      0, 0, 10, 10, g_tab, NULL, g_instance, NULL);
    g_tabs[3].label1 = create_label(g_tabs[3].panel, "Ideal Generators", 0, 0, 10, 10);
    g_tabs[3].edit1 = create_edit(g_tabs[3].panel, ID_IDEAL_GENS,
                                  ES_LEFT | ES_MULTILINE | ES_AUTOVSCROLL | WS_VSCROLL,
                                  0, 0, 10, 10,
                                  "a2^3=2*a1+1\r\na3^3=a1*a2+3");
    g_tabs[3].label2 = create_label(g_tabs[3].panel, "Polynomials", 0, 0, 10, 10);
    g_tabs[3].edit2 = create_edit(g_tabs[3].panel, ID_IDEAL_POLYS,
                                  ES_LEFT | ES_MULTILINE | ES_AUTOVSCROLL | WS_VSCROLL,
                                  0, 0, 10, 10,
                                  "a1^2+a2^2+a3^2-10\r\na3^3-a1*a2-3");
    g_tabs[3].label3 = create_label(g_tabs[3].panel, "Eliminate Variables", 0, 0, 10, 10);
    g_tabs[3].edit3 = create_edit(g_tabs[3].panel, ID_IDEAL_VARS, ES_LEFT | WS_TABSTOP,
                                  0, 0, 10, 10, "a3");
    g_tabs[3].label4 = create_label(g_tabs[3].panel, "Field", 0, 0, 10, 10);
    g_tabs[3].edit4 = create_edit(g_tabs[3].panel, ID_IDEAL_FIELD, ES_LEFT | WS_TABSTOP,
                                  0, 0, 10, 10, "257");

    g_tabs[4].panel = CreateWindowExA(0, "STATIC", "", WS_CHILD | WS_VISIBLE,
                                      0, 0, 10, 10, g_tab, NULL, g_instance, NULL);
    g_tabs[4].label1 = create_label(g_tabs[4].panel, "Degree List", 0, 0, 10, 10);
    g_tabs[4].edit1 = create_edit(g_tabs[4].panel, ID_RAND_DEGREES, ES_LEFT | WS_TABSTOP,
                                  0, 0, 10, 10, "[3,3,2]");
    g_tabs[4].label2 = create_label(g_tabs[4].panel, "Field", 0, 0, 10, 10);
    g_tabs[4].edit2 = create_edit(g_tabs[4].panel, ID_RAND_FIELD, ES_LEFT | WS_TABSTOP,
                                  0, 0, 10, 10, "257");
    g_tabs[4].label3 = create_label(g_tabs[4].panel, "Mode", 0, 0, 10, 10);
    g_tabs[4].edit3 = create_combobox(g_tabs[4].panel, ID_RAND_MODE, 0, 0, 10, 120);
    SendMessageA(g_tabs[4].edit3, CB_ADDSTRING, 0, (LPARAM) "Resultant");
    SendMessageA(g_tabs[4].edit3, CB_ADDSTRING, 0, (LPARAM) "Solve");
    SendMessageA(g_tabs[4].edit3, CB_ADDSTRING, 0, (LPARAM) "Complexity");
    SendMessageA(g_tabs[4].edit3, CB_SETCURSEL, DIXON_GUI_RANDOM_RESULTANT, 0);
    g_tabs[4].label4 = create_label(g_tabs[4].panel, "Omega", 0, 0, 10, 10);
    g_tabs[4].edit4 = create_edit(g_tabs[4].panel, ID_RAND_OMEGA, ES_LEFT | WS_TABSTOP,
                                  0, 0, 10, 10, "2.373");

    g_tabs[5].panel = CreateWindowExA(0, "STATIC", "", WS_CHILD | WS_VISIBLE,
                                      0, 0, 10, 10, g_tab, NULL, g_instance, NULL);
    // Save original window procedure and set our custom procedure
    g_original_panel_proc = (WNDPROC) SetWindowLongPtrA(g_tabs[5].panel, GWLP_WNDPROC, (LONG_PTR) panel_wnd_proc);
    
    g_music_play_button = create_button(g_tabs[5].panel, ID_MUSIC_PLAY, "Play", 0, 0, 10, 10);
    g_music_pause_button = create_button(g_tabs[5].panel, ID_MUSIC_PAUSE, "Pause", 0, 0, 10, 10);
    g_music_stop_button = create_button(g_tabs[5].panel, ID_MUSIC_STOP, "Stop", 0, 0, 10, 10);
    
    // Create progress bar with PBS_SMOOTH style for better visualization
    g_music_progress = CreateWindowExA(0, "msctls_progress32", "",
                                      WS_CHILD | WS_VISIBLE | PBS_SMOOTH,
                                      0, 0, 10, 10,
                                      g_tabs[5].panel, (HMENU) (INT_PTR) ID_MUSIC_PROGRESS, g_instance, NULL);
    // Initialize progress bar range and position
    SendMessageA(g_music_progress, PBM_SETRANGE, 0, MAKELPARAM(0, 100));
    SendMessageA(g_music_progress, PBM_SETPOS, 0, 0);
    // Set progress bar color (blue)
    SendMessageA(g_music_progress, PBM_SETBARCOLOR, 0, RGB(0, 122, 204));
    // Force initial redraw
    RedrawWindow(g_music_progress, NULL, NULL, RDW_FRAME | RDW_INVALIDATE | RDW_UPDATENOW);
    
    // 创建时间显示标签
    g_music_time_label = create_label(g_tabs[5].panel, "00:00 / 00:00", 0, 0, 10, 10);
    
    // Create lyrics display area
    g_music_lyrics = create_edit(g_tabs[5].panel, ID_MUSIC_LYRICS,
                                ES_LEFT | ES_MULTILINE | ES_AUTOVSCROLL | ES_READONLY | WS_VSCROLL | WS_TABSTOP,
                                0, 0, 10, 10,
                                "Click Play to start music...");
    
    // Show initial lyrics (first three lines) immediately when tab is created
    if (g_lyric_count > 0) {
        char lyric_text[1024] = "";
        
        // Show first three lines
        if (g_lyric_count >= 3) {
            strcat(lyric_text, g_lyrics[0].text);
            strcat(lyric_text, "\r\n");
            strcat(lyric_text, g_lyrics[1].text);
            strcat(lyric_text, "\r\n");
            strcat(lyric_text, g_lyrics[2].text);
        } else if (g_lyric_count == 2) {
            strcat(lyric_text, g_lyrics[0].text);
            strcat(lyric_text, "\r\n");
            strcat(lyric_text, g_lyrics[1].text);
        } else if (g_lyric_count == 1) {
            strcat(lyric_text, g_lyrics[0].text);
        }
        
        set_edit_text(g_music_lyrics, lyric_text);
    }

    // 初始化音乐
    init_music();

    layout_tab_pages();
    show_active_tab();
}

static void layout_main_window(HWND hwnd)
{
    RECT rc;
    int width;
    int height;
    int margin = 12;
    int path_y = 18;
    int method_y = 58;
    int row_h = 28;
    int browse_w = 112;
    int label_w = 62;
    int combo_w = 220;
    int input_h = 390;
    int button_y;
    int output_y;
    int output_h;
    int status_h = 24;
    int content_width;

    GetClientRect(hwnd, &rc);
    width = rc.right - rc.left;
    height = rc.bottom - rc.top;
    content_width = width - margin * 2;

    if (height < 760) input_h = 350;
    if (height < 680) input_h = 310;
    if (height > 920) input_h = 430;

    if (input_h > height - 286) input_h = height - 286;
    if (input_h < 230) input_h = 230;

    MoveWindow(g_path_label, margin, path_y + 4, label_w, 20, TRUE);
    MoveWindow(g_path_edit, margin + label_w + 8, path_y,
               content_width - label_w - browse_w - 24,
               row_h, TRUE);
    MoveWindow(GetDlgItem(hwnd, ID_PATH_BROWSE), width - margin - browse_w, path_y - 1,
               browse_w, row_h + 2, TRUE);

    MoveWindow(g_step1_label, margin, method_y + 4, 90, 20, TRUE);
    MoveWindow(g_step1_combo, margin + 94, method_y, combo_w, 240, TRUE);
    MoveWindow(g_step4_label, margin + 94 + combo_w + 18, method_y + 4, 90, 20, TRUE);
    MoveWindow(g_step4_combo, margin + 94 + combo_w + 18 + 94, method_y, combo_w, 240, TRUE);

    MoveWindow(g_input_group, margin, 98, content_width, input_h, TRUE);
    MoveWindow(g_tab, 12, 30, content_width - 24, input_h - 42, TRUE);
    layout_tab_pages();

    button_y = 98 + input_h + 14;
    MoveWindow(g_run_button, margin, button_y, 112, 30, TRUE);
    MoveWindow(GetDlgItem(hwnd, ID_CLEAR), margin + 124, button_y, 112, 30, TRUE);
    MoveWindow(g_save_button, margin + 248, button_y, 120, 30, TRUE);

    output_y = button_y + 42;
    output_h = height - output_y - status_h - margin;
    if (output_h < 140) output_h = 140;

    MoveWindow(g_output_group, margin, output_y, content_width, output_h, TRUE);
    {
        int inner_w = content_width - 24;
        int inner_h = output_h - 34;
        int split_gap = 10;
        int left_w = (inner_w - split_gap) / 2;
        int right_w = inner_w - split_gap - left_w;
        if (left_w < 120) left_w = 120;
        if (right_w < 120) right_w = 120;
        MoveWindow(g_output, 12, 22, left_w, inner_h, TRUE);
        MoveWindow(g_output_result, 12 + left_w + split_gap, 22, right_w, inner_h, TRUE);
    }

    SendMessage(g_status, WM_SIZE, 0, 0);
}

static void browse_for_dixon(void)
{
    OPENFILENAMEA ofn;
    char path[MAX_PATH] = "";

    ZeroMemory(&ofn, sizeof(ofn));
    ofn.lStructSize = sizeof(ofn);
    ofn.hwndOwner = g_main_window;
    ofn.lpstrFilter = "Executables\0*.exe\0All Files\0*.*\0";
    ofn.lpstrFile = path;
    ofn.nMaxFile = sizeof(path);
    ofn.Flags = OFN_FILEMUSTEXIST | OFN_PATHMUSTEXIST;

    if (GetOpenFileNameA(&ofn)) {
        set_edit_text(g_path_edit, path);
        set_status("Using selected drsolve.exe");
    }
}

static void save_output(void)
{
    OPENFILENAMEA ofn;
    char path[MAX_PATH] = "drsolve_win_output.txt";
    char *left_text;
    char *right_text;
    char *text;
    size_t left_len;
    size_t right_len;
    FILE *fp;

    left_text = get_edit_text(g_output);
    right_text = get_edit_text(g_output_result);
    left_len = left_text ? strlen(left_text) : 0;
    right_len = right_text ? strlen(right_text) : 0;
    if (left_len == 0 && right_len == 0) {
        free(left_text);
        free(right_text);
        MessageBoxA(g_main_window, "There is no output to save.", APP_TITLE, MB_OK | MB_ICONINFORMATION);
        return;
    }

    if (left_len > 0 && right_len > 0) {
        size_t sep_len = 4;
        text = (char *) malloc(left_len + sep_len + right_len + 1);
        if (!text) {
            free(left_text);
            free(right_text);
            MessageBoxA(g_main_window, "Out of memory while preparing output text.", APP_TITLE, MB_OK | MB_ICONERROR);
            return;
        }
        memcpy(text, left_text, left_len);
        memcpy(text + left_len, "\r\n\r\n", sep_len);
        memcpy(text + left_len + sep_len, right_text, right_len);
        text[left_len + sep_len + right_len] = '\0';
    } else if (left_len > 0) {
        text = left_text;
        left_text = NULL;
    } else {
        text = right_text;
        right_text = NULL;
    }

    ZeroMemory(&ofn, sizeof(ofn));
    ofn.lStructSize = sizeof(ofn);
    ofn.hwndOwner = g_main_window;
    ofn.lpstrFilter = "Text Files\0*.txt\0All Files\0*.*\0";
    ofn.lpstrFile = path;
    ofn.nMaxFile = sizeof(path);
    ofn.Flags = OFN_OVERWRITEPROMPT | OFN_PATHMUSTEXIST;

    if (GetSaveFileNameA(&ofn)) {
        fp = fopen(path, "wb");
        if (!fp) {
            MessageBoxA(g_main_window, "Failed to save the output file.", APP_TITLE, MB_OK | MB_ICONERROR);
        } else {
            fwrite(text, 1, strlen(text), fp);
            fclose(fp);
            set_status("Output saved");
        }
    }

    free(text);
    free(left_text);
    free(right_text);
}

static void init_music(void)
{
    // Initialize music-related controls
    parse_lyrics();
}

static void parse_lyrics(void)
{
    // Parse lyrics with timestamps
    g_lyric_count = 0;
    
    // Lyric data
    static char lyric_data[] = 
        "[00:00.001]Dixon Song\n"
        "[00:10.020](Verse1)\n"
        "[00:10.200]When the system's full of equations, hard to crack,\n"
        "[00:16.140]Groebner bases are powerful, yet slow to come back.\n"
        "[00:20.280]Zero-dimensional is all they know, that's such a shame,\n"
        "[00:26.790]FGLM can't help you when the unknowns aren't the same.\n"
        "[00:29.940](Chorus)\n"
        "[00:31.200]Oh, Dixon, Dixon, where'd you go?\n"
        "[00:33.600]Dixon, Dixon, steal the show.\n"
        "[00:36.060]Oh, Dixon, Dixon, where'd you go?\n"
        "[00:38.820]Dixon, Dixon, steal the show.\n"
        "[00:41.640](Verse2)\n"
        "[00:47.460]Sylvester's old and slow, one variable at a time,\n"
        "[00:53.040]For multivariate systems, that's a weary climb.\n"
        "[00:56.610]But Dixon builds one matrix, knocks out several at once,\n"
        "[01:01.890]Eliminates in one go, no more this clunky grind.\n"
        "[01:06.270](Chorus)\n"
        "[01:07.290]Oh, Dixon, Dixon, where'd you go?\n"
        "[01:09.960]Dixon, Dixon, steal the show.\n"
        "[01:12.420]Oh, Dixon, Dixon, where'd you go?\n"
        "[01:15.180]Dixon, Dixon, steal the show.\n"
        "[01:17.850](Verse3)\n"
        "[01:22.860]We count the lattice paths, we bound the matrix size,\n"
        "[01:27.780]Complexity gets tighter now, no more surprise.\n"
        "[01:32.940]Over-determined, under-determined, both we face,\n"
        "[01:38.550]Dixon fills the gap and finds its place.\n"
        "[01:42.090](Chorus)\n"
        "[01:43.620]Oh, Dixon, Dixon, where'd you go?\n"
        "[01:46.260]Dixon, Dixon, steal the show.\n"
        "[01:48.780]Oh, Dixon, Dixon, where'd you go?\n"
        "[01:51.540]Dixon, Dixon, steal the show.\n"
        "[01:53.670](Verse4)\n"
        "[01:59.100]We coded it in C, open-source and fast,\n"
        "[02:05.460]Submatrix selection makes the extraneous past.\n"
        "[02:08.970]Magma and msolve, they're strong, they're great,\n"
        "[02:14.520]But Dixon holds its own, don't underestimate.\n"
        "[02:18.810](Chorus)\n"
        "[02:20.760]Oh, Dixon, Dixon, where'd you go?\n"
        "[02:22.950]Dixon, Dixon, steal the show.\n"
        "[02:25.230]Oh, Dixon, Dixon, where'd you go?\n"
        "[02:27.240]Dixon, Dixon, steal the show.\n"
        "[02:29.970](Verse5)\n"
        "[02:37.440]From AO primitives to algebraic insight,\n"
        "[02:40.800]Algebraic attacks are coming, day and night.\n"
        "[02:45.030]We ran the Dixon, we counted every round,\n"
        "[02:50.940]Our estimates are lower, new bounds we've found.\n"
        "[02:55.020](Chorus)\n"
        "[02:56.820]Oh, Dixon, Dixon, where'd you go?\n"
        "[02:58.920]Dixon, Dixon, steal the show.\n"
        "[03:01.200]Oh, Dixon, Dixon, where'd you go?\n"
        "[03:04.140]Dixon, Dixon, steal~\n"
        "[03:07.170]Oh, Dixon, Dixon, now we know,\n"
        "[03:09.300]Dixon, Dixon, steals the show.\n"
        "[03:11.790]Oh, Dixon, Dixon, now we know,\n"
        "[03:14.580]Dixon, Dixon, steal~\n"
        "[03:17.220]Oh, Dixon, Dixon, where'd you go?\n"
        "[03:19.680]Dixon, Dixon, steals the show.\n"
        "[03:22.110]Oh, Dixon, Dixon, where'd you go?\n"
        "[03:24.930]Dixon, Dixon, steals the show.\n"
        "[03:30.000](Finale)\n"
        "[03:32.370]Oh, Dixon, Dixon, now we know,\n"
        "[03:35.130]Dixon, Dixon, steals the show.\n"
        "[03:37.680]Oh, Dixon, Dixon, now we know,\n"
        "[03:40.500]Dixon, Dixon, steals the show.";
    
    char* line = strtok(lyric_data, "\n");
    while (line && g_lyric_count < 100) {
        // Parse timestamp
        if (line[0] == '[') {
            int minutes, seconds, milliseconds;
            if (sscanf(line, "[%d:%d.%d]", &minutes, &seconds, &milliseconds) == 3) {
                g_lyrics[g_lyric_count].time = minutes * 60000 + seconds * 1000 + milliseconds;
                // Extract lyric content
                char* text_start = strchr(line, ']');
                if (text_start) {
                    strncpy(g_lyrics[g_lyric_count].text, text_start + 1, sizeof(g_lyrics[g_lyric_count].text) - 1);
                    g_lyric_count++;
                }
            }
        }
        line = strtok(NULL, "\n");
    }
}

static int get_music_position(void)
{
    // Get current music playback position (milliseconds)
    char buffer[256];
    
    // Directly get playback position using mciSendStringA
    if (mciSendStringA("status dixon_song position", buffer, sizeof(buffer), NULL) == 0) {
        return atoi(buffer);
    }
    
    return 0;
}

static void update_lyrics(void)
{
    // Update lyrics display based on current playback position
    if (!g_music_playing || g_lyric_count == 0) return;
    
    int position = get_music_position();
    
    // Update progress bar and time label
    if (g_music_progress && g_music_time_label) {
        // Get total music length
        char length_buffer[256];
        mciSendStringA("status dixon_song length", length_buffer, sizeof(length_buffer), NULL);
        int total_length = atoi(length_buffer);
        
        if (total_length > 0) {
            // Calculate progress percentage
            int progress = (position * 100) / total_length;
            // Update progress bar
            SendMessageA(g_music_progress, PBM_SETRANGE, 0, MAKELPARAM(0, 100));
            SendMessageA(g_music_progress, PBM_SETPOS, progress, 0);
            // Force redraw of progress bar
            RedrawWindow(g_music_progress, NULL, NULL, RDW_FRAME | RDW_INVALIDATE | RDW_UPDATENOW);
            
            // Update time label
            int current_minutes = position / 60000;
            int current_seconds = (position % 60000) / 1000;
            int total_minutes = total_length / 60000;
            int total_seconds = (total_length % 60000) / 1000;
            char time_buffer[32];
            sprintf(time_buffer, "%02d:%02d / %02d:%02d", current_minutes, current_seconds, total_minutes, total_seconds);
            SetWindowTextA(g_music_time_label, time_buffer);
        }
    }
    
    // Find current lyric to display
    int new_lyric = -1;
    for (int i = 0; i < g_lyric_count; i++) {
        if (g_lyrics[i].time <= position) {
            new_lyric = i;
        } else {
            break;
        }
    }
    
    // Always update display to ensure sync after seeking
    g_current_lyric = new_lyric;
    if (g_current_lyric >= 0 && g_current_lyric < g_lyric_count) {
        // Display five lines of lyrics, with current lyric in the middle
        char lyric_text[2048] = "";
        
        // Calculate start index
        int start_index = g_current_lyric - 2;
        if (start_index < 0) start_index = 0;
        
        // Calculate end index
        int end_index = g_current_lyric + 2;
        if (end_index >= g_lyric_count) end_index = g_lyric_count - 1;
        
        // Add empty lines if needed at the beginning
        int empty_lines_before = 2 - (g_current_lyric - start_index);
        for (int i = 0; i < empty_lines_before; i++) {
            strcat(lyric_text, "\r\n");
        }
        
        // Display lyrics lines
        for (int i = start_index; i <= end_index; i++) {
            strcat(lyric_text, g_lyrics[i].text);
            if (i < end_index) {
                strcat(lyric_text, "\r\n");
            }
        }
        
        // Add empty lines if needed at the end
        int empty_lines_after = 2 - (end_index - g_current_lyric);
        for (int i = 0; i < empty_lines_after; i++) {
            strcat(lyric_text, "\r\n");
        }
        
        set_edit_text(g_music_lyrics, lyric_text);
    }
}

static LRESULT CALLBACK panel_wnd_proc(HWND hwnd, UINT msg, WPARAM wparam, LPARAM lparam)
{
    switch (msg) {
        case WM_COMMAND:
            // Pass button click events to main window
            if (g_main_window) {
                SendMessageA(g_main_window, WM_COMMAND, wparam, lparam);
                return 0;
            }
            break;
        case WM_LBUTTONDOWN: {
            // Check if click is on progress bar
            POINT pt;
            pt.x = LOWORD(lparam);
            pt.y = HIWORD(lparam);
            
            // Get progress bar window rectangle relative to panel
            RECT progress_rect;
            GetWindowRect(g_music_progress, &progress_rect);
            MapWindowPoints(NULL, hwnd, (LPPOINT)&progress_rect, 2);
            
            // Check if click is within progress bar
            if (PtInRect(&progress_rect, pt)) {
                // Convert point to progress bar client coordinates
                MapWindowPoints(hwnd, g_music_progress, &pt, 1);
                
                // Get progress bar client area
                RECT rect;
                GetClientRect(g_music_progress, &rect);
                
                // Calculate progress percentage at click position
                int progress = (pt.x * 100) / (rect.right - rect.left);
                
                // Get total music length
                char length_buffer[256];
                mciSendStringA("status dixon_song length", length_buffer, sizeof(length_buffer), NULL);
                int total_length = atoi(length_buffer);
                
                if (total_length > 0) {
                    // Set time format to milliseconds before seeking
                    mciSendStringA("set dixon_song time format milliseconds", NULL, 0, NULL);
                    
                    // Calculate target position (in milliseconds)
                    int target_position = (progress * total_length) / 100;
                    
                    // Set music playback position
                    char command[256];
                    sprintf(command, "seek dixon_song to %d", target_position);
                    mciSendStringA(command, NULL, 0, NULL);
                    
                    // Resume playback if already playing, or start if not
                    mciSendStringA("play dixon_song", NULL, 0, NULL);
                    g_music_playing = 1;
                    SetTimer(g_main_window, ID_MUSIC_TIMER, 100, NULL);
                    
                    // Update lyrics immediately after seeking
                    g_current_lyric = -1; // Force update
                    update_lyrics();
                }
                
                return 0;
            }
            break;
        }
        default:
            break;
    }
    
    // Call original window procedure for all other messages
    if (g_original_panel_proc) {
        return CallWindowProcA(g_original_panel_proc, hwnd, msg, wparam, lparam);
    }
    return DefWindowProcA(hwnd, msg, wparam, lparam);
}

static void play_music(void)
{
    // Play music
    char music_path[MAX_PATH];
    char command[512];
    char error_buffer[256];
    
    // Output debug information
    OutputDebugStringA("Play button clicked\n");
    
    // Get current directory and build file path
    GetCurrentDirectoryA(MAX_PATH, music_path);
    OutputDebugStringA("Current directory: ");
    OutputDebugStringA(music_path);
    OutputDebugStringA("\n");
    
    strcat(music_path, "\\DixonSong.mp3");
    OutputDebugStringA("Music path: ");
    OutputDebugStringA(music_path);
    OutputDebugStringA("\n");
    
    // Check if file exists
    if (GetFileAttributesA(music_path) == INVALID_FILE_ATTRIBUTES) {
        OutputDebugStringA("File does not exist\n");
        MessageBoxA(g_main_window, "DixonSong.mp3 not found", "Error", MB_OK | MB_ICONERROR);
        return;
    } else {
        OutputDebugStringA("File exists\n");
    }
    
    // Check if music is already open but paused
    char status_buffer[256];
    mciSendStringA("status dixon_song mode", status_buffer, sizeof(status_buffer), NULL);
    
    if (strcmp(status_buffer, "paused") == 0) {
        // If music is already open and paused, resume playback
        OutputDebugStringA("Continuing paused music\n");
        MCIERROR err = mciSendStringA("resume dixon_song", NULL, 0, NULL);
        if (err != 0) {
            mciGetErrorStringA(err, error_buffer, sizeof(error_buffer));
            OutputDebugStringA("Resume error: ");
            OutputDebugStringA(error_buffer);
            OutputDebugStringA("\n");
            MessageBoxA(g_main_window, error_buffer, "Error resuming music", MB_OK | MB_ICONERROR);
            return;
        } else {
            OutputDebugStringA("Resume successful\n");
        }
    } else {
        // Stop any previously playing music
        mciSendStringA("stop dixon_song", NULL, 0, NULL);
        mciSendStringA("close dixon_song", NULL, 0, NULL);
        
        // Play MP3 file using mciSendString
        snprintf(command, sizeof(command), "open \"%s\" type mpegvideo alias dixon_song", music_path);
        OutputDebugStringA("Open command: ");
        OutputDebugStringA(command);
        OutputDebugStringA("\n");
        
        MCIERROR err = mciSendStringA(command, NULL, 0, NULL);
        if (err != 0) {
            mciGetErrorStringA(err, error_buffer, sizeof(error_buffer));
            OutputDebugStringA("Open error: ");
            OutputDebugStringA(error_buffer);
            OutputDebugStringA("\n");
            MessageBoxA(g_main_window, error_buffer, "Error opening music", MB_OK | MB_ICONERROR);
            return;
        } else {
            OutputDebugStringA("Open successful\n");
        }
        
        // Set time format to milliseconds
        mciSendStringA("set dixon_song time format milliseconds", NULL, 0, NULL);
        
        snprintf(command, sizeof(command), "play dixon_song");
        OutputDebugStringA("Play command: ");
        OutputDebugStringA(command);
        OutputDebugStringA("\n");
        
        err = mciSendStringA(command, NULL, 0, NULL);
        if (err != 0) {
            mciGetErrorStringA(err, error_buffer, sizeof(error_buffer));
            OutputDebugStringA("Play error: ");
            OutputDebugStringA(error_buffer);
            OutputDebugStringA("\n");
            MessageBoxA(g_main_window, error_buffer, "Error playing music", MB_OK | MB_ICONERROR);
            mciSendStringA("close dixon_song", NULL, 0, NULL);
            return;
        } else {
            OutputDebugStringA("Play successful\n");
        }
    }
    
    g_music_playing = 1; // Mark music as playing
    g_current_lyric = 0; // Start from first lyric
    
    // Show initial lyrics (first five lines) immediately when Play is clicked
    if (g_lyric_count > 0 && g_music_lyrics) {
        char lyric_text[2048] = "";
        
        // Display five lines of lyrics, with current lyric in the middle
        // Calculate start index
        int start_index = 0;
        // Calculate end index
        int end_index = g_current_lyric + 2;
        if (end_index >= g_lyric_count) end_index = g_lyric_count - 1;
        
        // Add empty lines at the beginning (2 lines)
        strcat(lyric_text, "\r\n\r\n");
        
        // Display lyrics lines
        for (int i = start_index; i <= end_index; i++) {
            strcat(lyric_text, g_lyrics[i].text);
            if (i < end_index) {
                strcat(lyric_text, "\r\n");
            }
        }
        
        // Add empty lines if needed at the end
        int empty_lines_after = 2 - (end_index - g_current_lyric);
        for (int i = 0; i < empty_lines_after; i++) {
            strcat(lyric_text, "\r\n");
        }
        
        set_edit_text(g_music_lyrics, lyric_text);
        // Force redraw of lyrics display
        RedrawWindow(g_music_lyrics, NULL, NULL, RDW_FRAME | RDW_INVALIDATE | RDW_UPDATENOW);
    }
    
    // Start lyrics update timer
    if (g_main_window) {
        SetTimer(g_main_window, ID_MUSIC_TIMER, 100, NULL);
        OutputDebugStringA("Timer started\n");
    } else {
        OutputDebugStringA("Main window is NULL\n");
    }
}

static void pause_music(void)
{
    // Pause music
    if (g_music_playing) {
        mciSendStringA("pause dixon_song", NULL, 0, NULL);
        g_music_playing = 0;
        
        // Stop lyrics update timer
        KillTimer(g_main_window, ID_MUSIC_TIMER);
    }
}

static void stop_music(void)
{
    // Stop music
    if (g_music_playing) {
        mciSendStringA("stop dixon_song", NULL, 0, NULL);
        mciSendStringA("close dixon_song", NULL, 0, NULL);
        g_music_playing = 0;
        
        // Stop lyrics update timer
        KillTimer(g_main_window, ID_MUSIC_TIMER);
    }
}

static void clear_all(void)
{
    set_edit_text(g_output, "");
    set_edit_text(g_output_result, "");
    set_status("Ready");
}

static DWORD WINAPI worker_thread_proc(LPVOID param)
{
    worker_result_t *result = (worker_result_t *) param;
    dixon_gui_response_init(&result->response);
    result->ok = dixon_gui_run_request(&result->request, &result->response,
                                       result->error, sizeof(result->error));
    PostMessageA(g_main_window, WM_DIXON_COMPLETE, 0, (LPARAM) result);
    return 0;
}

static void begin_request(void)
{
    int current;

    show_active_tab();
    current = TabCtrl_GetCurSel(g_tab);
    if (current < 0 || current >= TAB_COUNT) current = 0;
    worker_result_t *result;
    char *path;

    if (g_is_running) return;

    result = (worker_result_t *) calloc(1, sizeof(worker_result_t));
    if (!result) {
        MessageBoxA(g_main_window, "Out of memory while creating the worker request.", APP_TITLE, MB_OK | MB_ICONERROR);
        return;
    }

    path = get_edit_text(g_path_edit);
    result->request.dixon_path = path;
    g_requested_step1_method = get_method_combo_value(g_step1_combo);
    g_requested_step4_method = get_method_combo_value(g_step4_combo);

    switch (current) {
        case 0:
            result->request.mode = DIXON_GUI_MODE_SOLVE;
            result->request.polynomials = get_edit_text(g_tabs[0].edit1);
            result->request.field = get_edit_text(g_tabs[0].edit2);
            result->request.solve_verbose = 0;
            break;
        case 1:
            result->request.mode = DIXON_GUI_MODE_RESULTANT;
            result->request.polynomials = get_edit_text(g_tabs[1].edit1);
            result->request.eliminate_vars = get_edit_text(g_tabs[1].edit2);
            result->request.field = get_edit_text(g_tabs[1].edit3);
            break;
        case 2:
            result->request.mode = DIXON_GUI_MODE_COMPLEXITY;
            result->request.polynomials = get_edit_text(g_tabs[2].edit1);
            result->request.eliminate_vars = get_edit_text(g_tabs[2].edit2);
            result->request.field = get_edit_text(g_tabs[2].edit3);
            result->request.omega = get_edit_text(g_tabs[2].edit4);
            break;
        case 3:
            result->request.mode = DIXON_GUI_MODE_IDEAL;
            result->request.ideal_generators = get_edit_text(g_tabs[3].edit1);
            result->request.polynomials = get_edit_text(g_tabs[3].edit2);
            result->request.eliminate_vars = get_edit_text(g_tabs[3].edit3);
            result->request.field = get_edit_text(g_tabs[3].edit4);
            break;
        case 4:
            result->request.mode = DIXON_GUI_MODE_RANDOM;
            result->request.random_degrees = get_edit_text(g_tabs[4].edit1);
            result->request.field = get_edit_text(g_tabs[4].edit2);
            result->request.random_mode = (int) SendMessageA(g_tabs[4].edit3, CB_GETCURSEL, 0, 0);
            result->request.omega = get_edit_text(g_tabs[4].edit4);
            if (result->request.random_mode == CB_ERR) {
                result->request.random_mode = DIXON_GUI_RANDOM_RESULTANT;
            }
            break;
        default:
            free(path);
            free(result);
            return;
    }

    g_worker_thread = CreateThread(NULL, 0, worker_thread_proc, result, 0, NULL);
    if (!g_worker_thread) {
        MessageBoxA(g_main_window, "Failed to start the background worker.", APP_TITLE, MB_OK | MB_ICONERROR);
        free((void *) result->request.dixon_path);
        free((void *) result->request.polynomials);
        free((void *) result->request.eliminate_vars);
        free((void *) result->request.field);
        free((void *) result->request.ideal_generators);
        free((void *) result->request.omega);
        free((void *) result->request.random_degrees);
        free(result);
        return;
    }

    g_is_running = 1;
    EnableWindow(g_run_button, FALSE);
    EnableWindow(g_save_button, FALSE);
    set_status("Running drsolve.exe...");
}

static void finish_request(worker_result_t *result)
{
    if (result->ok) {
        set_edit_text(g_output, result->response.combined_output);
        update_result_output(result);
        if (result->response.exit_code == 0) {
            set_status("Finished successfully");
        } else {
            set_status("drsolve.exe returned a non-zero exit code");
        }
    } else {
        set_edit_text(g_output, result->error);
        set_edit_text(g_output_result, "");
        set_status("Failed to run drsolve.exe");
    }

    if (g_worker_thread) {
        CloseHandle(g_worker_thread);
        g_worker_thread = NULL;
    }

    g_is_running = 0;
    EnableWindow(g_run_button, TRUE);
    EnableWindow(g_save_button, TRUE);

    free((void *) result->request.dixon_path);
    free((void *) result->request.polynomials);
    free((void *) result->request.eliminate_vars);
    free((void *) result->request.field);
    free((void *) result->request.ideal_generators);
    free((void *) result->request.omega);
    free((void *) result->request.random_degrees);
    dixon_gui_response_clear(&result->response);
    free(result);
}

static LRESULT CALLBACK main_wnd_proc(HWND hwnd, UINT msg, WPARAM wparam, LPARAM lparam)
{
    switch (msg) {
        case WM_CREATE: {
            char default_path[MAX_PATH];
            INITCOMMONCONTROLSEX icc;
            LOGFONTA mono = {0};

            icc.dwSize = sizeof(icc);
            icc.dwICC = ICC_TAB_CLASSES | ICC_BAR_CLASSES;
            InitCommonControlsEx(&icc);

            g_ui_font = (HFONT) GetStockObject(DEFAULT_GUI_FONT);

            mono.lfHeight = -16;
            strcpy(mono.lfFaceName, "Consolas");
            g_mono_font = CreateFontIndirectA(&mono);
            if (!g_mono_font) g_mono_font = g_ui_font;

            g_path_label = create_label(hwnd, "CLI Path", 0, 0, 10, 10);
            g_path_edit = create_edit(hwnd, ID_PATH_EDIT, ES_LEFT | WS_TABSTOP, 0, 0, 10, 10, "");
            create_button(hwnd, ID_PATH_BROWSE, "Browse...", 0, 0, 10, 10);
            g_step1_label = create_label(hwnd, "Step 1", 0, 0, 10, 10);
            g_step1_combo = create_combobox(hwnd, ID_STEP1_METHOD, 0, 0, 10, 120);
            populate_method_combo(g_step1_combo);
            g_step4_label = create_label(hwnd, "Step 4", 0, 0, 10, 10);
            g_step4_combo = create_combobox(hwnd, ID_STEP4_METHOD, 0, 0, 10, 120);
            populate_method_combo(g_step4_combo);

            g_input_group = create_groupbox(hwnd, "Mode && Parameters");
            g_tab = CreateWindowExA(0, WC_TABCONTROLA, "",
                                    WS_CHILD | WS_VISIBLE | WS_CLIPSIBLINGS | WS_TABSTOP,
                                    0, 0, 10, 10,
                                    g_input_group, (HMENU) (INT_PTR) ID_TAB, g_instance, NULL);
            apply_font(g_tab, g_ui_font);
            g_tab_wndproc = (WNDPROC) SetWindowLongPtrA(g_tab, GWLP_WNDPROC, (LONG_PTR) tab_wnd_proc);
            create_tab_pages();

            g_run_button = create_button(hwnd, ID_RUN, "Run", 0, 0, 10, 10);
            create_button(hwnd, ID_CLEAR, "Clear Output", 0, 0, 10, 10);
            g_save_button = create_button(hwnd, ID_SAVE, "Save Output", 0, 0, 10, 10);

            g_output_group = create_groupbox(hwnd, "Output");
            g_output = create_edit(g_output_group, ID_OUTPUT,
                                   ES_LEFT | ES_MULTILINE | ES_AUTOVSCROLL | ES_READONLY | WS_VSCROLL | WS_TABSTOP,
                                   0, 0, 10, 10, "");
            g_output_result = create_edit(g_output_group, ID_OUTPUT_RESULT,
                                          ES_LEFT | ES_MULTILINE | ES_AUTOVSCROLL | ES_READONLY | WS_VSCROLL | WS_TABSTOP,
                                          0, 0, 10, 10, "");

            g_status = CreateWindowExA(0, STATUSCLASSNAMEA, "", WS_CHILD | WS_VISIBLE,
                                       0, 0, 0, 0, hwnd, (HMENU) (INT_PTR) ID_STATUS, g_instance, NULL);
            apply_font(g_status, g_ui_font);

            if (dixon_gui_find_default_cli(default_path, sizeof(default_path))) {
                set_edit_text(g_path_edit, default_path);
            }
            SetTimer(hwnd, ID_TAB_TIMER, 100, NULL);
            set_status("Ready");
            return 0;
        }
        case WM_SIZE:
            layout_main_window(hwnd);
            return 0;
        case WM_TIMER:
            if (wparam == ID_TAB_TIMER) {
                int tab_index = TabCtrl_GetCurSel(g_tab);
                if (tab_index < 0 || tab_index >= TAB_COUNT) tab_index = 0;
                if (tab_index != g_current_tab_index) show_active_tab();
                return 0;
            } else if (wparam == ID_MUSIC_TIMER) {
                update_lyrics();
                return 0;
            }
            break;
        case WM_COMMAND:
            switch (LOWORD(wparam)) {
                case ID_PATH_BROWSE:
                    browse_for_dixon();
                    return 0;
                case ID_RUN:
                    begin_request();
                    return 0;
                case ID_CLEAR:
                    clear_all();
                    return 0;
                case ID_SAVE:
                    save_output();
                    return 0;
                case ID_RAND_MODE:
                    if (HIWORD(wparam) == CBN_SELCHANGE) {
                        update_random_mode_controls();
                        return 0;
                    }
                    break;
                case ID_MUSIC_PLAY:
                    play_music();
                    return 0;
                case ID_MUSIC_PAUSE:
                    pause_music();
                    return 0;
                case ID_MUSIC_STOP:
                    stop_music();
                    return 0;
                default:
                    break;
            }
            break;
        case WM_NOTIFY:
            if (((LPNMHDR) lparam)->idFrom == ID_TAB && ((LPNMHDR) lparam)->code == TCN_SELCHANGE) {
                show_active_tab();
                return 0;
            }
            break;
        case WM_DIXON_COMPLETE:
            finish_request((worker_result_t *) lparam);
            return 0;
        case WM_DESTROY:
            KillTimer(hwnd, ID_TAB_TIMER);
            if (g_worker_thread) {
                WaitForSingleObject(g_worker_thread, INFINITE);
                CloseHandle(g_worker_thread);
                g_worker_thread = NULL;
            }
            if (g_mono_font && g_mono_font != g_ui_font) DeleteObject(g_mono_font);
            PostQuitMessage(0);
            return 0;
        default:
            break;
    }

    return DefWindowProcA(hwnd, msg, wparam, lparam);
}

int WINAPI WinMain(HINSTANCE instance, HINSTANCE prev_instance, LPSTR cmd_line, int show_cmd)
{
    WNDCLASSA wc;
    MSG msg;

    (void) prev_instance;
    (void) cmd_line;

    g_instance = instance;

    ZeroMemory(&wc, sizeof(wc));
    wc.lpfnWndProc = main_wnd_proc;
    wc.hInstance = instance;
    wc.lpszClassName = "DRsolveWinGuiMainWindow";
    wc.hCursor = LoadCursor(NULL, IDC_ARROW);
    wc.hbrBackground = (HBRUSH) (COLOR_WINDOW + 1);
    wc.hIcon = LoadIcon(NULL, IDI_APPLICATION);

    RegisterClassA(&wc);

    g_main_window = CreateWindowExA(
        0,
        wc.lpszClassName,
        APP_TITLE,
        WS_OVERLAPPEDWINDOW | WS_CLIPCHILDREN,
        CW_USEDEFAULT,
        CW_USEDEFAULT,
        1120,
        920,
        NULL,
        NULL,
        instance,
        NULL);

    if (!g_main_window) return 1;

    ShowWindow(g_main_window, show_cmd);
    UpdateWindow(g_main_window);

    while (GetMessageA(&msg, NULL, 0, 0) > 0) {
        TranslateMessage(&msg);
        DispatchMessageA(&msg);
    }

    return (int) msg.wParam;
}