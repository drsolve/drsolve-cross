#include <windows.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define DIXON_PATH_CAP 4096

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

    snprintf(buffer, size, "%s: %s", prefix, system_message);
}

static int get_module_dir(char *buffer, size_t size)
{
    DWORD len;
    char *slash;

    len = GetModuleFileNameA(NULL, buffer, (DWORD) size);
    if (len == 0 || len >= size) return 0;

    slash = strrchr(buffer, '\\');
    if (!slash) slash = strrchr(buffer, '/');
    if (!slash) return 0;

    *slash = '\0';
    return 1;
}

static int build_path(char *buffer, size_t size, const char *dir, const char *tail)
{
    return snprintf(buffer, size, "%s\\%s", dir, tail) < (int) size;
}

static char *dup_environment_value(const char *name)
{
    DWORD needed = GetEnvironmentVariableA(name, NULL, 0);
    char *value;

    if (needed == 0) {
        if (GetLastError() == ERROR_ENVVAR_NOT_FOUND) return NULL;
        return NULL;
    }

    value = (char *) malloc(needed);
    if (!value) return NULL;

    if (GetEnvironmentVariableA(name, value, needed) != needed - 1) {
        free(value);
        return NULL;
    }

    return value;
}

static const char *skip_program_name(const char *command_line)
{
    const char *p = command_line;

    while (*p == ' ' || *p == '\t') p++;

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
        while (*p && *p != ' ' && *p != '\t') p++;
    }

    while (*p == ' ' || *p == '\t') p++;
    return p;
}

static int launcher_owns_console(void)
{
    DWORD processes[8];
    DWORD count = GetConsoleProcessList(processes, 8);
    return count == 1;
}

static int configure_child_path(const char *dll_dir, char *error_message, size_t error_message_size)
{
    char *old_path = dup_environment_value("PATH");
    size_t old_len = old_path ? strlen(old_path) : 0;
    size_t dll_len = strlen(dll_dir);
    char *new_path = (char *) malloc(dll_len + (old_len ? old_len + 2 : 1));
    int ok;

    if (!new_path) {
        snprintf(error_message, error_message_size,
                 "Out of memory while preparing PATH for runtime DLLs.");
        free(old_path);
        return 0;
    }

    if (old_len > 0) {
        snprintf(new_path, dll_len + old_len + 2, "%s;%s", dll_dir, old_path);
    } else {
        snprintf(new_path, dll_len + 1, "%s", dll_dir);
    }

    ok = SetEnvironmentVariableA("PATH", new_path);
    if (!ok) {
        set_win32_error(error_message, error_message_size,
                        "Failed to configure PATH for runtime DLLs", GetLastError());
    }

    free(new_path);
    free(old_path);
    return ok ? 1 : 0;
}

static int launch_interactive_shell(void)
{
    STARTUPINFOA si;
    PROCESS_INFORMATION pi;
    char command_line[64] = "cmd.exe /K";

    ZeroMemory(&si, sizeof(si));
    GetStartupInfoA(&si);
    si.cb = sizeof(si);
    ZeroMemory(&pi, sizeof(pi));

    if (!CreateProcessA(NULL,
                        command_line,
                        NULL,
                        NULL,
                        TRUE,
                        0,
                        NULL,
                        NULL,
                        &si,
                        &pi)) {
        return 0;
    }

    WaitForSingleObject(pi.hProcess, INFINITE);
    CloseHandle(pi.hThread);
    CloseHandle(pi.hProcess);
    return 1;
}

int main(void)
{
    char root_dir[DIXON_PATH_CAP];
    char dll_dir[DIXON_PATH_CAP];
    char real_cli[DIXON_PATH_CAP];
    const char *tail = skip_program_name(GetCommandLineA());
    int launch_shell_after_exit = (*tail == '\0' && launcher_owns_console());
    char *command_line = NULL;
    size_t needed;
    STARTUPINFOA si;
    PROCESS_INFORMATION pi;
    DWORD exit_code = 1;
    DWORD creation_flags = 0;
    char error_message[512];

    if (!get_module_dir(root_dir, sizeof(root_dir))) {
        fprintf(stderr, "Failed to resolve the dixon.exe directory.\n");
        return 1;
    }

    if (!build_path(dll_dir, sizeof(dll_dir), root_dir, "dll")) {
        fprintf(stderr, "The DLL directory path is too long.\n");
        return 1;
    }

    if (!build_path(real_cli, sizeof(real_cli), root_dir, "bin\\dixon_cli_real.exe")) {
        fprintf(stderr, "The internal CLI path is too long.\n");
        return 1;
    }

    if (!configure_child_path(dll_dir, error_message, sizeof(error_message))) {
        fprintf(stderr, "%s\n", error_message);
        return 1;
    }

    SetEnvironmentVariableA("DIXON_DISPLAY_NAME", "dixon.exe");

    needed = strlen(real_cli) + strlen(tail) + 4;
    command_line = (char *) malloc(needed);
    if (!command_line) {
        fprintf(stderr, "Out of memory while building the launcher command line.\n");
        return 1;
    }

    if (*tail) {
        snprintf(command_line, needed, "\"%s\" %s", real_cli, tail);
    } else {
        snprintf(command_line, needed, "\"%s\"", real_cli);
    }

    ZeroMemory(&si, sizeof(si));
    GetStartupInfoA(&si);
    si.cb = sizeof(si);
    if (!(si.dwFlags & STARTF_USESTDHANDLES)) {
        si.dwFlags |= STARTF_USESTDHANDLES;
        si.hStdInput = GetStdHandle(STD_INPUT_HANDLE);
        si.hStdOutput = GetStdHandle(STD_OUTPUT_HANDLE);
        si.hStdError = GetStdHandle(STD_ERROR_HANDLE);
    }
    ZeroMemory(&pi, sizeof(pi));

    if (GetConsoleWindow() == NULL) {
        creation_flags |= CREATE_NO_WINDOW;
    }

    if (!CreateProcessA(real_cli,
                        command_line,
                        NULL,
                        NULL,
                        TRUE,
                        creation_flags,
                        NULL,
                        NULL,
                        &si,
                        &pi)) {
        set_win32_error(error_message, sizeof(error_message),
                        "Failed to start the internal CLI", GetLastError());
        fprintf(stderr, "%s\n", error_message);
        free(command_line);
        return 1;
    }

    free(command_line);

    WaitForSingleObject(pi.hProcess, INFINITE);
    if (!GetExitCodeProcess(pi.hProcess, &exit_code)) {
        set_win32_error(error_message, sizeof(error_message),
                        "Failed to read the internal CLI exit code", GetLastError());
        fprintf(stderr, "%s\n", error_message);
        exit_code = 1;
    }

    CloseHandle(pi.hThread);
    CloseHandle(pi.hProcess);

    if (launch_shell_after_exit) {
        launch_interactive_shell();
    }

    return (int) exit_code;
}
