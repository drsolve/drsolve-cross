#include <windows.h>
#include <commctrl.h>
#include <commdlg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "dixon_wrapper.h"

#define APP_TITLE "Dixon Windows GUI"
#define TAB_COUNT 5

#define ID_TAB 100
#define ID_PATH_EDIT 102
#define ID_PATH_BROWSE 103
#define ID_RUN 104
#define ID_CLEAR 105
#define ID_SAVE 106
#define ID_STATUS 107
#define ID_OUTPUT 108
#define ID_OUTPUT_RESULT 109

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

#define WM_DIXON_COMPLETE (WM_APP + 1)
#define ID_TAB_TIMER 1

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

static void show_active_tab(void);

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

    if (!join_result_path(result->request.dixon_path, saved_path, full_path, sizeof(full_path))) {
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
    int row_h = 28;
    int browse_w = 112;
    int label_w = 62;
    int input_h = 420;
    int button_y;
    int output_y;
    int output_h;
    int status_h = 24;
    int content_width;

    GetClientRect(hwnd, &rc);
    width = rc.right - rc.left;
    height = rc.bottom - rc.top;
    content_width = width - margin * 2;

    if (height < 760) input_h = 376;
    if (height < 680) input_h = 336;
    if (height > 920) input_h = 456;

    if (input_h > height - 250) input_h = height - 250;
    if (input_h < 260) input_h = 260;

    MoveWindow(g_path_label, margin, path_y + 4, label_w, 20, TRUE);
    MoveWindow(g_path_edit, margin + label_w + 8, path_y,
               content_width - label_w - browse_w - 24,
               row_h, TRUE);
    MoveWindow(GetDlgItem(hwnd, ID_PATH_BROWSE), width - margin - browse_w, path_y - 1,
               browse_w, row_h + 2, TRUE);

    MoveWindow(g_input_group, margin, 58, content_width, input_h, TRUE);
    MoveWindow(g_tab, 12, 30, content_width - 24, input_h - 42, TRUE);
    layout_tab_pages();

    button_y = 58 + input_h + 14;
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
        set_status("Using selected dixon.exe");
    }
}

static void save_output(void)
{
    OPENFILENAMEA ofn;
    char path[MAX_PATH] = "dixon_win_output.txt";
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
    set_status("Running dixon.exe...");
}

static void finish_request(worker_result_t *result)
{
    if (result->ok) {
        set_edit_text(g_output, result->response.combined_output);
        update_result_output(result);
        if (result->response.exit_code == 0) {
            set_status("Finished successfully");
        } else {
            set_status("dixon.exe returned a non-zero exit code");
        }
    } else {
        set_edit_text(g_output, result->error);
        set_edit_text(g_output_result, "");
        set_status("Failed to run dixon.exe");
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
    wc.lpszClassName = "DixonWinGuiMainWindow";
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
