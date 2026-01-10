/**
 * \file traceback.c
 * \brief Function definitions for traceback logging.
 *
 * \author Ching-Yin Ng
 */

#include <stdbool.h>
#include <stdio.h>
#include <string.h>

#include "c_traceback.h"
#include "internal/trace.h"
#include "internal/utils.h"

#if defined(_WIN32)
#include <io.h>
#define SAFE_WRITE(fd, buf, len) _write(fd, buf, len)
#define STDERR_FD 2
#else
#include <unistd.h>
#define SAFE_WRITE(fd, buf, len) write(fd, buf, len)
#define STDERR_FD STDERR_FILENO
#endif

typedef struct
{
    const char *reset;
    const char *bold;
    const char *error;
    const char *error_bold;
    const char *tb_text;
    const char *tb_counter;
    const char *tb_file;
    const char *tb_line;
    const char *tb_func;
    const char *tb_another_exception;
    const char *theme_bold;
    const char *theme;
} Theme;

static const char *LOGO_LINES[] = {
    "    %%%%%%%%%%%%    ",
    "  %%%%%%%%%%%%%%%%  ",
    " %%%%%%%%%%%%%%%%%% ",
    "%%%%%%%%%%*  %%%%%%%",
    "%%%%%%%*     %%%%%%%",
    "%%%%%*       %%%%%%%",
    "%%%%%%%*     %%%%%%%",
    "%%%%%%%%%%*  %%%%%%%",
    " %%%%%%%%%%%%%%%%%% ",
    "  %%%%%%%%%%%%%%%%  ",
    "    %%%%%%%%%%%%    "
};

/**
 * \brief Populates the Theme struct based on color availability.
 *
 * \param[in] use_color Whether to use color in the output.
 * \return The populated Theme struct.
 */
static Theme get_theme(bool use_color)
{
    if (!use_color)
    {
        return (Theme){"", "", "", "", "", "", "", "", "", "", "", ""};
    }

    return (Theme
    ){.reset = CTB_RESET_COLOR,
      .bold = CTB_THEME_BOLD_COLOR,
      .error = CTB_ERROR_COLOR,
      .error_bold = CTB_ERROR_BOLD_COLOR,
      .tb_text = CTB_TRACEBACK_TEXT_COLOR,
      .tb_counter = CTB_TRACEBACK_COUNTER_COLOR,
      .tb_file = CTB_TRACEBACK_FILE_COLOR,
      .tb_line = CTB_TRACEBACK_LINE_COLOR,
      .tb_func = CTB_TRACEBACK_FUNC_COLOR,
      .tb_another_exception = CTB_TRACEBACK_ANOTHER_EXCEPTION_TEXT_COLOR,
      .theme_bold = CTB_THEME_BOLD_COLOR,
      .theme = CTB_THEME_COLOR};
}

/**
 * \brief Returns the appropriate dash character (UTF-8 or ASCII).
 *
 * \param[in] stream The output stream.
 * \return The dash character as a string.
 */
static const char *get_dash(FILE *stream)
{
    return should_use_utf8(stream) ? "\u2500" : "-";
}

/**
 * \brief Helper function to print a single frame.
 *
 * \param[in] stream The output stream.
 * \param[in] index The index of the frame in the call stack.
 * \param[in] frame The frame to print.
 * \param[in] use_color Whether to use color in the output.
 * \param[in] theme The theme to use for coloring.
 * \param[in] source_code_color The color code for the source code.
 */
static void print_frame(
    FILE *stream,
    int index,
    const CTB_Frame_ *frame,
    const bool use_color,
    const Theme *theme,
    const char *source_code_color
)
{
    const int dir_len = get_parent_path_length(frame->filename);

    const char *source_code_color_code = use_color ? source_code_color : "";

    // clang-format off
    fprintf(
        stream,
        "  %s(#%02d)%s %sFile \"%s",
        theme->tb_counter, index, theme->reset,
        theme->tb_text, theme->reset
    );
    // clang-format on

    fprintf(stream, "%s%.*s%s", theme->tb_text, dir_len, frame->filename, theme->reset);
    fprintf(stream, "%s%s%s", theme->tb_file, frame->filename + dir_len, theme->reset);

    // clang-format off
    fprintf(
        stream,
        "%s\", line%s %s%d%s %sin%s %s%s%s:\n    %s%s%s\n",
        theme->tb_text, theme->reset,
        theme->tb_line, frame->line_number, theme->reset,
        theme->tb_text, theme->reset,
        theme->tb_func, frame->function_name, theme->reset,
        source_code_color_code, frame->source_code, theme->reset
    );
    // clang-format on
}

/**
 * \brief Helper function to print a horizontal rule with optional headers.
 *
 * \param[in] stream The output stream.
 * \param[in] use_color Whether to use color in the output.
 * \param[in] color_code The color code to use for the horizontal rule.
 * \param[in] header The header text to display in the middle of the rule.
 */
static void print_hrule_internal(
    FILE *stream,
    const bool use_color,
    const char *restrict color_code,
    const char *restrict header
)
{
    const int terminal_width = get_terminal_width(stream);
    const int max = CTB_HRULE_MAX_WIDTH;
    const int min = CTB_HRULE_MIN_WIDTH;

    int hrule_width = terminal_width;
    if (hrule_width < min)
    {
        hrule_width = min;
    }
    else if (hrule_width > max)
    {
        hrule_width = max;
    }

    const char *dash = get_dash(stream);
    const char *reset = use_color ? CTB_RESET_COLOR : "";
    const char *color = (use_color && color_code) ? color_code : "";

    int header_len = (header) ? (int)strlen(header) : 0;
    int left_width = hrule_width;
    int right_width = 0;

    if (header_len > 0)
    {
        const int padding = 2;
        const int available_space = hrule_width - header_len - padding;

        if (available_space > 0)
        {
            left_width = available_space / 2;
            right_width = available_space - left_width;
        }
        else
        {
            left_width = 2;
            right_width = 2;
        }
    }

    fprintf(stream, "%s", color);

    for (int i = 0; i < left_width; i++)
    {
        fputs(dash, stream);
    }

    if (header_len > 0)
    {
        fprintf(stream, " %s ", header);
    }

    for (int i = 0; i < right_width; i++)
    {
        fputs(dash, stream);
    }
    fprintf(stream, "%s\n", reset);
}

/**
 * \brief Print a horizontal rule without a header.
 *
 * \param[in] stream The output stream.
 * \param[in] use_color Whether to use color in the output.
 * \param[in] color_code The color code to use for the horizontal rule.
 */
static void
print_hrule(FILE *stream, const bool use_color, const char *restrict color_code)
{
    print_hrule_internal(stream, use_color, color_code, NULL);
}

/**
 * \brief Print a horizontal rule with a header.
 *
 * \param[in] stream The output stream.
 * \param[in] use_color Whether to use color in the output.
 * \param[in] color_code The color code to use for the horizontal rule.
 * \param[in] header The header text to display in the middle of the rule.
 */
static void print_hrule_with_header(
    FILE *stream,
    const bool use_color,
    const char *restrict color_code,
    const char *restrict header
)
{
    print_hrule_internal(stream, use_color, color_code, header);
}

void ctb_log_call_stack(void)
{
    const CTB_Context *context = get_context();
    FILE *const stream = stderr;
    const bool use_color = should_use_color(stream);
    const Theme theme = get_theme(use_color);

    const int num_frames = context->call_depth;
    const int num_frames_to_print =
        (num_frames > CTB_MAX_CALL_STACK_DEPTH) ? CTB_MAX_CALL_STACK_DEPTH : num_frames;

    print_hrule(stream, use_color, theme.theme);

    /* Print Header */
    fprintf(
        stream,
        "%sCall Stack (most recent call last):%s\n",
        theme.theme_bold,
        theme.reset
    );

    if (num_frames_to_print <= 0)
    {
        fputs("There is no recorded stack frame!\n", stream);
        print_hrule(stream, use_color, theme.theme);
        fflush(stream);
        return;
    }

    for (int i = 0; i < num_frames_to_print; i++)
    {
        print_frame(
            stream,
            i,
            &context->call_stack_frames[i],
            use_color,
            &theme,
            CTB_THEME_COLOR
        );
    }

    if (num_frames > CTB_MAX_CALL_STACK_DEPTH)
    {
        fprintf(
            stream,
            "\n%s[... Skipped %d frames ...]%s\n",
            theme.tb_text,
            num_frames - CTB_MAX_CALL_STACK_DEPTH,
            theme.reset
        );
    }

    print_hrule(stream, use_color, theme.theme);
    fflush(stream);
}

void ctb_log_traceback(void)
{
    const CTB_Context *context = get_context();
    FILE *const stream = stderr;
    const bool use_color = should_use_color(stream);
    const Theme theme = get_theme(use_color);

    const char *header_text = (CTB_TRACEBACK_HEADER && CTB_TRACEBACK_HEADER[0])
                                  ? CTB_TRACEBACK_HEADER
                                  : "Traceback";

    const int num_errors = context->num_errors;
    const int num_errors_to_print =
        (num_errors > CTB_MAX_NUM_ERROR) ? CTB_MAX_NUM_ERROR : num_errors;

    print_hrule(stream, use_color, CTB_ERROR_COLOR);

    if (num_errors_to_print <= 0)
    {
        fputs("There is no recorded error!\n", stream);
        print_hrule(stream, use_color, CTB_ERROR_COLOR);
        fflush(stream);
        return;
    }

    for (int e = 0; e < num_errors_to_print; e++)
    {
        const CTB_Error_Snapshot_ *snapshot = &context->error_snapshots[e];
        const int num_frames = snapshot->call_depth;
        const bool stack_frames_exceed_max = (num_frames > CTB_MAX_CALL_STACK_DEPTH);
        const int num_frames_to_print =
            stack_frames_exceed_max ? CTB_MAX_CALL_STACK_DEPTH : num_frames;

        /* Print Header */
        if (num_errors > 1)
        {
            fprintf(stream, "%s(#%02d)%s ", theme.error, e, theme.reset);
        }

        fprintf(
            stream,
            "%s%s%s %s(most recent call last):%s\n",
            theme.error_bold,
            header_text,
            theme.reset,
            theme.error,
            theme.reset
        );

        /* Print Stack Frames */
        for (int i = 0; i < num_frames_to_print; i++)
        {
            print_frame(
                stream,
                i,
                &snapshot->call_stack_frames[i],
                use_color,
                &theme,
                CTB_ERROR_COLOR
            );
        }

        if (stack_frames_exceed_max)
        {
            fprintf(
                stream,
                "\n      %s[... Skipped %d frames ...]%s\n\n",
                theme.tb_text,
                num_frames - CTB_MAX_CALL_STACK_DEPTH,
                theme.reset
            );
        }

        print_frame(
            stream,
            num_frames,
            &snapshot->error_frame,
            use_color,
            &theme,
            CTB_ERROR_COLOR
        );

        fprintf(stream, "%s%s", theme.error_bold, error_to_string(snapshot->error));
        if (snapshot->error_message[0])
        {
            fprintf(
                stream,
                ":%s %s%s%s",
                theme.reset,
                theme.error,
                snapshot->error_message,
                theme.reset
            );
        }
        else
        {
            fputs(theme.reset, stream);
        }
        fputs("\n", stream);

        if (e < (num_errors_to_print - 1))
        {
            fprintf(
                stream,
                "\n%sDuring handling of the above exception, another exception "
                "occurred:%s\n\n",
                theme.tb_another_exception,
                theme.reset
            );
        }
    }

    if (num_errors > CTB_MAX_NUM_ERROR)
    {
        fprintf(
            stream,
            "\n%s[... Truncated %d errors ...]%s\n",
            theme.error_bold,
            num_errors - CTB_MAX_NUM_ERROR,
            theme.reset
        );
    }

    print_hrule(stream, use_color, CTB_ERROR_COLOR);
    fflush(stream);
}

void ctb_dump_traceback(void)
{
    ctb_log_traceback();
    ctb_clear_error();
}

/**
 * \brief Helper function to print the left column of the compilation info.
 *
 * \param[in] stream The output stream.
 * \param[in] theme The theme to use for coloring.
 * \param[in] row_idx The current row index.
 * \param[in] label The label text.
 * \param[in] value The value text.
 */
static void print_compilation_info_row(
    FILE *stream,
    const Theme *theme,
    const int row_idx,
    const char *restrict label,
    const char *restrict value
)
{
    const int logo_height = sizeof(LOGO_LINES) / sizeof(LOGO_LINES[0]);
    const int left_padding = 2;
    const int gutter = 4;
    const int logo_width = strlen(LOGO_LINES[0]);

    if (row_idx < logo_height)
    {
        fprintf(
            stream,
            "%*s%s%s%s%*s",
            left_padding,
            "",
            theme->theme_bold,
            LOGO_LINES[row_idx],
            theme->reset,
            gutter,
            ""
        );
    }
    else
    {
        fprintf(stream, "%*s", left_padding + logo_width + gutter, "");
    }

    if (label)
    {
        if (value)
        {
            fprintf(stream, "%s%s%s%s", theme->theme_bold, label, theme->reset, value);
        }
        else
        {
            fprintf(stream, "%s%s%s", theme->theme_bold, label, theme->reset);
        }
    }

    fputs("\n", stream);
}

void ctb_print_compilation_info(void)
{
    FILE *const stream = stdout;
    const bool use_color = should_use_color(stream);
    const Theme theme = get_theme(use_color);
    const char *dash = get_dash(stream);
    const int logo_height = sizeof(LOGO_LINES) / sizeof(LOGO_LINES[0]);

    /* Determine OS String */
    const char *os_str = "Unknown";
#ifdef _WIN32
    os_str = "Windows";
#elif __APPLE__
    os_str = "MacOS";
#elif __linux__
    os_str = "Linux";
#endif

    /* Determine Compiler String */
    char compiler_str[64];
#ifdef _MSC_VER
    snprintf(compiler_str, sizeof(compiler_str), "MSVC (version: %d)", _MSC_VER);
#elif defined(__clang__)
    snprintf(
        compiler_str, sizeof(compiler_str), "Clang (version: %d)", __clang_major__
    );
#elif defined(__GNUC__)
    snprintf(compiler_str, sizeof(compiler_str), "GCC (version: %d)", __GNUC__);
#else
    snprintf(compiler_str, sizeof(compiler_str), "Unknown");
#endif

    /* Print Header */
    print_hrule_with_header(
        stream, use_color, CTB_THEME_COLOR, "C Traceback Compilation Info"
    );

    /* Prepare Config Values for printing */
    char buf_ver[128];
    char buf_date[128];
    char buf_stack[128];
    char buf_msg[128];
    char buf_err[128];
    char buf_term[128];
    char buf_file[128];
    char buf_hmax[128];
    char buf_hmin[128];

    snprintf(buf_ver, sizeof(buf_ver), "%s", CTB_VERSION);
    snprintf(buf_date, sizeof(buf_date), "%s %s", __DATE__, __TIME__);
    snprintf(buf_stack, sizeof(buf_stack), "%d", CTB_MAX_CALL_STACK_DEPTH);
    snprintf(buf_msg, sizeof(buf_msg), "%d", CTB_MAX_ERROR_MESSAGE_LENGTH);
    snprintf(buf_err, sizeof(buf_err), "%d", CTB_MAX_NUM_ERROR);
    snprintf(buf_term, sizeof(buf_term), "%d", CTB_DEFAULT_TERMINAL_WIDTH);
    snprintf(buf_file, sizeof(buf_file), "%d", CTB_DEFAULT_FILE_WIDTH);
    snprintf(buf_hmax, sizeof(buf_hmax), "%d", CTB_HRULE_MAX_WIDTH);
    snprintf(buf_hmin, sizeof(buf_hmin), "%d", CTB_HRULE_MIN_WIDTH);

    /* Construct Horizontal Line */
    char separator_line[64] = {0};
    for (int i = 0; i < 6; i++)
    {
        strcat(separator_line, dash);
    }

    /* Print Info Rows Linearly */
    int row = 0;
    print_compilation_info_row(stream, &theme, row++, "C Traceback Version: ", buf_ver);
    print_compilation_info_row(stream, &theme, row++, "Operating System: ", os_str);
    print_compilation_info_row(stream, &theme, row++, "Build Date: ", buf_date);
    print_compilation_info_row(stream, &theme, row++, "Compiler: ", compiler_str);
    print_compilation_info_row(stream, &theme, row++, "", NULL); /* Empty Row */
    print_compilation_info_row(stream, &theme, row++, "Config", NULL);
    print_compilation_info_row(stream, &theme, row++, separator_line, NULL);
    print_compilation_info_row(
        stream, &theme, row++, "Max Call Stack Depth: ", buf_stack
    );
    print_compilation_info_row(
        stream, &theme, row++, "Max Error Message Length: ", buf_msg
    );
    print_compilation_info_row(
        stream, &theme, row++, "Max Number of Errors: ", buf_err
    );
    print_compilation_info_row(
        stream, &theme, row++, "Default Terminal Width: ", buf_term
    );
    print_compilation_info_row(stream, &theme, row++, "Default File Width: ", buf_file);
    print_compilation_info_row(
        stream, &theme, row++, "Horizontal Rule Max Width: ", buf_hmax
    );
    print_compilation_info_row(
        stream, &theme, row++, "Horizontal Rule Min Width: ", buf_hmin
    );

    /* Fill remaining logo space */
    while (row < logo_height)
    {
        print_compilation_info_row(stream, &theme, row, NULL, NULL);
        row++;
    }

    /* Sample logging */
    fputs("\n", stream);
    fprintf(stream, "%sLogging (example)%s\n", theme.theme_bold, theme.reset);
    for (int d = 0; d < 17; d++)
    {
        fputs(dash, stream);
    }
    fputs("\n", stream);
    fflush(stream);
    LOG_ERROR(CTB_ERROR, "Sample error for compilation info");
    LOG_WARNING(CTB_USER_WARNING, "Sample warning for compilation info");
    LOG_MESSAGE("Sample info for compilation info");

    /* Sample Traceback */
    const int num_examples = 3;
    const CTB_Frame_ example_frames[3] = {
        {10, "example/example.c", "main", "hello_world();"},
        {25, "example/hello_world.c", "check_terminal", "data = compute(data)"},
        {50, "example/libs/utils.c", "compute", "recursion()"}
    };

    const CTB_Frame_ error_frame = {
        75, "example/libs/utils.c", "recursion", "<error thrown here>"
    };

    fputs("\n", stream);
    fprintf(stream, "%sTraceback (example)%s\n", theme.theme_bold, theme.reset);
    for (int d = 0; d < 19; d++)
    {
        fputs(dash, stream);
    }
    fputs("\n", stream);
    fflush(stream);

    const char *header_text = (CTB_TRACEBACK_HEADER && CTB_TRACEBACK_HEADER[0])
                                  ? CTB_TRACEBACK_HEADER
                                  : "Traceback";
    fprintf(
        stream,
        "%s%s%s %s(most recent call last):%s\n",
        theme.error_bold,
        header_text,
        theme.reset,
        theme.error,
        theme.reset
    );

    for (int i = 0; i < num_examples; i++)
    {
        print_frame(stream, i, &example_frames[i], use_color, &theme, CTB_ERROR_COLOR);
    }

    fprintf(
        stream,
        "\n      %s[... Skipped %d frames ...]%s\n\n",
        theme.tb_text,
        123,
        theme.reset
    );

    print_frame(stream, 127, &error_frame, use_color, &theme, CTB_ERROR_COLOR);
    fprintf(
        stream,
        "%s%s:%s %s%s%s\n",
        theme.error_bold,
        error_to_string(CTB_ERROR),
        theme.reset,
        theme.error,
        "Something went wrong!",
        theme.reset
    );

    print_hrule_with_header(stream, use_color, CTB_THEME_COLOR, "END");
    fflush(stream);
}

/**
 * \brief Async-signal-safe string writer.
 */
static void safe_print_str(const char *string)
{
    if (!string)
    {
        return;
    }
    size_t len = 0;
    while (string[len])
    {
        len++;
    }
    SAFE_WRITE(STDERR_FD, string, (unsigned int)len);
}

/**
 * \brief Async-signal-safe integer writer (converts int to string).
 */
static void safe_print_int(const int n)
{
    char buffer[64];
    int i = 0;
    int is_neg = 0;
    long long x = n;

    if (x == 0)
    {
        safe_print_str("0");
        return;
    }

    if (x < 0)
    {
        is_neg = 1;
        x = -x;
    }

    while (x > 0)
    {
        buffer[i++] = (x % 10) + '0';
        x /= 10;
    }

    if (is_neg)
    {
        buffer[i++] = '-';
    }

    // Reverse buffer to print
    for (int j = 0; j < i / 2; j++)
    {
        char temp = buffer[j];
        buffer[j] = buffer[i - j - 1];
        buffer[i - j - 1] = temp;
    }
    buffer[i] = '\0';

    SAFE_WRITE(STDERR_FD, buffer, (unsigned int)i);
}

/**
 * \brief Async-signal-safe helper function to print a single frame.
 *
 * \param[in] index The index of the frame in the call stack.
 * \param[in] frame The frame to print.
 */
static void safe_print_frame(int index, const CTB_Frame_ *frame)
{
    safe_print_str("  (#");
    if (index < 10)
    {
        safe_print_str("0");
    }
    safe_print_int(index);
    safe_print_str(") File \"");
    safe_print_str(frame->filename);
    safe_print_str("\", line ");
    safe_print_int(frame->line_number);
    safe_print_str(" in ");
    safe_print_str(frame->function_name);
    safe_print_str(":\n    ");
    safe_print_str(frame->source_code);
    safe_print_str("\n");
}

void ctb_dump_traceback_signal(const CTB_Error ctb_error)
{
    CTB_Context *context = get_context();
    if (!context)
    {
        safe_print_str("Critical Error: Could not access thread context.\n");
        return;
    }

    const char *header_text = (CTB_TRACEBACK_HEADER && CTB_TRACEBACK_HEADER[0])
                                  ? CTB_TRACEBACK_HEADER
                                  : "Traceback";

    const int num_errors = context->num_errors;
    const int num_errors_to_print =
        (num_errors > CTB_MAX_NUM_ERROR) ? CTB_MAX_NUM_ERROR : num_errors;

    safe_print_str("\n");
    // Red Bold
    for (int i = 0; i < CTB_DEFAULT_TERMINAL_WIDTH; i++)
    {
        safe_print_str("-");
    }
    safe_print_str("\n");

    for (int e = 0; e < num_errors_to_print; e++)
    {
        CTB_Error_Snapshot_ *snapshot = &context->error_snapshots[e];
        const int num_frames = snapshot->call_depth;
        const bool stack_frames_exceed_max = (num_frames > CTB_MAX_CALL_STACK_DEPTH);
        const int num_frames_to_print =
            stack_frames_exceed_max ? CTB_MAX_CALL_STACK_DEPTH : num_frames;

        /* Print Header */
        if (num_errors > 1)
        {
            safe_print_str("(#");
            if (e < 10)
            {
                safe_print_str("0");
            }
            safe_print_int(e);
            safe_print_str(") ");
        }

        safe_print_str(header_text);
        safe_print_str(" (most recent call last):\n");

        /* Print Stack Frames */
        for (int i = 0; i < num_frames_to_print; i++)
        {
            safe_print_frame(i, &snapshot->call_stack_frames[i]);
        }

        if (stack_frames_exceed_max)
        {
            safe_print_str("\n      [... Skipped ");
            safe_print_int(num_frames - CTB_MAX_CALL_STACK_DEPTH);
            safe_print_str(" frames ...]\n\n");
        }

        safe_print_frame(num_frames, &snapshot->error_frame);

        /* Print Error Message */
        safe_print_str(error_to_string(snapshot->error));
        if (snapshot->error_message[0])
        {
            safe_print_str(": ");
            safe_print_str(snapshot->error_message);
        }

        safe_print_str("\n\nDuring handling of the above exception, another exception "
                       "occurred:\n\n");
    }

    if (num_errors > CTB_MAX_NUM_ERROR)
    {
        safe_print_str("\n[... Truncated ");
        safe_print_int(num_errors - CTB_MAX_NUM_ERROR);
        safe_print_str(" errors ...]\n");
    }

    /* Print signal error*/
    if ((num_errors + 1) > 1)
    {
        safe_print_str("(#");
        if ((num_errors + 1) < 10)
        {
            safe_print_str("0");
        }
        safe_print_int(num_errors);
        safe_print_str(") ");
    }

    safe_print_str(header_text);
    safe_print_str(" (most recent call last):\n");

    /* Print Stack Frames */
    const int num_frames = context->call_depth;

    if (num_frames <= 0)
    {
        safe_print_str("  [No recorded stack frames]\n");
    }
    else
    {
        const bool stack_frames_exceed_max = (num_frames > CTB_MAX_CALL_STACK_DEPTH);
        const int num_frames_to_print =
            stack_frames_exceed_max ? CTB_MAX_CALL_STACK_DEPTH : num_frames;

        for (int i = 0; i < num_frames_to_print; i++)
        {
            safe_print_frame(i, &context->call_stack_frames[i]);
        }

        if (stack_frames_exceed_max)
        {
            safe_print_str("\n      [... Skipped ");
            safe_print_int(num_frames - CTB_MAX_CALL_STACK_DEPTH);
            safe_print_str(" frames ...]\n\n");
        }
    }

    /* Print Signal Error Message */
    safe_print_str(error_to_string(ctb_error));
    safe_print_str("\n");

    for (int i = 0; i < CTB_DEFAULT_TERMINAL_WIDTH; i++)
    {
        safe_print_str("-");
    }
    safe_print_str("\n");
}
