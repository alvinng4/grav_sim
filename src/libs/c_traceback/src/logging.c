/**
 * \file logging.c
 * \brief Function definitions for simple logging.
 *
 * \author Ching-Yin Ng
 */

#include <stdarg.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "c_traceback.h"
#include "internal/utils.h"

/**
 * \brief Helper for logging simple messages without the message body.
 *
 * \param[in] use_color Whether to use color in the output.
 * \param[in, out] stream The output stream.
 * \param[in] header_color The color code for the header.
 * \param[in] message_color The color code for the message.
 * \param[in] file_address The file address.
 * \param[in] line The line number.
 * \param[in] func The function name.
 * \param[in] header The header string.
 */
static void ctb_logging_core(
    const bool use_color,
    FILE *stream,
    const char *header_color,
    const char *message_color,
    const char *restrict file_address,
    const int line,
    const char *restrict func,
    const char *restrict header
)
{
    if (use_color)
    {
        const int dir_len = get_parent_path_length(file_address);

        // clang-format off
        fprintf(
            stream,
            "%s%s:%s %sFile \"%s",
            header_color, header, CTB_RESET_COLOR,
            CTB_TRACEBACK_TEXT_COLOR, CTB_RESET_COLOR
        );
        // clang-format on

        /* Print file address */
        if (dir_len > 0)
        {
            fprintf(
                stream,
                "%s%.*s%s",
                CTB_TRACEBACK_TEXT_COLOR,
                dir_len,
                file_address,
                CTB_RESET_COLOR
            );
            fprintf(
                stream,
                "%s%s%s",
                CTB_TRACEBACK_FILE_COLOR,
                file_address + dir_len,
                CTB_RESET_COLOR
            );
        }
        else
        {
            fprintf(
                stream,
                "%s%s%s",
                CTB_TRACEBACK_FILE_COLOR,
                file_address,
                CTB_RESET_COLOR
            );
        }

        // clang-format off
        fprintf(
            stream,
            "%s\", line%s %s%d%s %sin%s %s%s%s:\n   %s",
            CTB_TRACEBACK_TEXT_COLOR, CTB_RESET_COLOR,
            CTB_TRACEBACK_LINE_COLOR, line, CTB_RESET_COLOR,
            CTB_TRACEBACK_TEXT_COLOR, CTB_RESET_COLOR,
            CTB_TRACEBACK_FUNC_COLOR, func, CTB_RESET_COLOR,
            message_color
        );
        // clang-format on
    }
    else
    {
        fprintf(
            stream,
            "%s: File \"%s\", line %d in %s:\n    ",
            header,
            file_address,
            line,
            func
        );
    }
}

/**
 * \brief Log message.
 *
 * \param[in, out] stream The output stream.
 * \param[in] header_color The color code for the header.
 * \param[in] message_color The color code for the message.
 * \param[in] file_address The file address.
 * \param[in] line The line number.
 * \param[in] func The function name.
 * \param[in] header The header string.
 * \param[in] msg The message string.
 */
static void ctb_log(
    FILE *stream,
    const char *header_color,
    const char *message_color,
    const char *restrict file_address,
    const int line,
    const char *restrict func,
    const char *restrict header,
    const char *restrict msg
)
{
    const bool use_color = should_use_color(stream);
    ctb_logging_core(
        use_color, stream, header_color, message_color, file_address, line, func, header
    );

    if (use_color)
    {
        fprintf(stream, "%s%s\n", msg, CTB_RESET_COLOR);
    }
    else
    {
        fprintf(stream, "%s\n", msg);
    }

    fflush(stream);
}

/**
 * \brief Log message with variadic arguments.
 *
 * \param[in, out] stream The output stream.
 * \param[in] header_color The color code for the header.
 * \param[in] message_color The color code for the message.
 * \param[in] file_address The file address.
 * \param[in] line The line number.
 * \param[in] func The function name.
 * \param[in] header The header string.
 * \param[in] msg The message format string.
 * \param[in] args The variadic arguments.
 */
static void ctb_log_fmt(
    FILE *stream,
    const char *header_color,
    const char *message_color,
    const char *restrict file_address,
    const int line,
    const char *restrict func,
    const char *restrict header,
    const char *restrict msg,
    va_list args
)
{
    const bool use_color = should_use_color(stream);
    ctb_logging_core(
        use_color, stream, header_color, message_color, file_address, line, func, header
    );

    if (use_color)
    {
        vfprintf(stream, msg, args);
        fprintf(stream, "%s\n", CTB_RESET_COLOR);
    }
    else
    {
        vfprintf(stream, msg, args);
        fputc('\n', stream);
    }
    fflush(stream);
}

void ctb_log_error(
    const char *restrict file,
    const int line,
    const char *restrict func,
    const CTB_Error error,
    const char *restrict msg
)
{
    FILE *stream = stderr;
    ctb_log(
        stream,
        CTB_ERROR_BOLD_COLOR,
        CTB_ERROR_COLOR,
        file,
        line,
        func,
        error_to_string(error),
        msg
    );
}

void ctb_log_warning(
    const char *restrict file,
    const int line,
    const char *restrict func,
    const CTB_Warning warning,
    const char *restrict msg
)
{
    FILE *stream = stderr;
    ctb_log(
        stream,
        CTB_WARNING_BOLD_COLOR,
        CTB_WARNING_COLOR,
        file,
        line,
        func,
        warning_to_string(warning),
        msg
    );
}

void ctb_log_message(
    const char *restrict file,
    const int line,
    const char *restrict func,
    const char *restrict msg
)
{
    FILE *stream = stdout;
    ctb_log(
        stream,
        CTB_NORMAL_BOLD_COLOR,
        CTB_NORMAL_COLOR,
        file,
        line,
        func,
        "Message",
        msg
    );
}

void ctb_log_error_fmt(
    const char *restrict file,
    const int line,
    const char *restrict func,
    const CTB_Error error,
    const char *restrict msg,
    ...
)
{
    FILE *stream = stderr;
    va_list args;
    va_start(args, msg);
    ctb_log_fmt(
        stream,
        CTB_ERROR_BOLD_COLOR,
        CTB_ERROR_COLOR,
        file,
        line,
        func,
        error_to_string(error),
        msg,
        args
    );
    va_end(args);
}

void ctb_log_warning_fmt(
    const char *restrict file,
    const int line,
    const char *restrict func,
    const CTB_Warning warning,
    const char *restrict msg,
    ...
)
{
    FILE *stream = stderr;
    va_list args;
    va_start(args, msg);
    ctb_log_fmt(
        stream,
        CTB_WARNING_BOLD_COLOR,
        CTB_WARNING_COLOR,
        file,
        line,
        func,
        warning_to_string(warning),
        msg,
        args
    );
    va_end(args);
}

void ctb_log_message_fmt(
    const char *restrict file,
    const int line,
    const char *restrict func,
    const char *restrict msg,
    ...
)
{
    FILE *stream = stdout;
    va_list args;
    va_start(args, msg);
    ctb_log_fmt(
        stream,
        CTB_NORMAL_BOLD_COLOR,
        CTB_NORMAL_COLOR,
        file,
        line,
        func,
        "Message",
        msg,
        args
    );
    va_end(args);
}
