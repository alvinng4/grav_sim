/**
 * \file error.c
 * \brief Implementation of exceptions raising and handling functions.
 *
 * \author Ching-Yin Ng
 */

#include <stdarg.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>

#include "internal/trace.h"

/**
 * \brief Setup error snapshot core information without message.
 *
 * \param[in,out] context The CTB_Context pointer.
 * \param[in,out] error_snapshot The CTB_Error_Snapshot_ pointer to setup.
 * \param[in] error The error type.
 * \param[in] file File where the error is throwd.
 * \param[in] line Line number where the error is throwd.
 * \param[in] func Function name where the error is throwd.
 */
static void ctb_setup_error_snapshot_core(
    CTB_Context *restrict context,
    CTB_Error_Snapshot_ *restrict error_snapshot,
    CTB_Error error,
    const char *restrict file,
    const int line,
    const char *restrict func
)
{
    error_snapshot->error = error;
    error_snapshot->call_depth = context->call_depth;
    error_snapshot->error_frame.filename = file;
    error_snapshot->error_frame.line_number = line;
    error_snapshot->error_frame.function_name = func;
    error_snapshot->error_frame.source_code = "<Error thrown here>";

    const int min_depth = (context->call_depth < CTB_MAX_CALL_STACK_DEPTH)
                              ? context->call_depth
                              : CTB_MAX_CALL_STACK_DEPTH;
    if (min_depth > 0)
    {
        memcpy(
            error_snapshot->call_stack_frames,
            context->call_stack_frames,
            sizeof(CTB_Frame_) * min_depth
        );
    }
}

void ctb_throw_error(
    CTB_Error error,
    const char *restrict file,
    const int line,
    const char *restrict func,
    const char *restrict msg
)
{
    CTB_Context *context = get_context();

    const int num_errors = context->num_errors;
    if (num_errors >= 0 && num_errors < CTB_MAX_NUM_ERROR)
    {
        CTB_Error_Snapshot_ *error_snapshot = &(context->error_snapshots[num_errors]);
        ctb_setup_error_snapshot_core(context, error_snapshot, error, file, line, func);
        if (msg != NULL)
        {
            snprintf(
                error_snapshot->error_message, CTB_MAX_ERROR_MESSAGE_LENGTH, "%s", msg
            );
        }
        else
        {
            error_snapshot->error_message[0] = '\0';
        }
    }

    (context->num_errors)++;
}

void ctb_throw_error_fmt(
    CTB_Error error,
    const char *restrict file,
    const int line,
    const char *restrict func,
    const char *restrict msg,
    ...
)
{
    CTB_Context *context = get_context();

    const int num_errors = context->num_errors;
    if (num_errors >= 0 && num_errors < CTB_MAX_NUM_ERROR)
    {
        CTB_Error_Snapshot_ *error_snapshot = &(context->error_snapshots[num_errors]);
        ctb_setup_error_snapshot_core(context, error_snapshot, error, file, line, func);

        va_list args;
        va_start(args, msg);
        vsnprintf(
            error_snapshot->error_message, CTB_MAX_ERROR_MESSAGE_LENGTH, msg, args
        );
        va_end(args);
    }

    (context->num_errors)++;
}

bool ctb_check_error(void)
{
    return get_context()->num_errors > 0;
}

void ctb_clear_error(void)
{
    get_context()->num_errors = 0;
}

void ctb_clear_context(void)
{
    CTB_Context *context = get_context();
    context->num_errors = 0;
    context->call_depth = 0;
}
