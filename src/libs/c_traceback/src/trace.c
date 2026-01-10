/**
 * \file context.c
 * \brief Declaration of C Traceback global context, and function definitions for C
 * Traceback library.
 *
 * \author Ching-Yin Ng
 */

#include "internal/trace.h"

static ctb_thread_local CTB_Context ctb_traceback_context = {0};

CTB_Context *get_context(void)
{
    return &ctb_traceback_context;
}

void ctb_push_call_stack_frame(
    const char *file, const char *func, const int line, const char *source_code
)
{
    CTB_Context *context = get_context();
    const int call_depth = context->call_depth;
    int frame_index = call_depth;

    if (call_depth < 0)
    {
        return;
    }

    if (call_depth >= CTB_MAX_CALL_STACK_DEPTH)
    {
        frame_index = CTB_MAX_CALL_STACK_DEPTH - 1;
    }

    CTB_Frame_ *frame = &context->call_stack_frames[frame_index];
    frame->filename = file;
    frame->function_name = func;
    frame->line_number = line;
    frame->source_code = source_code;
    context->call_depth++;
}

void ctb_pop_call_stack_frame(void)
{
    CTB_Context *context = get_context();
    if (context->call_depth > 0)
    {
        (context->call_depth)--;
    }
}
