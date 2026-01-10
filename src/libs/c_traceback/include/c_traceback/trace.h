/**
 * \file trace.h
 * \brief Header file for call stack tracing macros and functions.
 *
 * \author Ching-Yin Ng
 */

#ifndef C_TRACEBACK_TRACE_H
#define C_TRACEBACK_TRACE_H

/**
 * \brief Wrapper macro for expression to automatically manage call stack frames without
 * checking for errors.
 *
 * \param[in] expr The expression to be traced.
 */
#define TRACE(expr)                                                                    \
    do                                                                                 \
    {                                                                                  \
        ctb_push_call_stack_frame(__FILE__, __func__, __LINE__, #expr);                \
        (expr);                                                                        \
        ctb_pop_call_stack_frame();                                                    \
    } while (0)

/**
 * \brief Wrapper macro for a code block to automatically manage call stack frames
 * without checking for errors.
 *
 * \param[in] ... The block of code to be traced.
 */
#define TRACE_BLOCK(...)                                                               \
    do                                                                                 \
    {                                                                                  \
        ctb_push_call_stack_frame(__FILE__, __func__, __LINE__, #__VA_ARGS__);         \
        __VA_ARGS__                                                                    \
        ctb_pop_call_stack_frame();                                                    \
    } while (0)

/**
 * \brief Wrapper macro for expression to automatically manage call stack frames and
 * check for errors.
 *
 * \param[in] expr The expression to be traced.
 *
 * \return Whether the expression executed without error.
 */
#define TRY(expr)                                                                      \
    (ctb_push_call_stack_frame(__FILE__, __func__, __LINE__, #expr),                   \
     (expr),                                                                           \
     ctb_pop_call_stack_frame(),                                                       \
     !ctb_check_error())

/**
 * \brief Wrapper for an expression. If an error occurs after the expression, jump to
 * label.
 *
 * \param[in] expr The expression to be traced.
 * \param[in] label The label to jump to on error.
 */
#define TRY_GOTO(expr, label)                                                          \
    do                                                                                 \
    {                                                                                  \
        ctb_push_call_stack_frame(__FILE__, __func__, __LINE__, #expr);                \
        (expr);                                                                        \
        ctb_pop_call_stack_frame();                                                    \
        if (ctb_check_error())                                                         \
        {                                                                              \
            goto label;                                                                \
        }                                                                              \
    } while (0)

/**
 * \brief Wrapper for a block of code. If an error occurs after the block executes, jump
 * to label.
 *
 * \param[in] label The label to jump to on error.
 * \param[in] ...   The block of code to be traced.
 */
#define TRY_BLOCK_GOTO(label, ...)                                                     \
    do                                                                                 \
    {                                                                                  \
        ctb_push_call_stack_frame(__FILE__, __func__, __LINE__, #__VA_ARGS__);         \
        __VA_ARGS__                                                                    \
        ctb_pop_call_stack_frame();                                                    \
        if (ctb_check_error())                                                         \
        {                                                                              \
            goto label;                                                                \
        }                                                                              \
    } while (0)

/**
 * \brief Push a new call stack frame.
 * \param[in] file File where the function is called.
 * \param[in] func Function name where the function is called.
 * \param[in] line Line number where the function is called.
 * \param[in] source_code Source code of the function call.
 */
void ctb_push_call_stack_frame(
    const char *file, const char *func, const int line, const char *source_code
);

/**
 * \brief Pop the top call stack frame.
 */
void ctb_pop_call_stack_frame(void);

#endif /* C_TRACEBACK_TRACE_H */
