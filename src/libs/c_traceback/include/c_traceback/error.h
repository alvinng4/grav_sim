/**
 * \file error.c
 * \brief Header file for exceptions raising and handling functions.
 *
 * \author Ching-Yin Ng
 */

#ifndef C_TRACEBACK_ERROR_H
#define C_TRACEBACK_ERROR_H

#include <stdbool.h>

#include "c_traceback/error_codes.h"

/**
 * \brief Wrapper for throwing an error with the current call stack.
 *
 * \param[in] ctb_error The error type.
 * \param[in] msg Error message.
 */
#define THROW(ctb_error, msg)                                                          \
    do                                                                                 \
    {                                                                                  \
        ctb_throw_error(ctb_error, __FILE__, __LINE__, __func__, msg);                 \
    } while (0)

/**
 * \brief Wrapper for throwing an error with formatted message with the current call
 * stack.
 *
 * \param[in] ctb_error The error type.
 * \param[in] msg Error message.
 * \param[in] ... Additional arguments for formatting the message.
 */
#define THROW_FMT(ctb_error, msg, ...)                                                 \
    do                                                                                 \
    {                                                                                  \
        ctb_throw_error_fmt(                                                           \
            ctb_error, __FILE__, __LINE__, __func__, msg, __VA_ARGS__                  \
        );                                                                             \
    } while (0)

/**
 * \brief Throw an error with the current call stack.
 *
 * \param[in] error The error type.
 * \param[in] file File where the message is sent.
 * \param[in] line Line number where the message is sent.
 * \param[in] func Function where the message is sent.
 * \param[in] msg Error message.
 * \param[in] ... Additional arguments for formatting the message.
 */
void ctb_throw_error(
    CTB_Error error,
    const char *restrict file,
    const int line,
    const char *restrict func,
    const char *restrict msg
);

/**
 * \brief Throw an error with formatted message with the current call stack.
 *
 * \param[in] error The error type.
 * \param[in] file File where the message is sent.
 * \param[in] line Line number where the message is sent.
 * \param[in] func Function where the message is sent.
 * \param[in] msg Error message.
 * \param[in] ... Additional arguments for formatting the message.
 */
void ctb_throw_error_fmt(
    CTB_Error error,
    const char *restrict file,
    const int line,
    const char *restrict func,
    const char *restrict msg,
    ...
);

/**
 * \brief Check if any error has occurred.
 *
 * \return true if an error has occurred, false otherwise.
 */
bool ctb_check_error(void);

/**
 * \brief Clear all recorded errors.
 */
void ctb_clear_error(void);

/**
 * \brief Clear the C Traceback context (including errors and call stack).
 */
void ctb_clear_context(void);

#endif // C_TRACEBACK_ERROR_H
