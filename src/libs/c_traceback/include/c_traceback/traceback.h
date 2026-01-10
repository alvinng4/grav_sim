/**
 * \file traceback.c
 * \brief Header file for traceback-related macros and functions.
 *
 * \author Ching-Yin Ng
 */

#ifndef C_TRACEBACK_TRACEBACK_H
#define C_TRACEBACK_TRACEBACK_H

/**
 * \brief Log the current call stack to stderr.
 */
void ctb_log_call_stack(void);

/**
 * \brief Log the traceback of all recorded errors to stderr.
 */
void ctb_log_traceback(void);

/**
 * \brief Dump the traceback of all recorded errors to stderr and clear errors.
 */
void ctb_dump_traceback(void);

#endif /* C_TRACEBACK_TRACEBACK_H */
