/**
 * \file signal_handler.h
 * \brief Header for signal handlers in C Traceback.
 *
 * \author Ching-Yin Ng
 */

#ifndef C_TRACEBACK_SIGNAL_HANDLER_H
#define C_TRACEBACK_SIGNAL_HANDLER_H

/**
 * \brief Install signal handlers for C Traceback.
 */
void ctb_install_signal_handler(void);

/**
 * \brief Dump the traceback to stderr on signal error.
 *
 * \param[in] ctb_error The signal error type.
 */
void ctb_dump_traceback_signal(const CTB_Error ctb_error);

#endif /* C_TRACEBACK_SIGNAL_HANDLER_H */
