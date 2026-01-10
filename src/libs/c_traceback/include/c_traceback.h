/**
 * \file c_traceback.h
 * \brief Library header for c_traceback library
 *
 * \author Ching-Yin Ng
 */

#ifndef C_TRACEBACK_H
#define C_TRACEBACK_H

#include "c_traceback/color_codes.h"
#include "c_traceback/error.h"
#include "c_traceback/error_codes.h"
#include "c_traceback/logging.h"
#include "c_traceback/signal_handler.h"
#include "c_traceback/trace.h"
#include "c_traceback/traceback.h"

#ifndef CTB_VERSION
#define CTB_VERSION "Unknown"
#endif

/**
 * Traceback header
 * (i.e. the title printed at the top of the traceback output
 * before "(most recent call last)"). It will be highlighted in
 * bold.
 *
 * You can change this value to customize the header, e.g.
 * "MyApp Traceback"
 *
 * If its value is NULL or an empty string, "Traceback" will be used.
 */
#define CTB_TRACEBACK_HEADER "grav_sim Traceback"

// Maximum number of call stack frames
#define CTB_MAX_CALL_STACK_DEPTH 32

// Maximum number of simultaneous errors
#define CTB_MAX_NUM_ERROR 16

// Maximum length of error message
#define CTB_MAX_ERROR_MESSAGE_LENGTH 256

// Terminal width when it cannot be determined
#define CTB_DEFAULT_TERMINAL_WIDTH 80

// Default output width when printing to file
#define CTB_DEFAULT_FILE_WIDTH 120

// Horizontal rule width
#define CTB_HRULE_MAX_WIDTH 120
#define CTB_HRULE_MIN_WIDTH 50

/**
 * \brief Print compilation information such as compiler version,
 *        compilation date, and configurations.
 */
void ctb_print_compilation_info(void);

#endif /* C_TRACEBACK_H */
