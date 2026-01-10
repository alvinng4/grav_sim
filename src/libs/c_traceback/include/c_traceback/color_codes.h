/**
 * \file color_codes.h
 * \brief ANSI color codes for C traceback and error messages.
 *
 * \author Ching-Yin NG
 */

#ifndef C_TRACEBACK_COLOR_CODES_H
#define C_TRACEBACK_COLOR_CODES_H

// clang-format off
#define CTB_RESET_COLOR "\033[0m"                           // Reset to normal text

/* Theme color */
#define CTB_THEME_COLOR "\033[0;36m"                        // Cyan
#define CTB_THEME_BOLD_COLOR "\033[1;36m"                   // Cyan, Bold

/* Text color */
#define CTB_TEXT_PRIMARY_COLOR ""                           // Text color (default)
#define CTB_TEXT_PRIMARY_BOLD_COLOR "\033[1m"               // Default text color, Bold
#define CTB_TEXT_SECONDARY_COLOR "\033[38;5;240m"           // Dark gray

/* Error and warning */
#define CTB_ERROR_BOLD_COLOR "\033[1;31m"                   // Bright Red, Bold
#define CTB_ERROR_COLOR "\033[0;31m"                        // Red
#define CTB_WARNING_BOLD_COLOR "\033[1;33m"                 // Bright Yellow, Bold
#define CTB_WARNING_COLOR "\033[0;33m"                      // Yellow

/* Normal text / messages */
#define CTB_NORMAL_BOLD_COLOR "\033[1;95m"                  // Bright Purple, Bold
#define CTB_NORMAL_COLOR "\033[0;35m"                       // Purple

/* Traceback */
#define CTB_TRACEBACK_TEXT_COLOR CTB_TEXT_SECONDARY_COLOR                   // Dark gray
#define CTB_TRACEBACK_COUNTER_COLOR CTB_TEXT_SECONDARY_COLOR                // Dark gray
#define CTB_TRACEBACK_ANOTHER_EXCEPTION_TEXT_COLOR CTB_TEXT_PRIMARY_COLOR   // White / text color

/* File / Line / Func */
#define CTB_TRACEBACK_FILE_PARENT_COLOR CTB_TRACEBACK_TEXT_COLOR    // Dark gray
#define CTB_TRACEBACK_FILE_COLOR CTB_THEME_BOLD_COLOR               // Cyan, Bold
#define CTB_TRACEBACK_LINE_COLOR CTB_THEME_BOLD_COLOR               // Cyan, Bold
#define CTB_TRACEBACK_FUNC_COLOR CTB_THEME_BOLD_COLOR               // Cyan, Bold
// clang-format on

#endif /* C_TRACEBACK_COLOR_CODES_H */
