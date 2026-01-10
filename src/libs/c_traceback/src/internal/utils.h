/**
 * \file utils.h
 * \brief Utility functions for c_traceback library.
 *
 * \author Ching-Yin Ng
 */

#ifndef C_TRACEBACK_INTERNAL_UTILS_H
#define C_TRACEBACK_INTERNAL_UTILS_H

#include <stdbool.h>

/**
 * \brief Determine if we should use UTF-8 encoding for the given output stream.
 *
 * \param[in] stream The output stream (e.g., stdout, stderr).
 * \return true if UTF-8 encoding should be used, false otherwise.
 */
bool should_use_utf8(FILE *stream);

/**
 * \brief Determine if ANSI color codes should be used for the given output stream.
 *
 * \param[in] stream The output stream (e.g., stdout, stderr).
 * \return true if ANSI color codes should be used, false otherwise.
 */
bool should_use_color(FILE *stream);

/**
 * \brief Get the length of the parent directory path from the given file path.
 *
 * \param[in] path The file path.
 * \return The length of the parent directory path. Returns 0 if there is no file
 * separator.
 */
int get_parent_path_length(const char *restrict path);

/**
 * \brief Get the width of the terminal for the given output stream.
 *
 * \param[in] stream The output stream (e.g., stdout, stderr), or NULL for files.
 * \return The width of the terminal in characters.
 */
int get_terminal_width(FILE *stream);

#endif /* C_TRACEBACK_INTERNAL_UTILS_H */
