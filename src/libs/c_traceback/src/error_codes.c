/**
 * \file error_codes.c
 * \brief Implementation of functions related to C Traceback Error Codes
 *
 * \author Ching-Yin Ng
 */

#include "c_traceback/error_codes.h"

// clang-format off
const char *error_to_string(CTB_Error error)
{
    switch (error)
    {
        case CTB_SUCCESS:                       return "Success";
        case CTB_UNKNOWN_ERROR:                 return "Unknown Error";
        
        /* Base */
        case CTB_ERROR:                         return "Error";
        case CTB_SYSTEM_ERROR:                  return "SystemError";
        case CTB_RESOURCE_ERROR:                return "ResourceError";
        case CTB_VALUE_ERROR:                   return "ValueError";

        /* Signal */
        case CTB_SIGNAL_ERROR:                  return "SignalError";
        case CTB_SIGNAL_ABORT:                  return "SignalAbort";
        case CTB_SIGNAL_SEGMENTATION_FAULT:     return "SignalSegmentationFault";
        case CTB_SIGNAL_INVALID_INSTRUCTION:    return "SignalInvalidInstruction";
        case CTB_SIGNAL_TERMINATION:            return "SignalTermination";
        case CTB_SIGNAL_FLOATING_POINT_EXCEPTION: return "SignalFloatingPointException";
        case CTB_SIGNAL_KEYBOARD_INTERRUPT:     return "SignalKeyboardInterrupt";

        /* Memory */
        case CTB_MEMORY_ERROR:                  return "MemoryError";
        case CTB_OUT_OF_MEMORY_ERROR:           return "OutOfMemoryError";
        case CTB_BUFFER_ERROR:                  return "BufferError";

        /* Pointer */
        case CTB_POINTER_ERROR:                 return "PointerError";
        case CTB_NULL_POINTER_ERROR:            return "NullPointerError";

        /* Runtime */
        case CTB_RUNTIME_ERROR:                 return "RuntimeError";
        case CTB_NOT_IMPLEMENTED_ERROR:         return "NotImplementedError";

        /* Lookup */
        case CTB_LOOKUP_ERROR:                  return "LookupError";
        case CTB_INDEX_ERROR:                   return "IndexError";
        case CTB_KEY_ERROR:                     return "KeyError";
        
        /* Math */
        case CTB_MATH_ERROR:                    return "MathError";
        case CTB_MATH_DOMAIN_ERROR:             return "MathDomainError";
        case CTB_MATH_OVERFLOW_ERROR:           return "MathOverflowError";
        case CTB_MATH_UNDERFLOW_ERROR:          return "MathUnderflowError";
        case CTB_ZERO_DIVISION_ERROR:           return "ZeroDivisionError";

        /* OS */
        case CTB_OS_ERROR:                      return "OSError";
        case CTB_BLOCKING_IO_ERROR:             return "BlockingIOError";
        case CTB_CHILD_PROCESS_ERROR:           return "ChildProcessError";
        case CTB_FILE_EXISTS_ERROR:             return "FileExistsError";
        case CTB_FILE_NOT_FOUND_ERROR:          return "FileNotFoundError";
        case CTB_INTERRUPTED_ERROR:             return "InterruptedError";
        case CTB_IS_DIRECTORY_ERROR:            return "IsDirectoryError";
        case CTB_NOT_DIRECTORY_ERROR:           return "NotDirectoryError";
        case CTB_PERMISSION_ERROR:              return "PermissionError";
        case CTB_PROCESS_LOOKUP_ERROR:          return "ProcessLookupError";
        case CTB_TIMEOUT_ERROR:                 return "TimeoutError";

        /* Network */
        case CTB_NETWORK_ERROR:                 return "NetworkError";
        case CTB_CONNECTION_FAILED_ERROR:       return "ConnectionFailedError";
        case CTB_HOST_UNREACHABLE_ERROR:        return "HostUnreachableError";

        /* Concurrency */
        case CTB_CONCURRENCY_ERROR:             return "ConcurrencyError";
        case CTB_THREAD_CREATION_ERROR:         return "ThreadCreationError";

        /* Parse */
        case CTB_PARSE_ERROR:                   return "ParseError";
        case CTB_ENCODING_ERROR:                return "EncodingError";
        case CTB_INVALID_FORMAT_ERROR:          return "InvalidFormatError";
        case CTB_SYNTAX_ERROR:                  return "SyntaxError";
        case CTB_UNEXPECTED_EOF_ERROR:          return "UnexpectedEOFError";
        case CTB_CHECKSUM_MISMATCH_ERROR:       return "ChecksumMismatchError";

        /* Dynamic Link */
        case CTB_DYNAMIC_LINK_ERROR:            return "DynamicLinkError";
        case CTB_LIBRARY_LOAD_ERROR:            return "LibraryLoadError";

        default:                                return "Unrecognized Error Code";
    }
}

const char *warning_to_string(CTB_Warning warning)
{
    switch (warning)
    {
        case CTB_NO_WARNING:            return "No Warning";
        case CTB_WARNING:               return "Warning";
        case CTB_DEPRECATION_WARNING:   return "DeprecationWarning";
        case CTB_PERFORMANCE_WARNING:   return "PerformanceWarning";
        case CTB_RESOURCE_WARNING:      return "ResourceWarning";
        case CTB_PRECISION_WARNING:     return "PrecisionWarning";
        case CTB_COMPATIBILITY_WARNING: return "CompatibilityWarning";
        case CTB_SECURITY_WARNING:      return "SecurityWarning";
        case CTB_PARSE_WARNING:         return "ParseWarning";
        case CTB_USER_WARNING:          return "UserWarning";
        case CTB_UNKNOWN_WARNING:
        default:                        return "Unknown Warning";
    }
}
// clang-format on
