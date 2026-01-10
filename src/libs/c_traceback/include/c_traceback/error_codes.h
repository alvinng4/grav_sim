/**
 * \file c_traceback_errors.h
 * \brief C Traceback Error and Warning Hierarchy
 *
 * \author Ching-Yin Ng
 */

#ifndef C_TRACEBACK_ERROR_CODES_H
#define C_TRACEBACK_ERROR_CODES_H

/**
    C Traceback Error Hierarchy:

    Error
    ├── SystemError
    ├── ResourceError
    ├── ValueError
    ├── SignalError
    │    ├── SignalAbort
    │    ├── SignalSegmentationFault
    │    ├── SignalInvalidInstruction
    │    ├── SignalTermination
    │    ├── SignalFloatingPointException
    │    └── SignalKeyboardInterrupt
    ├── MemoryError
    │    ├── OutOfMemoryError
    │    └── BufferError
    ├── PointerError
    │    └── NullPointerError
    ├── RuntimeError
    │    └── NotImplementedError
    ├── LookupError
    │    ├── IndexError
    │    └── KeyError
    ├── MathError
    │    ├── MathDomainError
    │    ├── MathOverflowError
    │    ├── MathUnderflowError
    │    └── ZeroDivisionError
    ├── OSError
    │    ├── BlockingIOError
    │    ├── ChildProcessError
    │    ├── FileExistsError
    │    ├── FileNotFoundError
    │    ├── InterruptedError
    │    ├── IsADirectoryError
    │    ├── NotADirectoryError
    │    ├── PermissionError
    │    ├── ProcessLookupError
    │    └── TimeoutError
    ├── NetworkError
    │    ├── ConnectionFailedError
    │    └── HostUnreachableError
    ├── ConcurrencyError
    │    └── ThreadCreationError
    ├── ParseError
    │    ├── EncodingError
    │    ├── InvalidFormatError
    │    ├── SyntaxError
    │    ├── UnexpectedEOFError
    │    └── ChecksumMismatchError
    └── DynamicLinkError
         └── LibraryLoadError

    C Traceback Warning Hierarchy:

    Warning
    ├── DeprecationWarning
    ├── PerformanceWarning
    ├── ResourceWarning
    ├── PrecisionWarning
    ├── CompatibilityWarning
    ├── SecurityWarning
    ├── ParseWarning
    └── UserWarning

*/

typedef enum CTB_Error
{

    CTB_SUCCESS = 0,
    CTB_UNKNOWN_ERROR = -1,

    /* --- Base Error & System Errors (1-49) --- */
    CTB_ERROR = 1,
    CTB_SYSTEM_ERROR,
    CTB_RESOURCE_ERROR,
    CTB_VALUE_ERROR,

    /* --- Signal Errors (50–99) --- */
    CTB_SIGNAL_ERROR = 50,
    CTB_SIGNAL_ABORT,
    CTB_SIGNAL_SEGMENTATION_FAULT,
    CTB_SIGNAL_INVALID_INSTRUCTION,
    CTB_SIGNAL_TERMINATION,
    CTB_SIGNAL_FLOATING_POINT_EXCEPTION,
    CTB_SIGNAL_KEYBOARD_INTERRUPT,

    /* --- Memory Errors (100–199) --- */
    CTB_MEMORY_ERROR = 100,
    CTB_OUT_OF_MEMORY_ERROR,
    CTB_BUFFER_ERROR,

    /* --- Pointer Errors (200–249) --- */
    CTB_POINTER_ERROR = 200,
    CTB_NULL_POINTER_ERROR,

    /* --- Runtime Errors (250–299) --- */
    CTB_RUNTIME_ERROR = 250,
    CTB_NOT_IMPLEMENTED_ERROR,

    /* --- Lookup Errors (300–399) --- */
    CTB_LOOKUP_ERROR = 300,
    CTB_INDEX_ERROR,
    CTB_KEY_ERROR,

    /* --- Math Errors (400–499) --- */
    CTB_MATH_ERROR = 400,
    CTB_MATH_DOMAIN_ERROR,
    CTB_MATH_OVERFLOW_ERROR,
    CTB_MATH_UNDERFLOW_ERROR,
    CTB_ZERO_DIVISION_ERROR,

    /* --- OS Errors (500–599) --- */
    CTB_OS_ERROR = 500,
    CTB_BLOCKING_IO_ERROR,
    CTB_CHILD_PROCESS_ERROR,
    CTB_FILE_EXISTS_ERROR,
    CTB_FILE_NOT_FOUND_ERROR,
    CTB_INTERRUPTED_ERROR,
    CTB_IS_DIRECTORY_ERROR,
    CTB_NOT_DIRECTORY_ERROR,
    CTB_PERMISSION_ERROR,
    CTB_PROCESS_LOOKUP_ERROR,
    CTB_TIMEOUT_ERROR,

    /* --- Network Errors (600–699) --- */
    CTB_NETWORK_ERROR = 600,
    CTB_CONNECTION_FAILED_ERROR,
    CTB_HOST_UNREACHABLE_ERROR,

    /* --- Concurrency Errors (700–799) --- */
    CTB_CONCURRENCY_ERROR = 700,
    CTB_THREAD_CREATION_ERROR,

    /* --- Parse Errors (800–899) --- */
    CTB_PARSE_ERROR = 800,
    CTB_ENCODING_ERROR,
    CTB_INVALID_FORMAT_ERROR,
    CTB_SYNTAX_ERROR,
    CTB_UNEXPECTED_EOF_ERROR,
    CTB_CHECKSUM_MISMATCH_ERROR,

    /* --- Dynamic Link Errors (900–999) --- */
    CTB_DYNAMIC_LINK_ERROR = 900,
    CTB_LIBRARY_LOAD_ERROR,

    /* --- 1000+: application-specific --- */

} CTB_Error;

typedef enum CTB_Warning
{
    CTB_UNKNOWN_WARNING = -1,
    CTB_NO_WARNING = 0,

    CTB_WARNING,
    CTB_DEPRECATION_WARNING,
    CTB_PERFORMANCE_WARNING,
    CTB_RESOURCE_WARNING,
    CTB_PRECISION_WARNING,
    CTB_COMPATIBILITY_WARNING,
    CTB_SECURITY_WARNING,
    CTB_PARSE_WARNING,
    CTB_USER_WARNING
} CTB_Warning;

/**
 * \brief Convert error type to its corresponding string representation.
 *
 * \param[in] error The error type to convert.
 * \return A constant character pointer to the string representation of the error code.
 */
const char *error_to_string(CTB_Error error);

/**
 * \brief Convert warning type to its corresponding string representation.
 *
 * \param[in] warning The warning type to convert.
 * \return A constant character pointer to the string representation of the error code.
 */
const char *warning_to_string(CTB_Warning warning);

#endif /* C_TRACEBACK_ERROR_CODES_H */
