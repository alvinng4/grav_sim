/**
 * \file signal_handler.c
 * \brief Signal handlers for C Traceback.
 *
 * \author Ching-Yin Ng
 */

#include <signal.h>
#include <stdio.h>
#include <stdlib.h>

#include "c_traceback.h"

#ifdef _WIN32

static void ctb_internal_signal_handler(int sig)
{
    CTB_Error ctb_sig = CTB_SIGNAL_ERROR;

    // clang-format off
    switch (sig)
    {
        case SIGABRT: ctb_sig = CTB_SIGNAL_ABORT; break;
        case SIGSEGV: ctb_sig = CTB_SIGNAL_SEGMENTATION_FAULT; break;
        case SIGILL:  ctb_sig = CTB_SIGNAL_INVALID_INSTRUCTION; break;
        case SIGTERM: ctb_sig = CTB_SIGNAL_TERMINATION; break;
        case SIGFPE:  ctb_sig = CTB_SIGNAL_FLOATING_POINT_EXCEPTION; break;
        case SIGINT:  ctb_sig = CTB_SIGNAL_KEYBOARD_INTERRUPT; break;
        default:      ctb_sig = CTB_SIGNAL_ERROR; break;
    }
    // clang-format on

    ctb_dump_traceback_signal(ctb_sig);

    signal(sig, SIG_DFL);
    raise(sig);
}

void ctb_install_signal_handler(void)
{
    if (signal(SIGABRT, ctb_internal_signal_handler) == SIG_ERR)
    {
        LOG_WARNING(CTB_WARNING, "Failed to set signal handler for SIGABRT");
    }
    if (signal(SIGSEGV, ctb_internal_signal_handler) == SIG_ERR)
    {
        LOG_WARNING(CTB_WARNING, "Failed to set signal handler for SIGSEGV");
    }
    if (signal(SIGILL, ctb_internal_signal_handler) == SIG_ERR)
    {
        LOG_WARNING(CTB_WARNING, "Failed to set signal handler for SIGILL");
    }
    if (signal(SIGTERM, ctb_internal_signal_handler) == SIG_ERR)
    {
        LOG_WARNING(CTB_WARNING, "Failed to set signal handler for SIGTERM");
    }
    if (signal(SIGFPE, ctb_internal_signal_handler) == SIG_ERR)
    {
        LOG_WARNING(CTB_WARNING, "Failed to set signal handler for SIGFPE");
    }
    if (signal(SIGINT, ctb_internal_signal_handler) == SIG_ERR)
    {
        LOG_WARNING(CTB_WARNING, "Failed to set signal handler for SIGINT");
    }
}

#else
#include <unistd.h>

#define CTB_ALT_STACK_SIZE (64 * 1024)
static char ctb_stack_memory[CTB_ALT_STACK_SIZE];

static void ctb_internal_signal_handler(int sig, siginfo_t *info, void *context)
{
    (void)context;
    (void)info;

    CTB_Error ctb_sig = CTB_SIGNAL_ERROR;

    // clang-format off
    switch (sig)
    {
        case SIGABRT: ctb_sig = CTB_SIGNAL_ABORT; break;
        case SIGSEGV: ctb_sig = CTB_SIGNAL_SEGMENTATION_FAULT; break;
        case SIGILL:  ctb_sig = CTB_SIGNAL_INVALID_INSTRUCTION; break;
        case SIGTERM: ctb_sig = CTB_SIGNAL_TERMINATION; break;
        case SIGFPE:  ctb_sig = CTB_SIGNAL_FLOATING_POINT_EXCEPTION; break;
        case SIGINT:  ctb_sig = CTB_SIGNAL_KEYBOARD_INTERRUPT; break;
#ifdef SIGBUS
        case SIGBUS:  ctb_sig = CTB_SIGNAL_SEGMENTATION_FAULT; break;
#endif /* SIGBUS */
        default:      ctb_sig = CTB_SIGNAL_ERROR; break;
    }
    // clang-format on

    ctb_dump_traceback_signal(ctb_sig);

    /* Restore default handler and re-raise signal */
    struct sigaction sa_dfl;
    sa_dfl.sa_handler = SIG_DFL;
    sigemptyset(&sa_dfl.sa_mask);
    sa_dfl.sa_flags = 0;
    sigaction(sig, &sa_dfl, NULL);
    raise(sig);
}

void ctb_install_signal_handler(void)
{
    stack_t ss;
    ss.ss_sp = ctb_stack_memory;
    ss.ss_size = CTB_ALT_STACK_SIZE;
    ss.ss_flags = 0;

    if (sigaltstack(&ss, NULL) == -1)
    {
        LOG_WARNING(CTB_WARNING, "Failed to set alternate signal stack");
    }

    struct sigaction sa;

    sa.sa_flags = SA_SIGINFO | SA_ONSTACK;
    sa.sa_sigaction = ctb_internal_signal_handler;
    sigfillset(&sa.sa_mask);

    // Register Signals
    if (sigaction(SIGABRT, &sa, NULL) == -1)
    {
        LOG_WARNING(CTB_WARNING, "Failed to set signal handler for SIGABRT");
    }
    if (sigaction(SIGSEGV, &sa, NULL) == -1)
    {
        LOG_WARNING(CTB_WARNING, "Failed to set signal handler for SIGSEGV");
    }
    if (sigaction(SIGILL, &sa, NULL) == -1)
    {
        LOG_WARNING(CTB_WARNING, "Failed to set signal handler for SIGILL");
    }
    if (sigaction(SIGTERM, &sa, NULL) == -1)
    {
        LOG_WARNING(CTB_WARNING, "Failed to set signal handler for SIGTERM");
    }
    if (sigaction(SIGFPE, &sa, NULL) == -1)
    {
        LOG_WARNING(CTB_WARNING, "Failed to set signal handler for SIGFPE");
    }
    if (sigaction(SIGINT, &sa, NULL) == -1)
    {
        LOG_WARNING(CTB_WARNING, "Failed to set signal handler for SIGINT");
    }
#ifdef SIGBUS
    if (sigaction(SIGBUS, &sa, NULL) == -1)
    {
        LOG_WARNING(CTB_WARNING, "Failed to set signal handler for SIGBUS");
    }
#endif /* SIGBUS */
}

#endif /* _WIN32 */
