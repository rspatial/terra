/*
 * Check for user interruption in C++ code interrupting execution of the current code
 *
 * This code has been written by Simon Urbanek
 * I took it from the R-devel mailing list
 * in the thread "[Rd] Interrupting C++ code execution"
 * The mail is dated: 25 april 2011
 *
 * It allows to check for user interruption without
 * leaving the c++ function that calls it.
 *
 * Potential drawbacks according to its author:
 * The problem with it is that it will eat all errors, even if they were not yours
 * (e.g. those resulting from events triggered the event loop), so I would not recommend it for general use.
 *
 */
#ifndef _RcppProgress_INTERRUPTS_HPP
#define _RcppProgress_INTERRUPTS_HPP

// N.B: seems the right way to include Rinternals.h in a Rcpp package
#include <Rcpp.h>


static void chkIntFn(void *dummy) {
	R_CheckUserInterrupt();
}

// this will call the above in a top-level context so it won't longjmp-out of your context
inline bool checkInterrupt() {
	return (R_ToplevelExec(chkIntFn, NULL) == FALSE);
}

// fix a bug because of the length macro (in sglOptim_1.0.122.1)
#undef length
#endif
