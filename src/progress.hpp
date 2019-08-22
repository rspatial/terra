/*
 * progress.hpp
 *
 * A Front-end class for InterruptableProgressMonitor.
 *
 * Author: karl.forner@gmail.com
 *
 */
#ifndef _RcppProgress_PROGRESS_HPP
#define _RcppProgress_PROGRESS_HPP

#include "progress_monitor_interruptable.hpp"
#include "progress_bar_simple.hpp"

// e.g. for  Rf_error
#include <R_ext/Error.h>


class Progress {
public:
	/**
	 *
	 * Main constructor
	 * @param max the expected number of tasks to perform
	 * @param display_progress whether to display a progress bar in the console
	 * @param pb  the ProgressBar instance to use

	 */
	Progress(
	  unsigned long max,
	  bool display_progress = true,
    ProgressBar& pb = default_progress_bar()
  ) {
    //if ( monitor_singleton() != 0) { // something is wrong, two simultaneous Progress monitoring
      //Rf_error("ERROR: there is already an InterruptableProgressMonitor instance defined");
    //}
    monitor_singleton() = new InterruptableProgressMonitor(max, display_progress, pb);
	}

	~Progress() { cleanup(); 	}

public: // ==== USER INTERFACE =====
	/**
	 * cleanup
	 *
	 *  should normally not be called, unless a something bad happens (
	 *  a process/thread that crashes).
	 *
	 */
	void cleanup() {
    if (monitor_singleton() != 0) delete monitor_singleton();
    monitor_singleton() = 0;
	}

	/**
	 * increment the current progress.
	 *
	 * This method should preferably be used intead of update in a OpenMP context.
	 *
	 * Iff called by the master thread, it will also update the display if needed
	 *
	 * @param amount the number of newly performed tasks to report
	 *
	 * @return false iff the computation is aborted
	 */
	bool increment(unsigned long amount=1) { return monitor().increment(amount); }

	/**
	 * set the current progress indicator
	 *
	 * Iff called by the master thread, it will also update the display if needed
	 *
	 * @param current the total number of performed tasks so far (by all threads)
	 *
	 * @return false iff the computation is aborted
	 */
	bool update(unsigned long current) { return monitor().update(current); }

	/**
	 * return if the computation has been aborted.
	 * N.B: do not perform any check by itselfd
	 */
	bool is_aborted() const { return monitor().is_aborted(); }

	/**
	 * check that the no interruption has been requested and return the current status
	 *
	 * Iff called by the master thread, it will check for R-user level interruption.
	 *
	 * @return true iff the computation is aborted
	 */
	static bool check_abort() { return monitor().check_abort(); }

private:
	static InterruptableProgressMonitor*& monitor_singleton() {
		static InterruptableProgressMonitor* p = 0;
		return p;
	}

	// trick to provide a default static member in a header file
	static SimpleProgressBar& default_progress_bar() {
	  static SimpleProgressBar pb;
	  pb.reset();
	  return pb;
	}

public: // ==== OTHER PUBLIC INTERFACE =====
	static InterruptableProgressMonitor& monitor() { return *monitor_singleton(); }

};

#endif
