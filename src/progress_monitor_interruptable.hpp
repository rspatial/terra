/*
 * interruptable_progress_monitor.hpp
 *
 * A class that monitors the progress of computations:
 *   - can display a progress bar
 *   - can handle user interrupt (R user level) or programmatic abort
 *   - can be used in OpenMP loops
 *
 * Author: karl.forner@gmail.com
 *
 */
#ifndef _RcppProgress_INTERRUPTABLE_PROGRESS_MONITOR_HPP
#define _RcppProgress_INTERRUPTABLE_PROGRESS_MONITOR_HPP

#include "interrupts.hpp"
#include "progress_bar.hpp"
//#include "fetch_raw_gwas_bar.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

class InterruptableProgressMonitor {
public: // ====== LIFECYCLE =====

	/**
	 * Main constructor
	 *
	 * @param max the expected number of tasks to perform
	 * @param display_progress whether to display a progress bar in the console
	 * @param pb    the ProgressBar instance to use
	 */
	InterruptableProgressMonitor(
	  unsigned long max,
	  bool display_progress,
	  ProgressBar& pb
  ) : _progress_bar(pb)
{
		reset(max, display_progress);
		if (is_display_on()) {
		  _progress_bar.display();
		}
	}

	~InterruptableProgressMonitor() {
		if (is_display_on() && !is_aborted()) _progress_bar.end_display();
	}

public: // ===== ACCESSORS/SETTERS =====
	void set_display_status(bool on) { _display_progress = on; 	}
	bool is_display_on() const { return _display_progress; }
	unsigned long get_max() const { return _max; }
	bool is_aborted() const { return _abort; }


public: // ===== PBLIC MAIN INTERFACE =====
	/**
	 * increment the current progress.
	 *
	 * Iff called by the master thread, it will also update the display if needed
	 *
	 * @param amount the number of newly performed tasks to report
	 *
	 * @return false iff the computation is aborted
	 */
	bool increment(unsigned long amount=1) {
		if ( is_aborted() )
			return false;
		return is_master() ? update_master(_current + amount) : atomic_increment(amount);
	}

	/**
	 * set the current progress indicator
	 *
	 * Iff called by the master thread, it will also update the display if needed
	 *
	 * @param current the total number of performed tasks so far (by all threads)
	 *
	 * @return false iff the computation is aborted
	 */
	bool update(unsigned long current) {
		if ( is_aborted() )
			return false;
		return is_master() ? update_master(current) : atomic_update(current);
	}

	/**
	 * check that the no interruption has been requested and return the current status
	 *
	 * Iff called by the master thread, it will check for R-user level interruption.
	 *
	 * @return true iff the computation is aborted
	 */
	bool check_abort() {
		if ( is_aborted() )
			return true;

		if ( is_master() )  {
			check_user_interrupt_master();
		}
		return is_aborted();
	}

	/**
	 * request computation abortion
	 */
	void abort() {
#ifdef _OPENMP
#pragma omp critical
#endif
		_abort = true;

	}

	/**
	 * return true iff the thread is the master.
	 * In case of non-OpenMP loop, always return true
	 */
	bool is_master() const {
#ifdef _OPENMP
		return omp_get_thread_num() == 0;
#else
		return true;
#endif
	}

public: // ===== methods for MASTER thread =====

	/**
	 * set the current progress indicator and update the progress bar display if needed.
	 *
	 *
	 * @param current the total number of performed tasks
	 *
	 * @return false iff the computation is aborted
	 */
	bool update_master(unsigned long current) {
		_current = current;
		if (is_display_on()) _progress_bar.update(progress(current));
		return ! is_aborted();
	}

	void check_user_interrupt_master() {
		if ( !is_aborted() && checkInterrupt() ) {
			abort();
		}
	}

public: // ===== methods for non-MASTER threads =====

	bool atomic_increment(unsigned long amount=1) {
#ifdef _OPENMP
#pragma omp atomic
#endif
		_current+=amount;
		return ! is_aborted();
	}

	bool atomic_update(unsigned long current) {
#ifdef _OPENMP
#pragma omp critical
#endif
		_current=current;
		return ! is_aborted();
	}

protected: // ==== other instance methods =====
	// convert current value to [0-1] progress
	double progress(unsigned long current) {
		return double(current) / double(_max);
	}

	/**
	 * reset the monitor.
	 *
	 * Currently not really useful
	 *
	 * @param max the expected number of tasks to perform
	 * @param display_progress whether to display a progress bar in the console
	 *
	 */
	void reset(unsigned long max = 1, bool display_progress = true) {
		_max = max;
		if ( _max <= 0 )
			_max = 1;
		_current = 0;
		_display_progress = display_progress;
		_abort = false;
	}


private: // ===== INSTANCE VARIABLES ====
  ProgressBar& _progress_bar;
	unsigned long _max; 			// the nb of tasks to perform
	unsigned long _current; 		// the current nb of tasks performed

	bool _abort;					// whether the process should abort
	bool _display_progress;			// whether to display the progress bar

};

#endif
