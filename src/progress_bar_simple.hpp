/*
 * simple_progress_bar.hpp
 *
 * A class that display a progress bar
 *
 * Author: karl.forner@gmail.com
 *
 */
#ifndef _RcppProgress_SIMPLE_PROGRESS_BAR_HPP
#define _RcppProgress_SIMPLE_PROGRESS_BAR_HPP

#include "progress_bar.hpp"

#include <R_ext/Print.h>

// for unices only
#if !defined(WIN32) && !defined(__WIN32) && !defined(__WIN32__)
#include <Rinterface.h>
#endif

class SimpleProgressBar: public ProgressBar{
  public: // ====== LIFECYCLE =====

    /**
    * Main constructor
    */
    SimpleProgressBar()  { reset(); }

    ~SimpleProgressBar() {}

  public: // ===== main methods =====

    void display() {
//      REprintf("0%%   10   20   30   40   50   60   70   80   90   100%%\n");
//      REprintf("[----|----|----|----|----|----|----|----|----|----|\n");
//      REprintf("[----|----|----|----|----|----|----|----|----|----|");
		REprintf("=====|====|====|====|====|====|====|====|====|====");
		REprintf("\r");
		flush_console();
    }

    // will finalize display if needed
    void update(float progress) {
      _update_ticks_display(progress);
      if (_ticks_displayed >= _max_ticks)
        _finalize_display();
    }

    void end_display() {
      update(1);
      reset();
    }

    void reset() {
      _max_ticks = 50;
      _ticks_displayed = 0;
      _finalized = false;
    }


  protected: // ==== other instance methods =====

    // update the ticks display corresponding to progress
    void _update_ticks_display(float progress) {
      int nb_ticks = _compute_nb_ticks(progress);
      int delta = nb_ticks - _ticks_displayed;
      if (delta > 0) {
        _display_ticks(delta);
        _ticks_displayed = nb_ticks;
      }

    }

    void _finalize_display() {
      if (_finalized) return;

  //    REprintf("|\n");
      REprintf("\n");
      flush_console();
      _finalized = true;
    }

    int _compute_nb_ticks(float progress) {
      return int(progress * _max_ticks);
    }

    void _display_ticks(int nb) {
      for (int i = 0; i < nb; ++i) {
      //  REprintf("*");
	    REprintf("-");
        R_FlushConsole();
      }
    }

    // N.B: does nothing on windows
    void flush_console() {
#if !defined(WIN32) && !defined(__WIN32) && !defined(__WIN32__)
       R_FlushConsole();
#endif
    }

  private:
    int _max_ticks;   		// the total number of ticks to print
    int _ticks_displayed; 	// the nb of ticks already displayed
    bool _finalized;

};

#endif
