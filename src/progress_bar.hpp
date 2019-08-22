/*
 * progress_bar.hpp
 *
 * An abstract class for classes that display a progress bar
 *
 * Author: karl.forner@gmail.com
 *
 */
#ifndef _RcppProgress_PROGRESS_BAR_HPP
#define _RcppProgress_PROGRESS_BAR_HPP

class ProgressBar {
  public:

    // subclasses should not rely on the destructor to finalize the display
    virtual ~ProgressBar() = 0;

    // start the display. It will be updated by subsequent calls to update()
    virtual void display() = 0;

    // update if needed the display
    virtual void update(float progress) = 0;

    // finalize the display
    virtual void end_display() = 0;
};

inline ProgressBar::~ProgressBar() {}

#endif
