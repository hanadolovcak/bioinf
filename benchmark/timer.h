// Author: Hana Dolovčak
#ifndef BAMBOO_TIMER_H
#define BAMBOO_TIMER_H

#include <chrono>

namespace bamboo {

// Simple wall-clock timer based on std::chrono::steady_clock.
class Timer {
public:
    // Starts the timer.
    Timer() { reset(); }

    // Restarts the timer from now.
    void reset() {
        start_ = std::chrono::steady_clock::now();
    }

    // Returns elapsed seconds since the timer was constructed or last reset.
    double seconds() const {
        auto now = std::chrono::steady_clock::now();
        std::chrono::duration<double> d = now - start_;
        return d.count();
    }

private:
    std::chrono::steady_clock::time_point start_;
};

}

#endif
