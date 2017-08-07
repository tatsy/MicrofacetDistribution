#ifndef _COMMON_H_
#define _COMMON_H_

#include <cmath>

static const double Pi = 4.0 * std::atan(1.0);

static inline double sign(double x) {
    if (x == 0.0) {
        return 0.0;
    } else {
        return x / std::abs(x);
    }
}

inline double clamp(double x, double lower = 0.0, double higher = 1.0) {
    return std::max(lower, std::min(x, higher));
}

inline uint8_t toInt(double x) { return (uint8_t)(std::pow(clamp(x), 1.0 / 2.2) * 255.0 + 0.5); }


#endif  // _COMMON_H_
