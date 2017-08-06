#ifndef _RANDOM_H_
#define _RANDOM_H_

#include <random>

//! Random number generator with C++11 Mersenne twister
class Random {
public:
    Random(uint32_t seed = 0u)
        : mt(seed) {
    }

    double next() {
        return dist(mt);
    }

private:
    std::mt19937 mt;
    std::uniform_real_distribution<double> dist{0.0, 1.0};
};

static thread_local Random random;

#endif  // _RANDOM_H_
