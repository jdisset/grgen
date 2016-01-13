#ifndef COMMON_H
#define COMMON_H
#include <iostream>
#include <random>
#include "json/json.hpp"

#define INIT_CONCENTRATION 0.5

static std::default_random_engine grnRand =
    std::default_random_engine(std::chrono::system_clock::now().time_since_epoch().count());
enum class ProteinType { input = 0, regul = 1, output = 2 };
template <typename T> T mix(const T &a, const T &b, const double v) {
	double r = v > 1.0 ? 1.0 : (v < 0.0 ? 0.0 : v);
	return (a * (1.0 - r)) + (r * b);
}


#endif
