
#ifndef AMORTI_HPP
#define AMORTI_HPP
#include <iostream>
#include <array>
#include <vector>
#include "common.h"
#include "protein.hpp"

using namespace std;
#define K 0.05
#define MIN_INERTIA 0.7
#define MAX_INERTIA 100.0
#define DAMP_CUTOFF 0.05

struct Damp {
	// we use 4 coordinates proteins
	typedef Protein<4> Protein;
	// we need only one parameter (beta)
	static const unsigned int nbParams = 1;
	// and we produce 3 dimensional signatures (enhnance, inhibit, damp)
	static const unsigned int nbSignatureParams = 3;

	static const array<pair<double, double>, nbParams> paramsLimits() {
		return {{pair<double, double>{0.0, 20.0}}};
	}

	// returns the influence (enhance, inhibit & damp coef) of protein p0 onto p1
	static array<double, nbSignatureParams>
		getInfluence(const Protein &p0, const Protein &p1,
				const array<double, nbParams> &params) {
			return {{exp(-params[0] * abs(p0.coords.at(1) - p1.coords.at(0))),
				exp(-params[0] * abs(p0.coords.at(2) - p1.coords.at(0))),
				exp(-params[0] * abs(p0.coords.at(3) - p1.coords.at(0)))}};
		}

	// computes the new concentration for a protein p
	// and the influences (signature * concentration) of nbP proteins
	static double
		computeNewConcentration(const Protein &p,
				const array<double, nbSignatureParams> &influences,
				const array<double, nbParams> &params, const int nbP) {

			double d = max(influences[2] / (double)nbP - DAMP_CUTOFF, 0.0);
			double v = p.c - p.prevc;
			double inertia = mix(MIN_INERTIA, MAX_INERTIA, p.coords.at(3));
			double damp = 2.0 * d * 2.0 * sqrt(influences[2] * K);
			double forces = ((influences[0] - influences[1]) / (double)nbP) +
				(K * (1.0 - 2.0 * p.c)) - (2.0 * damp * v);
			return max(0.0, min(p.c + v + (forces / inertia), 1.0));
		}
};
#endif
