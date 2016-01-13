#ifndef CLASSIC_HPP
#define CLASSIC_HPP
#include <iostream>
#include <array>
#include <vector>
#include "common.h"
#include "protein.hpp"

using namespace std;

struct Classic {
	// we use 3 coordinates proteins (id, enh, inh)
	typedef Protein<3> Protein;
	// we need 2 parameters (beta, alpha)
	static const unsigned int nbParams = 2;
	// and we produce 2 dimensional signatures (enhnance, inhibit)
	static const unsigned int nbSignatureParams = 2;

	static const array<pair<double, double>, nbParams> paramsLimits() {
		return {{{pair<double, double>{0.0, 20.0}}, {pair<double, double>{0.0, 1.0}}}};
	}

	// returns the influence (enhance, inhibit, ...) of protein p0 onto p1
	static array<double, nbSignatureParams> getInfluence(
	    const Protein &p0, const Protein &p1, const array<double, nbParams> &params) {
		return {{
		    exp(-params[0] * abs(p0.coords.at(1) - p1.coords.at(0))),
		    exp(-params[0] * abs(p0.coords.at(2) - p1.coords.at(0))),
		}};
	}

	// computes the new concentration for a protein p
	// and the influences (signature * concentration) of nbP proteins
	static double computeNewConcentration(
	    const Protein &p, const array<double, nbSignatureParams> &influences,
	    const array<double, nbParams> &params, const int nbP) {
		double forces = ((influences[0] - influences[1]) / (double)nbP) * params[1];
		return max(0.0, min(p.c + forces, 1.0));
	}
};
#endif
