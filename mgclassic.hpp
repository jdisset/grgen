#ifndef MGCLASSIC_HPP
#define MGCLASSIC_HPP
#include <array>
#include <iostream>
#include <unordered_map>
#include <utility>
#include <vector>
#include "common.h"
#include "protein.hpp"

using namespace std;

struct MGClassic {
 public:
	// we use 3 coordinates proteins (id, enh, inh)
	static constexpr int IDSIZE = 32;
	using Protein_t = HiProtein<3, int, 0, IDSIZE>;

	// we need 2 parameters (beta, alpha)
	static constexpr unsigned int nbParams = 2;
	// and we produce 2 dimensional signatures (enhnance, inhibit)
	static constexpr unsigned int nbSignatureParams = 2;

	static const array<pair<double, double>, nbParams> paramsLimits() {
		return {{{0.5, 2.0}, {0.5, 2.0}}};
	}

	// helpers for proteins coords
	static inline int& getId(Protein_t& p) { return p.coords[0]; }
	static inline int& getEnh(Protein_t& p) { return p.coords[1]; }
	static inline int& getInh(Protein_t& p) { return p.coords[2]; }

	using InfluenceVec = std::array<double, nbSignatureParams>;

	template <typename GRN>
	static InfluenceVec getInfluenceVec(Protein_t* p, Protein_t* p1, const GRN& grn) {
		// how p is influenced by p1
		return {{exp(-grn.params[0] * (abs(getId(*p) - getEnh(*p1)))),
		         exp(-grn.params[0] * (abs(getId(*p) - getInh(*p1))))}};
	}

	template <typename GRN> static void step(GRN& grn, unsigned int nbSteps) {
		for (unsigned int step = 0; step < nbSteps; ++step) {
			double sumConcentration = 0.0;
			for (auto& s : grn.signatures) {
				auto& p = s.first;
				double enh = 0.0, inh = 0.0;
				for (auto& infl : s.second) {
					enh += infl.first->c * infl.second[0];
					inh += infl.first->c * infl.second[1];
				}
				p->prevc = max(
				    0.0,
				    p->c + (grn.params[1] / static_cast<double>(s.second.size())) * (enh - inh));
				sumConcentration += p->prevc;
			}
			if (sumConcentration == 0) sumConcentration = 1;
			for (auto& s : grn.signatures) {
				auto tmp = s.first->prevc;
				s.first->prevc = s.first->c;
				s.first->c = tmp / sumConcentration;
			}
		}
	}
};
#endif
