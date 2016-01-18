#ifndef AMORTI_HPP
#define AMORTI_HPP
#include <iostream>
#include <array>
#include <vector>
#include <unordered_map>
#include <utility>
#include "common.h"
#include "protein.hpp"

using namespace std;
#define K 0.05
#define MIN_INERTIA 0.7
#define MAX_INERTIA 100.0
#define DAMP_CUTOFF 0.05

struct Damp {
	// we use 4 coordinates proteins
	using Protein_t = Protein<4>;
	// we need only one parameter (beta)
	static const unsigned int nbParams = 1;
	// and we produce 3 dimensional signatures (enhnance, inhibit, damp)
	static const unsigned int nbSignatureParams = 3;
	// aliases for ProteinType
	static constexpr ProteinType pinput = ProteinType::input;
	static constexpr ProteinType pregul = ProteinType::regul;
	static constexpr ProteinType poutput = ProteinType::output;

	static const array<pair<double, double>, nbParams> paramsLimits() {
		return {{{0.0, 20.0}}};
	}

	// returns the influence (enhance, inhibit & damp coef) of protein p0 onto p1

	// computes the new concentration for a protein p
	// and the influences (signature * concentration) of nbP proteins

	template <typename GRN> void updateSignatures(GRN& grn) {
		grn.signatures.clear();
		eachP<Protein_t>({pinput, pregul, poutput}, grn.proteins, [&](Protein_t& p0) {
			std::unordered_map<Protein_t*, array<double, nbSignatureParams>> buffer;
			eachP<Protein_t>({pinput, pregul, poutput}, grn.proteins, [&](Protein_t& p1) {
				buffer[&p1] = {{exp(-grn.params[0] * abs(p0.coords.at(1) - p1.coords.at(0))),
				                exp(-grn.params[0] * abs(p0.coords.at(2) - p1.coords.at(0))),
				                exp(-grn.params[0] * abs(p0.coords.at(3) - p1.coords.at(0)))}};
			});
			grn.signatures[&p0] = buffer;
		});
	}

	template <typename GRN> void step(GRN& grn, unsigned int nbSteps) {
		for (auto s = 0u; s < nbSteps; ++s) {
			std::array<std::unordered_map<Protein_t*, double>, 3>
			    buffer;  // replacement proteins
			buffer[1].reserve(grn.getProteinSize(ProteinType::regul));
			buffer[2].reserve(grn.getProteinSize(ProteinType::output));
			for (auto t = 1u; t < 3; ++t) {
				for (auto& pr0 : grn.proteins[t]) {
					Protein_t* p0 = &pr0.second;
					array<double, nbSignatureParams> influence{};
					double n = 0;
					for (auto t1 = 0u; t1 < 2; ++t1) {
						for (auto& pr1 : grn.proteins[t1]) {
							Protein_t* p1 = &pr1.second;
							for (size_t i = 0; i < influence.size(); ++i) {
								influence[i] += p1->c * grn.signatures.at(p1).at(p0).at(i);
							}
							++n;
						}
					}
					double nbp = static_cast<double>(grn.getNbProteins());
					double d = std::max((influence[2] / nbp) - DAMP_CUTOFF, 0.0);
					double v = p0->c - p0->prevc;
					double inertia = mix(MIN_INERTIA, MAX_INERTIA, p0->coords.at(3));
					double damp = 2.0 * d * 2.0 * sqrt(influence[2] * K);
					double forces = ((influence[0] - influence[1]) / nbp) +
					                (K * (1.0 - 2.0 * p0->c)) - (2.0 * damp * v);
					buffer[t].emplace(p0,
					                  std::max(0.0, std::min(p0->c + v + (forces / inertia), 1.0)));
				}
			}
			for (auto t = 1u; t < 3; ++t) {
				for (auto& pr : grn.proteins[t]) {
					auto& p = pr.second;
					p.prevc = p.c;
					p.c = buffer[t].at(&p);
				}
			}
		}
	}
};
#endif
