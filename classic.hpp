#ifndef CLASSIC_HPP
#define CLASSIC_HPP
#include <iostream>
#include <array>
#include <vector>
#include <utility>
#include <unordered_map>
#include "common.h"
#include "protein.hpp"

using namespace std;

struct Classic {
	// we use 3 coordinates proteins (id, enh, inh)
	static constexpr int IDSIZE = 16;
	using Protein_t = Protein<3, int, 0, IDSIZE>;

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

	// aliases for ProteinType
	static constexpr ProteinType pinput = ProteinType::input;
	static constexpr ProteinType pregul = ProteinType::regul;
	static constexpr ProteinType poutput = ProteinType::output;

	int maxEnhance = 0, maxInhibit = 0;

	Classic() {}

	template <typename GRN> void updateSignatures(GRN& grn) {
		grn.signatures.clear();
		eachP<Protein_t>({pinput, pregul, poutput}, grn.proteins, [&](Protein_t& p0) {
			std::unordered_map<Protein_t*, array<double, nbSignatureParams>> buffer;
			eachP<Protein_t>({pinput, pregul, poutput}, grn.proteins, [&](Protein_t& p1) {
				buffer[&p1] = {{static_cast<double>(IDSIZE - abs(getEnh(p0) - getId(p1))),
				                static_cast<double>(IDSIZE - abs(getInh(p0) - getId(p1)))}};
				if (buffer[&p1][0] > maxEnhance) maxEnhance = buffer[&p1][0];
				if (buffer[&p1][1] > maxInhibit) maxInhibit = buffer[&p1][1];
			});
			grn.signatures[&p0] = buffer;
		});
		eachP<Protein_t>({pinput, pregul, poutput}, grn.proteins, [&](Protein_t& p0) {
			eachP<Protein_t>({pinput, pregul, poutput}, grn.proteins, [&](Protein_t& p1) {
				grn.signatures[&p0][&p1] = {
				    {exp(grn.params[0] * grn.signatures[&p0][&p1][0] - maxEnhance),
				     exp(grn.params[1] * grn.signatures[&p0][&p1][1] - maxInhibit)}};
			});
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
								influence[i] += p1->c * grn.signatures.at(p0).at(p1).at(i);
							}
							++n;
						}
					}
					buffer[t].emplace(p0, max(0.0, p0->c +
					                                   grn.params[1] / grn.getNbProteins() *
					                                       (influence[0] - influence[1])));
				}
			}
			// Normalizing regul & output proteins concentrations
			double sumConcentration = 0.0;
			for (auto t = 1u; t < 3; ++t) {
				for (auto& pr : grn.proteins[t]) {
					auto& p = pr.second;
					p.prevc = p.c;
					p.c = buffer[t].at(&p);
					sumConcentration += p.c;
				}
			}
			if (sumConcentration > 0) {
				for (auto t = 1u; t < 3; ++t) {
					for (auto& pr : grn.proteins[t]) {
						pr.second.c /= sumConcentration;
					}
				}
			}
		}
	}
};
#endif
