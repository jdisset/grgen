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
	static constexpr int IDSIZE = 32;
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
		grn.signatures.resize(grn.actualProteins.size());
		for (size_t i = 0; i < grn.actualProteins.size(); ++i) {
			grn.signatures[i].resize(grn.actualProteins.size());
			for (size_t j = 0; j < grn.actualProteins.size(); ++j) {
				auto& p0 = grn.actualProteins[i];
				auto& p1 = grn.actualProteins[j];
				grn.signatures[i][j] = {
				    {static_cast<double>(IDSIZE - abs(getEnh(p0) - getId(p1))),
				     static_cast<double>(IDSIZE - abs(getInh(p0) - getId(p1)))}};
				if (grn.signatures[i][j][0] > maxEnhance) maxEnhance = grn.signatures[i][j][0];
				if (grn.signatures[i][j][1] > maxEnhance) maxEnhance = grn.signatures[i][j][1];
			}
		}
		for (size_t i = 0; i < grn.actualProteins.size(); ++i) {
			for (size_t j = 0; j < grn.actualProteins.size(); ++j) {
				grn.signatures[i][j] = {
				    {exp(grn.params[0] * grn.signatures[i][j][0] - maxEnhance),
				     exp(grn.params[1] * grn.signatures[i][j][1] - maxInhibit)}};
			}
		}
	}

	template <typename GRN> void step(GRN& grn, unsigned int nbSteps) {
		// std::cerr << " grn.signatures" << std::endl;
		// for (size_t i = 0; i < grn.actualProteins.size(); ++i) {
		// for (size_t j = 0; j < grn.actualProteins.size(); ++j) {
		// std::cerr << "[" << i << "][" << j << "] = {" << grn.signatures[i][j][0] << ", "
		//<< grn.signatures[i][j][1] << "} ";
		//}
		// std::cerr << std::endl;
		//}
		for (auto s = 0u; s < nbSteps; ++s) {
			std::vector<double> nextProteins;  // only reguls & outputs concentrations
			nextProteins.reserve(grn.getNbProteins() - grn.getProteinSize(ProteinType::input));
			for (size_t j = grn.getFirstRegulIndex(); j < grn.getNbProteins(); ++j) {
				array<double, nbSignatureParams> influence{};
				for (size_t k = 0; k < grn.getFirstOutputIndex(); ++k) {
					for (size_t i = 0; i < influence.size(); ++i) {
						influence[i] += grn.actualProteins[k].c * grn.signatures[k][j][i];
					}
				}
				
				nextProteins.push_back(max(0.0, grn.actualProteins[j].c +
				                                    (grn.params[1] / grn.getNbProteins()) *
				                                        (influence[0] - influence[1])));
				std::cerr << "Influence = " << influence[0] << ", " << influence[1]
				          << ", pushed back " << nextProteins[nextProteins.size() - 1]
				          << std::endl;
			}
			// Normalizing regul & output proteins concentrations
			double sumConcentration = 0.0;
			for (auto i : nextProteins) {
				sumConcentration += i;
			}
			std::cerr << "sumc = " << sumConcentration << std::endl;
			if (sumConcentration > 0) {
				for (auto& i : nextProteins) {
					i /= sumConcentration;
				}
			}
			for (size_t i = grn.getFirstRegulIndex(); i < grn.getNbProteins(); ++i) {
				grn.actualProteins[i].c = nextProteins[i];
			}
			for (auto& i : nextProteins) {
				std::cerr << "np =" << i << std::endl;
			}
		}
		for (auto& t : grn.proteinsRefs) {
			for (auto& p : t) {
				std::cerr << " p " << p.first << " [" << p.second
				          << "] = " << grn.actualProteins[p.second].c << std::endl;
			}
		}
	}
};
#endif
