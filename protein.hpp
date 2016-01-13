#ifndef PROTEIN_HPP
#define PROTEIN_HPP

#include <assert.h>
#include <array>
#include <vector>
#include <set>
#include "common.h"

using namespace std;

#define MULTIPLE_MUTATION_PROBA 0.1
template <unsigned int nbCoords> struct Protein {
	using json = nlohmann::json;
	array<double, nbCoords> coords{};
	double c = INIT_CONCENTRATION;
	double prevc = INIT_CONCENTRATION;

	Protein(const Protein &p) : coords(p.coords), c(p.c), prevc(p.prevc){};
	Protein() {
		uniform_real_distribution<double> distribution(0.0, 1.0);
		for (auto &i : coords) {
			i = distribution(grnRand);
		}
	}
	void reset() {
		c = INIT_CONCENTRATION;
		prevc = c;
	}
	explicit Protein(const json &o) {
		assert(o.count("coords"));
		json coordsArray = o.at("coords");
		assert(coordsArray.size() == nbCoords);
		size_t i = 0;
		for (auto &co : coordsArray) {
			sscanf(co.get<string>().c_str(), "%lf", &coords[i++]);
		}
	}
	void mutate() {
		// we want to mutate at least 1 coord, maybe more.
		uniform_real_distribution<double> dReal(0.0, 1.0);
		uniform_int_distribution<int> dInt(0, nbCoords);
		set<size_t> toMutate;
		toMutate.insert(dInt(grnRand));
		for (size_t i = 0; i < nbCoords; ++i) {
			if (dReal(grnRand) < MULTIPLE_MUTATION_PROBA) toMutate.insert(i);
		}
		for (auto &i : toMutate) {
			coords[i] = dReal(grnRand);
		}
	}

	json toJSON() const {
		json o;
		json coordsArray;
		for (auto &co : coords) {
			char buf[50];
			snprintf(buf, sizeof(buf), "%a", co);
			coordsArray.push_back(buf);
		}
		o["coords"] = coordsArray;
		return o;
	}

	double getDistanceWith(const Protein<nbCoords> &p) {
		double sum = 0;
		for (size_t i = 0; i < nbCoords; ++i) {
			sum += pow(coords[i] - p.coords.at(i), 2);
		}
		return sqrt(sum);
	}
};
#endif
