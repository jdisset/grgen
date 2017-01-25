#ifndef PROTEIN_HPP
#define PROTEIN_HPP

#include <assert.h>
#include <array>
#include <random>
#include <string>
#include <type_traits>
#include <unordered_set>
#include <utility>
#include <vector>
#include "common.h"

#define MULTIPLE_MUTATION_PROBA 0.1
template <unsigned int nbCoords, typename CoordsType = double, int minCoord = 0,
          int maxCoord = 1>
struct Protein {
	using json = nlohmann::json;

	std::array<CoordsType, nbCoords>
	    coords{};                       // proteins coords (id, enh, inh for example)
	double c = INIT_CONCENTRATION;      // current concentration
	double prevc = INIT_CONCENTRATION;  // previous concentration

	bool operator==(const Protein &b) const {
		for (size_t i = 0; i < nbCoords; ++i)
			if (coords[i] != b.coords[i]) return false;
		return c == b.c;
	}
	bool operator!=(const Protein &b) const { return !(*this == b); }

	// switching between integral or real random distribution
	template <typename T = CoordsType>
	typename std::enable_if<!std::is_integral<T>::value, T>::type getRandomCoord() {
		std::uniform_real_distribution<double> distribution(static_cast<double>(minCoord),
		                                                    static_cast<double>(maxCoord));
		return static_cast<CoordsType>(distribution(grnRand));
	}
	template <typename T = CoordsType>
	typename std::enable_if<std::is_integral<T>::value, T>::type getRandomCoord() {
		std::uniform_int_distribution<int> distribution(minCoord, maxCoord);
		return static_cast<CoordsType>(distribution(grnRand));
	}

	Protein(const decltype(coords) &co, double conc) : coords(co), c(conc), prevc(conc) {
		for (auto &coo : coords) {
			coo = std::max(static_cast<CoordsType>(minCoord),
			               std::min(static_cast<CoordsType>(maxCoord), coo));
		}
	}
	Protein(const Protein &p) : coords(p.coords), c(p.c), prevc(p.prevc){};
	Protein() {
		// Constructs a protein with random coords
		for (auto &i : coords) i = getRandomCoord();
	}

	explicit Protein(const json &o) {
		// constructs a protein from a json object
		assert(o.count("coords"));
		assert(o.count("c"));
		c = o.at("c");
		if (o.count("pc")) prevc = o.at("pc");
		auto vcoords = o.at("coords").get<std::vector<CoordsType>>();
		assert(vcoords.size() == coords.size());
		for (size_t i = 0; i < vcoords.size(); ++i) coords[i] = vcoords[i];
	}

	void reset() {
		c = INIT_CONCENTRATION;
		prevc = c;
	}

	void mutate() {
		std::uniform_int_distribution<int> dInt(0, nbCoords);
		int mutated = dInt(grnRand);
		coords[mutated] = getRandomCoord();
	}

	json toJSON() const {
		json o;
		o["coords"] = coords;
		o["c"] = c;
		o["pc"] = prevc;
		return o;
	}

	static constexpr double getMaxDistance() {
		return sqrt(std::pow((maxCoord - minCoord), 2) * nbCoords);
	}

	// double getDistanceWith(const Protein &p) const {
	// double sum = 0;
	// for (size_t i = 0; i < nbCoords; ++i) {
	// sum += std::pow(
	// static_cast<double>(coords.at(i)) - static_cast<double>(p.coords.at(i)), 2);
	//}
	// return sqrt(sum) / getMaxDistance();
	//}
	double getDistanceWith(const Protein &p) const {
		return (0.75 * abs(p.coords.at(0) - coords.at(0)) +
		        0.125 * (abs(p.coords.at(1) - coords.at(1)) +
		                 abs(p.coords.at(2) - coords.at(2)))) /
		       static_cast<double>(maxCoord - minCoord);
	}
};
#endif
