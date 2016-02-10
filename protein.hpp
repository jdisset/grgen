#ifndef PROTEIN_HPP
#define PROTEIN_HPP

#include <assert.h>
#include <array>
#include <vector>
#include <unordered_set>
#include <string>
#include <random>
#include <type_traits>
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

	static bool areSame(const Protein &p0, const Protein &p1) {
		for (size_t i = 0; i < nbCoords; ++i)
			if (p0.coords[i] != p1.coords[i]) return false;
		return true;
	}

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

	Protein(const Protein &p) : coords(p.coords), c(p.c), prevc(p.prevc){};
	Protein() {
		// Constructs a protein with random coords
		for (auto &i : coords) i = getRandomCoord();
	}

	explicit Protein(const json &o) {
		// constructs a protein from a json object
		assert(o.count("coords"));
		json coordsArray = o.at("coords");
		assert(coordsArray.size() == nbCoords);
		size_t i = 0;
		for (auto &co : coordsArray) coords[i++] = co.get<CoordsType>();
	}

	void reset() {
		c = INIT_CONCENTRATION;
		prevc = c;
	}

	void mutate() {
		std::uniform_real_distribution<double> dReal(0.0, 1.0);  // just a dice roll
		std::uniform_int_distribution<int> dInt(0, nbCoords);
		size_t mutated = dInt(grnRand);
		coords[mutated] = getRandomCoord();
		// we want to mutate at least 1 coord, maybe more.
		for (size_t i = 0; i < nbCoords; ++i) {
			if (i != mutated && dReal(grnRand) < MULTIPLE_MUTATION_PROBA)
				coords[i] = getRandomCoord();
		}
	}

	json toJSON() const {
		json o;
		o["coords"] = coords;
		return o;
	}

	static constexpr double getMaxDistance() {
		return sqrt(std::pow((maxCoord - minCoord), 2) * nbCoords);
	}

	double getDistanceWith(const Protein &p) {
		double sum = 0;
		for (size_t i = 0; i < nbCoords; ++i) {
			sum += std::pow(
			    static_cast<double>(coords.at(i)) - static_cast<double>(p.coords.at(i)), 2);
		}
		return sqrt(sum) / getMaxDistance();
	}
};
#endif
