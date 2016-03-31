#include <algorithm>
#include "../common.h"
#include "../json/json.hpp"
#include "../mgrn.hpp"
#include "catch.hpp"

template <typename T> bool inVector(std::vector<T> vec, T elem) {
	return std::find(vec.begin(), vec.end(), elem) != vec.end();
}

template <typename G> void atLeastOneInputOneOutput(G &g) {
	bool i = false;
	bool o = false;
	for (auto &p : g.actualProteins) {
		if (p.input) i = true;
		if (p.output) i = true;
	}
	REQUIRE(i);
	REQUIRE(o);
}

template <typename T> bool ptrInVector(std::vector<T> vec, T *adr) {
	for (auto &e : vec)
		if (&e == adr) return true;
	return false;
}

template <typename K, typename T> bool inMap(std::map<K, T> m, T elem) {
	for (auto &e : m)
		if (&e.second == elem) return true;
	return false;
}

template <typename P> void checkProteinLimits(P &p) {
	for (auto &c : p.coords) {
		REQUIRE(c < P::max_coord);
		REQUIRE(c > P::min_coord);
	}
	REQUIRE(p.c > 0);
}

template <typename G> void checkProteinsLimits(G &g) {
	for (auto &p : g.actualProteins) checkProteinLimits(p);
}

template <typename G> void checkMGRNIntegrity(G &g, G *topLevel) {
	REQUIRE(g.master == topLevel);
	atLeastOneInputOneOutput(g);
	checkProteinsPtrIntegrity(g);
	checkProteinsLimits(g);
	for (auto &sg : g.subNets) {
		REQUIRE(sg.parent == &sg);
		checkMGRNIntegrity(sg, topLevel);
	}
}

template <typename G> void checkProteinsPtrIntegrity(G &g) {
	size_t s = 0;
	size_t i = 0;
	size_t o = 0;
	bool topLevel = g.master == &g;
	for (auto &p : g.actualProteins) {
		REQUIRE(inVector(g.allProteinsPtr, {&p, true}));
		if (topLevel) {
			if (p.input) {
				REQUIRE(inMap(g.inputProteins, s));
				++i;
			}
			if (p.output) {
				REQUIRE(inMap(g.outputProteins, s));
				++o;
			}
		}
		++s;
	}
	for (auto &sg : g.subNets) {
		for (auto &p : sg.actualProteins) {
			if (p.input || p.output) {
				REQUIRE(inVector(g.allProteinsPtr, {&p, false}));
				++s;
			}
		}
	}
	REQUIRE(g.allProteinsPtr.size() == s);
	REQUIRE(g.inputProteins.size() == i);
	REQUIRE(g.outputProteins.size() == o);
}

template <typename T> MGRN<T> constructRandomMGRN() {
	MGRN<T> mg0;
	using Protein = typename T::Protein_t;
	mg0.addRandomProtein(ProteinType::input, "i0");
	mg0.addRandomProtein(ProteinType::input, "i1");
	mg0.addRandomProtein(ProteinType::input, "i2");
	mg0.addRandomProtein(ProteinType::output, "o0");
	mg0.addRandomProtein(ProteinType::output, "o1");
	REQUIRE(mg0.actualProteins.size() == 5);
	REQUIRE(mg0.getNbProteins() == 5);
	REQUIRE(mg0.allProteinsPtr.size() == 5);
	REQUIRE(mg0.subNets.size() == 0);
	checkMGRNIntegrity(mg0, &mg0);

	auto mg1 = mg0;
	REQUIRE(mg1.getParams().size() == mg1.getParams().size());
	auto params0 = mg0.getParams();
	auto params1 = mg1.getParams();
	for (size_t i = 0; i < params0.size(); ++i) {
		REQUIRE(params0[i] == params1[i]);
	}
	REQUIRE(mg1.actualProteins.size() == mg0.actualProteins.size());
	REQUIRE(mg1.subNets.size() == mg0.subNets.size());
	CheckMGRNIntegrity(mg1, &mg1);
	REQUIRE(mg1 == mg0);
	mg0.addSubNet(mg1);
	REQUIRE(mg1 != mg0);
	REQUIRE(mg0.subNets.size() == 1);
	REQUIRE(mg0.subNets[0] == mg1);

	auto p = Protein();
	mg0.subNets[0].addProtein(p);
	REQUIRE(mg0.subNets[0].actualProteins.size() == mg0.actualProteins.size() + 1u);
	checkMGRNIntegrity(mg0, &mg0);

	mg0.addRandomSubNet(1);  // add a random subnet with one protein
	REQUIRE(mg0.subNets.size() == 2);
	checkMGRNIntegrity(mg0, &mg0);
	mg0.subNets[1].addRandomSubNet(15);
	checkMGRNIntegrity(mg0, &mg0);
	mg0.addRandomProtein(30);  // adding 30 completely random proteins
	mg0.subNets[1].addRandomSubNet(5);
	mg0.subNets[1].addRandomSubNet(0);
	mg0.subNets[0].addRandomSubNet(3);
	mg0.subNets[0].subNets[0].addRandomSubNet(12);
	REQUIRE(mg0.subNets[1].subNets.size() == 3);
	REQUIRE(mg0.subNets[0].subNets.size() == 2);
	REQUIRE(mg0.subNets.size() == 2);
	checkMGRNIntegrity(mg0, &mg0);
	checkMGRNIntegrity(mg1, &mg1);  // mg1 should still be the same, unaffected
	return mg0;
};

template <typename T> void deterministicGRN() {
	for (int n = 0; n < 100; ++n) {
		auto grn0 = constructRandomMGRN<T>();
		auto grn1 = grn0;
		grn0.reset();
		grn1.reset();
		REQUIRE(grn0 == grn1);
		grn0.step(10);
		grn1.step(10);
		REQUIRE(grn0 == grn1);
		auto grn2 = grn1;
		REQUIRE(grn0 == grn2);
		grn0.step(100);
		grn1.step(100);
		grn2.step(100);
		REQUIRE(grn0 == grn1);
		REQUIRE(grn0 == grn2);
	}
}

template <typename T> void testMGRN() {
	// construction
	auto mgrn = constructRandomMGRN<T>();

	// deletion
	size_t nbP = mgrn.getNbOwnProteins();
	REQUIRE(nbP == mgrn.actualProteins.size());
	size_t nbIn = mgrn.getNbOwnProteins(ProteinType::input);
	size_t nbOut = mgrn.getNbOwnProteins(ProteinType::output);
	size_t N = mgrn.getNbOwnProteins(ProteinType::regul) * 0.9;
	for (size_t i = 0; i < N; ++i) mgrn.deleteRandomRegul();
	checkMGRNIntegrity(mgrn, &mgrn);
	REQUIRE(mgrn.getNbOwnProteins() == nbP - N);
	REQUIRE(mgrn.getNbOwnProteins(ProteinType::input) == nbIn);
	REQUIRE(mgrn.getNbOwnProteins(ProteinType::output) == nbOut);
	REQUIRE(mgrn.getNbOwnProteins(ProteinType::regul) == nbP - nbIn - nbOut - N);
	for (int i = 0; i < 500; ++i) mgrn.subNets[0].deleteRandomProtein();
	checkMGRNIntegrity(mgrn, &mgrn);

	// serialization / deserialization
	auto jsonStr = mgrn.toJSON();
	MGRN<T> mcopy(jsonStr);
	checkMGRNIntegrity(mcopy);
	REQUIRE(mcopy == mgrn);
	auto copyStr = mcopy.toJSON();
	REQUIRE(copyStr == jsonStr);

	deterministicGRN<T>();

	// mutation
	auto otherCopy = mgrn;
	REQUIRE(otherCopy == mgrn);
	REQUIRE(MGRN<T>::relativeDistance(mgrn, otherCopy) == 0);
	for (int i = 0; i < 200; ++i) {
		auto tmp = mgrn;
		mgrn.mutate();
		double dist = MGRN<T>::relativeDistance(tmp, mgrn);
		REQUIRE(dist > 0.0);
		REQUIRE(dist < 1.0);
	}
	REQUIRE(otherCopy != mgrn);
	checkMGRNIntegrity(mgrn, &mgrn);

	// crossover
	int nbDifferentOffspring = 0;
	int nbCrossovers = 150;
	double avgDist0 = 0.0;
	double avgDist1 = 0.0;
	double avgDistDiff = 0.0;
	for (int i = 0; i < nbCrossovers; ++i) {
		auto offspring = MGRN<T>::crossover(otherCopy, mgrn);
		checkMGRNIntegrity(offspring, &offspring);
		double dist0 = MGRN<T>::relativeDistance(mgrn, offspring);
		double dist1 = MGRN<T>::relativeDistance(otherCopy, offspring);
		double distDiff = dist0 - dist1;
		avgDist0 += dist0;
		avgDist1 += dist1;
		avgDistDiff += distDiff;
	}
	avgDist0 /= (double)nbCrossovers;
	avgDist1 /= (double)nbCrossovers;
	avgDistDiff /= (double)nbCrossovers;
	REQUIRE(avgDist0 > 0);
	REQUIRE(avgDist1 > 0);
	REQUIRE(std::abs(avgDistDiff) < 0.05);
	REQUIRE(nbDifferentOffspring > nbCrossovers * 0.5);
}

TEST_CASE("MGRN declaration, init & serialization", "[mgrn]") {
	testMGRN<Classic>();
	testMGRN<RealCircle>();
}
