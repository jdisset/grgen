#include <algorithm>
#include "../common.h"
#include "../json/json.hpp"
#include "../mgclassic.hpp"
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
		if (p.output) o = true;
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
	for (auto &e : m) {
		if (e.second == elem) return true;
	}

	return false;
}

template <typename G> void checkProteinsLimits(G &g) {
	for (auto &p : g.actualProteins) {
		for (auto &c : p.coords) {
			REQUIRE(c <= G::Protein::getMaxCoord());
			REQUIRE(c >= G::Protein::getMinCoord());
		}
		REQUIRE(p.c >= 0);
		REQUIRE(p.prevc >= 0);
	}
}

template <typename G> void checkMGRNIntegrity(G &g, G *topLevel) {
	REQUIRE(g.master == topLevel);
	for (auto &sg : g.subNets) REQUIRE(sg.parent == &g);
	atLeastOneInputOneOutput(g);
	checkProteinsPtrIntegrity(g);
	checkProteinsLimits(g);
	for (auto &sg : g.subNets) checkMGRNIntegrity(sg, topLevel);
}

template <typename G> void checkProteinsPtrIntegrity(G &g) {
	size_t s = 0;
	size_t i = 0;
	size_t o = 0;
	bool topLevel = g.master == &g;
	REQUIRE(topLevel == g.isMaster());
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
	if (topLevel) {
		REQUIRE(g.inputProteins.size() == i);
		REQUIRE(g.outputProteins.size() == o);
	}
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
	REQUIRE(mg1.params.size() == mg1.params.size());
	auto params0 = mg0.params;
	auto params1 = mg1.params;
	for (size_t i = 0; i < params0.size(); ++i) {
		REQUIRE(params0[i] == params1[i]);
	}
	REQUIRE(mg1.actualProteins.size() == mg0.actualProteins.size());
	REQUIRE(mg1.subNets.size() == mg0.subNets.size());
	checkMGRNIntegrity(mg1, &mg1);
	REQUIRE(mg1 == mg0);
	mg0.addSubNet(mg1);
	REQUIRE(mg1 != mg0);
	REQUIRE(mg0.subNets.size() == 1);
	REQUIRE(mg0.subNets[0] == mg1);
	REQUIRE(mg0.subNets[0].parent == &mg0);
	REQUIRE(mg0.getListOfAllGRNs().size() == 2);
	REQUIRE(mg0.getListOfAllProteins().size() == 10);

	auto p = Protein();
	mg0.subNets[0].addProtein(p);
	REQUIRE(mg0.subNets[0].actualProteins.size() == mg0.actualProteins.size() + 1u);
	REQUIRE(mg0.subNets[0].subNets.size() == 0);
	checkMGRNIntegrity(mg0, &mg0);
	mg0.addRandomSubNet(1);  // add a random subnet with one protein
	REQUIRE(mg0.subNets.size() == 2);
	checkMGRNIntegrity(mg0, &mg0);
	mg0.subNets[1].addRandomSubNet(15);
	REQUIRE(mg0.subNets[1].subNets.size() == 1);
	REQUIRE(mg0.subNets[1].subNets[0].actualProteins.size() == 15);
	REQUIRE(mg0.subNets.size() == 2);
	checkMGRNIntegrity(mg0, &mg0);
	mg0.addRandomProtein(30);  // adding 30 completely random proteins
	mg0.subNets[1].addRandomSubNet(5);
	mg0.subNets[1].addRandomSubNet();
	mg0.subNets[0].addRandomSubNet(3);
	mg0.subNets[0].subNets[0].addRandomSubNet(12);
	REQUIRE(mg0.subNets[1].subNets.size() == 3);
	REQUIRE(mg0.subNets[0].subNets.size() == 1);
	REQUIRE(mg0.subNets[0].subNets[0].subNets.size() == 1);
	REQUIRE(mg0.subNets.size() == 2);
	checkMGRNIntegrity(mg0, &mg0);
	checkMGRNIntegrity(mg1, &mg1);  // mg1 should still be the same, unaffected

	return mg0;
}

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
	auto mgrn = constructRandomMGRN<T>();

	SECTION("deletion") {
		std::cerr << "deletion" << std::endl;
		size_t nbP = mgrn.getNbOwnProteins();
		REQUIRE(nbP == mgrn.actualProteins.size());
		size_t nbIn = mgrn.getNbOwnProteins(ProteinType::input);
		size_t nbOut = mgrn.getNbOwnProteins(ProteinType::output);
		size_t N =
		    static_cast<size_t>(floor(mgrn.getNbOwnProteins(ProteinType::regul) * 0.9f));
		for (size_t i = 0; i < N; ++i) mgrn.deleteRandomRegul();
		checkMGRNIntegrity(mgrn, &mgrn);
		REQUIRE(mgrn.getNbOwnProteins() == nbP - N);
		REQUIRE(mgrn.getNbOwnProteins(ProteinType::input) == nbIn);
		REQUIRE(mgrn.getNbOwnProteins(ProteinType::output) == nbOut);
		REQUIRE(mgrn.getNbOwnProteins(ProteinType::regul) == nbP - nbIn - nbOut - N);
		for (int i = 0; i < 500; ++i) mgrn.subNets[0].deleteRandomProtein();
		checkMGRNIntegrity(mgrn, &mgrn);
	}

	SECTION("serialization") {
		std::cerr << "serialization" << std::endl;
		auto jsonStr = mgrn.serialize();
		MGRN<T> mcopy(jsonStr);
		REQUIRE(mcopy == mgrn);
		checkMGRNIntegrity(mcopy, &mcopy);
		auto copyStr = mcopy.serialize();
		REQUIRE(copyStr == jsonStr);
		mcopy.setProteinConcentration("i0", ProteinType::input, 0.1234654567891253443456789);
		REQUIRE(mcopy != mgrn);
		MGRN<T> mcopy2(mcopy.serialize());
		REQUIRE(mcopy2 == mcopy);
	}

	SECTION("deterministic") {
		std::cerr << "determinism" << std::endl;
		deterministicGRN<T>();
	}

	SECTION("mutation & crossover") {
		std::cerr << "mutation" << std::endl;
		auto otherCopy = mgrn;
		REQUIRE(otherCopy == mgrn);
		REQUIRE(MGRN<T>::relativeDistance(mgrn, otherCopy) == 0);
		double avgD = 0;
		double nbMuts = 200;
		for (double i = 0; i < nbMuts; ++i) {
			auto tmp = mgrn;
			mgrn.mutate();
			double dist = MGRN<T>::relativeDistance(tmp, mgrn);
			avgD += dist;
			REQUIRE(dist < 1.0);
			checkMGRNIntegrity(mgrn, &mgrn);
		}
		avgD /= nbMuts;
		REQUIRE(avgD > 0.0);
		REQUIRE(avgD < 0.8);
		REQUIRE(otherCopy != mgrn);
		checkMGRNIntegrity(mgrn, &mgrn);

		std::cerr << "crossover" << std::endl;
		int nbCrossovers = 200;
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
		REQUIRE(std::abs(avgDistDiff) < 0.15);
	}
}

TEST_CASE("MGRN declaration, init & serialization", "[mgrn]") { testMGRN<MGClassic>(); }
