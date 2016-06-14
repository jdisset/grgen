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

template <typename T, typename U> bool inSignature(T *p, U &sig) {
	for (auto &s : sig)
		if (p == s.first) return true;
	return false;
}

template <typename T, typename U> auto getInfluencePool(T *p, U &g) {
	for (auto &s : g.signatures)
		if (p == s.first) return s.second;
	return typename decltype(g.signatures)::value_type::second_type();
}

template <typename T> void checkSignatures(T &g) {
	// signature should concern every updatable internal proteins
	for (auto &p : g.actualProteins) {
		if (!g.isMaster() || !p.input) {
			// updatable protein
			REQUIRE(inSignature(&p, g.signatures));
			auto pool = getInfluencePool(&p, g);
			if (p.input && !p.output && !g.isMaster()) {
				// its a pure input (= an output of the parent sub)
				// only influenced by parent
				size_t psize = 0;
				for (auto &p0 : g.parent->allProteinsPtr) {
					if ((p0.first->input && p0.second) || (!p0.first->input && !p0.first->output) ||
					    (p0.first->input && p0.first->output) || (p0.first->output && !p0.second)) {
						REQUIRE(inSignature(p0.first, pool));
						++psize;
					}
				}
				REQUIRE(psize == pool.size());
			}
			if ((p.output && !p.input) || (!p.input && !p.output)) {
				// only influenced by this level
				size_t psize = 0;
				for (auto &p0 : g.allProteinsPtr) {
					if ((p0.first->input && p0.second) || (!p0.first->input && !p0.first->output) ||
					    (p0.first->input && p0.first->output) || (p0.first->output && !p0.second)) {
						REQUIRE(inSignature(p0.first, pool));
						++psize;
					}
				}
				REQUIRE(psize == pool.size());
			}
			if ((p.input && !g.isMaster() && p.output)) {
				// influenced by both
				size_t psize = 0;
				for (auto &p0 : g.allProteinsPtr) {
					if ((p0.first->input && p0.second) || (!p0.first->input && !p0.first->output) ||
					    (p0.first->input && p0.first->output) || (p0.first->output && !p0.second)) {
						REQUIRE(inSignature(p0.first, pool));
						++psize;
					}
				}
				for (auto &p0 : g.parent->allProteinsPtr) {
					if ((p0.first->input && p0.second) || (!p0.first->input && !p0.first->output) ||
					    (p0.first->input && p0.first->output) || (p0.first->output && !p0.second)) {
						REQUIRE(inSignature(p0.first, pool));
						++psize;
					}
				}
				REQUIRE(psize == pool.size());
			}
		}
	}
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
	// checkSignatures(g);
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

	SECTION("copy & move") {
		// std::cerr << "copy & move" << std::endl;
		auto mgcopy = mgrn;
		REQUIRE(mgrn.getListOfAllProteins().size() == mgcopy.getListOfAllProteins().size());
		REQUIRE(mgcopy == mgrn);
		MGRN<T> othercopy(mgrn);
		REQUIRE(mgrn.getListOfAllProteins().size() ==
		        othercopy.getListOfAllProteins().size());
		REQUIRE(othercopy == mgrn);
		vector<MGRN<T>> vec0 = {{mgrn, mgcopy}};
		auto vec1 = vec0;
		auto vec2 = vec0;
		vec1.clear();
		vec1 = vec2;
		for (size_t i = 0; i < vec1.size(); ++i) REQUIRE(vec1[i] == vec0[i]);
	}

	SECTION("deletion") {
		// std::cerr << "deletion" << std::endl;
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
		// std::cerr << "serialization" << std::endl;
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
		// std::cerr << "determinism" << std::endl;
		// deterministicGRN<T>();
	}

	SECTION("mutation & crossover") {
		// std::cerr << "mutation" << std::endl;
		// auto otherCopy = mgrn;
		// REQUIRE(otherCopy == mgrn);
		// REQUIRE(MGRN<T>::relativeDistance(mgrn, otherCopy) == 0);
		// double avgD = 0;
		// double nbMuts = 200;
		// for (double i = 0; i < nbMuts; ++i) {
		// auto tmp = mgrn;
		// mgrn.mutate();
		// double dist = MGRN<T>::relativeDistance(tmp, mgrn);
		// avgD += dist;
		// REQUIRE(dist < 1.0);
		// checkMGRNIntegrity(mgrn, &mgrn);
		//}
		// avgD /= nbMuts;
		// REQUIRE(avgD > 0.0);
		// REQUIRE(avgD < 0.8);
		// REQUIRE(otherCopy != mgrn);
		// checkMGRNIntegrity(mgrn, &mgrn);

		//// std::cerr << "crossover" << std::endl;
		// int nbCrossovers = 200;
		// double avgDist0 = 0.0;
		// double avgDist1 = 0.0;
		// double avgDistDiff = 0.0;
		// for (int i = 0; i < nbCrossovers; ++i) {
		// auto offspring = MGRN<T>::crossover(otherCopy, mgrn);
		// checkMGRNIntegrity(offspring, &offspring);
		// double dist0 = MGRN<T>::relativeDistance(mgrn, offspring);
		// double dist1 = MGRN<T>::relativeDistance(otherCopy, offspring);
		// double distDiff = dist0 - dist1;
		// avgDist0 += dist0;
		// avgDist1 += dist1;
		// avgDistDiff += distDiff;
		//}
		// avgDist0 /= (double)nbCrossovers;
		// avgDist1 /= (double)nbCrossovers;
		// avgDistDiff /= (double)nbCrossovers;
		// REQUIRE(avgDist0 > 0);
		// REQUIRE(avgDist1 > 0);
		// REQUIRE(std::abs(avgDistDiff) < 0.15);
	}
}

TEST_CASE("MGRN random declaration, init & serialization", "[mgrn]") {
	for (int i = 0; i < 0; ++i) testMGRN<MGClassic>();
}

template <typename P, typename G>
std::pair<bool, typename G::InfluenceVec> getSignature(P *a, P *b, const G &g) {
	for (auto &s : g.signatures) {
		if (s.first == a) {
			for (auto &sa : s.second) {
				if (sa.first == b) {
					return {true, sa.second};
				}
			}
			break;
		}
	}
	return {false, typename G::InfluenceVec()};
}

template <typename P, typename G>
std::vector<std::pair<typename G::Protein *, typename G::InfluenceVec>> getInfluencesTo(
    P *a, const G &g) {
	for (auto &s : g.signatures)
		if (s.first == a) return s.second;
	return {};
}

template <class T> void printInfluences(const T &I) {
	std::cerr << "  --------- " << std::endl;
	for (const auto &i : I)
		std::cerr << " - with " << i.first->toJSON().dump() << " : " << i.second[0]
		          << std::endl;
}

template <typename T> void scenario1() {
	MGRN<T> g0;
	g0.params = {{12.0, 1.0}};
	g0.addProtein(ProteinType::input, "A", {{{0.69, 0.1, 0.2}}, 0.5, true, false, false});
	g0.addProtein({{{0.15, 0.8, 0.9}}, 0.5, false, false, false});  // B
	g0.addProtein({{{0.5, 0.5, 0.5}}, 0.5, false, false, false});   // C
	g0.addProtein(ProteinType::output, "D", {{{0.1, 0.2, 0.3}}, 0.5, false, true, false});
	MGRN<T> g1tmp;
	g1tmp.params = {{7.0, 1.2}};
	g1tmp.addProtein({{{0.7, 0.2, 0.33}}, 0.5, true, false, false});    // E
	g1tmp.addProtein({{{0.9, 0.1, 0.65}}, 0.5, false, true, false});    // F
	g1tmp.addProtein({{{0.75, 0.12, 0.3}}, 0.5, false, false, false});  // G
	MGRN<T> g2tmp;
	g2tmp.params = {{12.5, 1.5}};
	g2tmp.addProtein({{{0.48, 0.95, 0.15}}, 0.5, true, true, false});  // H
	MGRN<T> g3tmp;
	g3tmp.params = {{19.0, 0.6}};
	g3tmp.addProtein({{{0.1, 0.87, 0.25}}, 0.5, true, true, false});  // H

	g1tmp.addSubNet(g2tmp);
	g1tmp.addSubNet(g3tmp);
	g0.addSubNet(g1tmp);

	REQUIRE(g0.isMaster());
	REQUIRE(g0.subNets.size() == 1);
	auto &g1 = g0.subNets[0];
	REQUIRE(g1.subNets.size() == 2);
	auto &g2 = g1.subNets[0];
	auto &g3 = g1.subNets[1];
	REQUIRE(!g1.isMaster());
	REQUIRE(!g2.isMaster());
	REQUIRE(!g3.isMaster());
	REQUIRE(g1.master == &g0);
	REQUIRE(g2.master == &g0);
	REQUIRE(g3.master == &g0);

	// signatures
	const double eps = 0.000001;

	REQUIRE(g0.signatures.size() == 3);

	auto *A = &g0.actualProteins[0];
	auto *B = &g0.actualProteins[1];
	auto *C = &g0.actualProteins[2];
	auto *D = &g0.actualProteins[3];
	auto *E = &g1.actualProteins[0];
	auto *F = &g1.actualProteins[1];
	auto *G = &g1.actualProteins[2];
	auto *H = &g2.actualProteins[0];
	auto *I = &g3.actualProteins[0];

	std::cerr << "A = " << A->toJSON().dump() << std::endl;
	std::cerr << "B = " << B->toJSON().dump() << std::endl;
	std::cerr << "C = " << C->toJSON().dump() << std::endl;
	std::cerr << "D = " << D->toJSON().dump() << std::endl;
	std::cerr << "E = " << E->toJSON().dump() << std::endl;
	std::cerr << "F = " << F->toJSON().dump() << std::endl;
	std::cerr << "G = " << G->toJSON().dump() << std::endl;
	std::cerr << "H = " << H->toJSON().dump() << std::endl;
	std::cerr << "I = " << I->toJSON().dump() << std::endl;

	auto sA = getInfluencesTo(A, g0);
	auto sB = getInfluencesTo(B, g0);
	auto sC = getInfluencesTo(C, g0);
	auto sD = getInfluencesTo(D, g0);
	auto sE = getInfluencesTo(E, g0);
	auto sF = getInfluencesTo(F, g0);
	auto sG = getInfluencesTo(G, g0);
	auto sH = getInfluencesTo(H, g0);
	auto sI = getInfluencesTo(I, g0);

	std::cerr << "sB = " << std::endl;
	printInfluences(sB);
	REQUIRE(sA.size() == 0);
	REQUIRE(sB.size() == 4);
	REQUIRE(sC.size() == 4);
	REQUIRE(sD.size() == 4);
	REQUIRE(sE.size() == 3);
	REQUIRE(sF.size() == 4);
	REQUIRE(sG.size() == 4);
	REQUIRE(sH.size() == 3);

	// REQUIRE(sBA.second[3] == Approx(0.548811636094027).epsilon(eps));
}

TEST_CASE("Scenarios", "[mgrn]") { scenario1<MGClassic>(); }
