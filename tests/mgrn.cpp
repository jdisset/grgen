#include "../grn.hpp"
#include "../json/json.hpp"
#include "../mgclassic.hpp"
#include "../mgrn.hpp"
#include "../real.hpp"
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

template <typename G> void checkInOutProteins(G &g) {
	// check inputs & outputs
	if (g.isMaster()) {
		REQUIRE(g.inputProteins.size() > 0);
		REQUIRE(g.outputProteins.size() > 0);
	} else {
		REQUIRE(g.inputProteins.size() == 0);
		REQUIRE(g.outputProteins.size() == 0);
	}
	for (auto &i : g.inputProteins) {
		REQUIRE(g.actualProteins.size() > i.second);
		REQUIRE(!g.actualProteins[i.second].modifiable);
		REQUIRE(g.actualProteins[i.second].input);
	}
	for (auto &i : g.outputProteins) {
		REQUIRE(g.actualProteins.size() > i.second);
		REQUIRE(!g.actualProteins[i.second].modifiable);
		REQUIRE(g.actualProteins[i.second].output);
	}
}

template <typename G> void checkMGRNIntegrity(G &g, G *topLevel) {
	REQUIRE(g.master == topLevel);
	for (auto &sg : g.subNets) REQUIRE(sg.parent == &g);
	atLeastOneInputOneOutput(g);
	checkInOutProteins(g);
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

template <typename T> std::string printSignatures(const T &g, unsigned int padd = 0) {
	std::ostringstream o;
	o << std::endl;
	for (unsigned int i = 0; i < padd; ++i) o << " ";
	o << "|-" << std::endl;
	for (auto &s : g.signatures) {
		for (unsigned int i = 0; i < padd; ++i) o << " ";
		o << s.first << "[ ";
		for (auto &s2 : s.second) {
			o << "{" << s2.first << " -> " << s2.second[0] << "} ";
		}
		o << "]" << std::endl;
	}
	for (auto &s : g.subNets) o << printSignatures(s, padd + 2);
	return o.str();
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
		INFO("grn0 =  " << grn0.serialize());
		INFO("grn1 =  " << grn1.serialize());
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
		INFO("copy & move");
		auto mgcopy = mgrn;
		REQUIRE(mgrn.getListOfAllProteins().size() == mgcopy.getListOfAllProteins().size());
		REQUIRE(mgcopy == mgrn);
		// same signatures as well
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
		INFO("deletion");
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
		INFO("serialization");
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
		INFO("determinism");
		deterministicGRN<T>();
	}

	SECTION("mutation & crossover") {
		INFO("mutation");
		auto otherCopy = mgrn;
		REQUIRE(otherCopy == mgrn);
		REQUIRE(MGRN<T>::relativeDistance(mgrn, otherCopy) == 0);
		double avgD = 0;
		double nbMuts = 300;
		for (double i = 0; i < nbMuts; ++i) {
			auto tmp = mgrn;
			auto mut = mgrn.mutate();
			INFO("mutation type : " << mut);
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

		INFO("crossover");
		int nbCrossovers = 300;
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

template <typename P, typename G>
std::pair<bool, typename G::InfluenceVec> getSignature(P *a, P *b, const G &g) {
	for (const auto &s : g.signatures) {
		if (s.first == a) {
			for (const auto &sa : s.second) {
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
	for (const auto &s : g.signatures)
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
	for (int i = 0; i < 50; ++i) {
		MGRN<T> g0;
		g0.params = {{12.0, 1.0}};
		g0.addProtein(ProteinType::input, "A", {{{0.69, 0.1, 0.2}}, 0.5, true, false, false});
		g0.addProtein(ProteinType::output, "D", {{{0.1, 0.2, 0.3}}, 0.5, false, true, false});
		MGRN<T> g1tmp;
		g1tmp.params = {{7.0, 1.2}};
		g1tmp.addProtein({{{0.7, 0.2, 0.33}}, 0.5, true, false, false});  // E
		g1tmp.addProtein({{{0.9, 0.1, 0.65}}, 0.5, false, true, false});  // F
		g0.addSubNet(g1tmp);
		auto g1 = g0;
		REQUIRE(g0.allProteinsPtr.size() == g1.allProteinsPtr.size());
		REQUIRE(g0.subNets[0].allProteinsPtr.size() == g1.subNets[0].allProteinsPtr.size());
		REQUIRE(g0.signatures.size() == g1.signatures.size());
		for (size_t k = 0; k < g0.signatures.size(); ++k)
			REQUIRE(g0.signatures[k].second.size() == g1.signatures[k].second.size());
		for (size_t k = 0; k < g0.subNets[0].signatures.size(); ++k)
			REQUIRE(g0.subNets[0].signatures[k].second.size() ==
			        g1.subNets[0].signatures[k].second.size());
		g0.reset();
		g1.reset();
		REQUIRE(g0 == g1);
		g0.step(10);
		g1.step(10);
		INFO("G0 signatures = " << printSignatures(g0));
		INFO("G1 signatures = " << printSignatures(g1));
		REQUIRE(g0 == g1);
		auto g2 = g1;
		REQUIRE(g0 == g2);
		g0.step(100);
		g1.step(100);
		g2.step(100);
		REQUIRE(g0 == g1);
		REQUIRE(g0 == g2);
	}
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

	INFO("A = " << A->toJSON().dump());
	INFO("B = " << B->toJSON().dump());
	INFO("C = " << C->toJSON().dump());
	INFO("D = " << D->toJSON().dump());
	INFO("E = " << E->toJSON().dump());
	INFO("F = " << F->toJSON().dump());
	INFO("G = " << G->toJSON().dump());
	INFO("H = " << H->toJSON().dump());
	INFO("I = " << I->toJSON().dump());

	auto sA = getInfluencesTo(A, g0);
	auto sB = getInfluencesTo(B, g0);
	auto sC = getInfluencesTo(C, g0);
	auto sD = getInfluencesTo(D, g0);
	auto sE = getInfluencesTo(E, g1);
	auto sF = getInfluencesTo(F, g1);
	auto sG = getInfluencesTo(G, g1);
	auto sH = getInfluencesTo(H, g2);
	auto sI = getInfluencesTo(I, g3);

	REQUIRE(sA.size() == 0);
	REQUIRE(sB.size() == 4);
	REQUIRE(sC.size() == 4);
	REQUIRE(sD.size() == 4);
	REQUIRE(sE.size() == 3);
	REQUIRE(sF.size() == 4);
	REQUIRE(sG.size() == 4);
	REQUIRE(sH.size() == 4);
	REQUIRE(sI.size() == 4);

	auto sBA = getSignature(B, A, g0);
	REQUIRE(sBA.first);
	auto sBB = getSignature(B, B, g0);
	REQUIRE(sBB.first);
	auto sBC = getSignature(B, C, g0);
	REQUIRE(sBC.first);
	auto sBF = getSignature(B, F, g0);
	REQUIRE(sBF.first);

	auto sCA = getSignature(C, A, g0);
	REQUIRE(sCA.first);
	auto sCB = getSignature(C, B, g0);
	REQUIRE(sCB.first);
	auto sCC = getSignature(C, C, g0);
	REQUIRE(sCC.first);
	auto sCF = getSignature(C, F, g0);
	REQUIRE(sCF.first);

	auto sDA = getSignature(D, A, g0);
	REQUIRE(sDA.first);
	auto sDB = getSignature(D, B, g0);
	REQUIRE(sDB.first);
	auto sDC = getSignature(D, C, g0);
	REQUIRE(sDC.first);
	auto sDF = getSignature(D, F, g0);
	REQUIRE(sDF.first);

	auto sEA = getSignature(E, A, g1);
	REQUIRE(sEA.first);
	auto sEB = getSignature(E, B, g1);
	REQUIRE(sEB.first);
	auto sEC = getSignature(E, C, g1);
	REQUIRE(sEC.first);

	auto sFE = getSignature(F, E, g1);
	REQUIRE(sFE.first);
	auto sFG = getSignature(F, G, g1);
	REQUIRE(sFG.first);
	auto sFH = getSignature(F, H, g1);
	REQUIRE(sFH.first);
	auto sFI = getSignature(F, I, g1);
	REQUIRE(sFI.first);

	auto sGE = getSignature(G, E, g1);
	REQUIRE(sGE.first);
	auto sGG = getSignature(G, G, g1);
	REQUIRE(sGG.first);
	auto sGH = getSignature(G, H, g1);
	REQUIRE(sGH.first);
	auto sGI = getSignature(G, I, g1);
	REQUIRE(sGI.first);

	auto sHE = getSignature(H, E, g2);
	REQUIRE(sHE.first);
	auto sHG = getSignature(H, G, g2);
	REQUIRE(sHG.first);
	auto sHH = getSignature(H, H, g2);
	REQUIRE(sHH.first);
	auto sHI = getSignature(H, I, g2);
	REQUIRE(sHI.first);

	auto sIE = getSignature(I, E, g3);
	REQUIRE(sIE.first);
	auto sIG = getSignature(I, G, g3);
	REQUIRE(sIG.first);
	auto sIH = getSignature(I, H, g3);
	REQUIRE(sIH.first);
	auto sII = getSignature(I, I, g3);
	REQUIRE(sII.first);

	const double eps = 1e-12;
	REQUIRE(sBA.second[0] == Approx(0.548811636094027).epsilon(eps));
	REQUIRE(sBB.second[0] == Approx(0.000409734978979786).epsilon(eps));
	REQUIRE(sBC.second[0] == Approx(0.0149955768204777).epsilon(eps));
	REQUIRE(sBF.second[0] == Approx(0.548811636094027).epsilon(eps));

	REQUIRE(sCA.second[0] == Approx(0.00822974704902002).epsilon(eps));
	REQUIRE(sCB.second[0] == Approx(0.0273237224472925).epsilon(eps));
	REQUIRE(sCC.second[0] == Approx(1).epsilon(eps));
	REQUIRE(sCF.second[0] == Approx(0.00822974704902002).epsilon(eps));

	REQUIRE(sDA.second[0] == Approx(1).epsilon(eps));
	REQUIRE(sDB.second[0] == Approx(0.000224867324178848).epsilon(eps));
	REQUIRE(sDC.second[0] == Approx(0.00822974704902002).epsilon(eps));
	REQUIRE(sDF.second[0] == Approx(1).epsilon(eps));

	REQUIRE(sEA.second[0] == Approx(0.0149955768204777).epsilon(eps));
	REQUIRE(sEB.second[0] == Approx(0.496585303791409).epsilon(eps));
	REQUIRE(sEC.second[0] == Approx(0.246596963941607).epsilon(eps));

	REQUIRE(sFE.second[0] == Approx(0.00744658307092434).epsilon(eps));
	REQUIRE(sFG.second[0] == Approx(0.00425355574481513).epsilon(eps));
	REQUIRE(sFH.second[0] == Approx(0.704688089718714).epsilon(eps));
	REQUIRE(sFI.second[0] == Approx(0.810584245970187).epsilon(eps));

	REQUIRE(sGE.second[0] == Approx(0.0212797364383772).epsilon(eps));
	REQUIRE(sGG.second[0] == Approx(0.0121551783299149).epsilon(eps));
	REQUIRE(sGH.second[0] == Approx(0.246596963941607).epsilon(eps));
	REQUIRE(sGI.second[0] == Approx(0.43171052342908).epsilon(eps));

	REQUIRE(sHE.second[0] == Approx(0.0301973834223185).epsilon(eps));
	REQUIRE(sHG.second[0] == Approx(0.0111089965382423).epsilon(eps));
	REQUIRE(sHH.second[0] == Approx(0.00280879419452551).epsilon(eps));
	REQUIRE(sHI.second[0] == Approx(0.00763509421885996).epsilon(eps));

	REQUIRE(sIE.second[0] == Approx(0.149568619222635).epsilon(eps));
	REQUIRE(sIG.second[0] == Approx(0.683861409212356).epsilon(eps));
	REQUIRE(sIH.second[0] == Approx(0.000000096859922509254).epsilon(eps));
	REQUIRE(sII.second[0] == Approx(0.000000442865378096327).epsilon(eps));
}

template <typename G> std::string printClassicSignatures(const G &g) {
	std::ostringstream os;
	for (const auto &p : g.getSignatures()) {
		os << " -|";
		for (const auto &s : p) {
			os << " {" << s[0] << ", " << s[1] << "} ";
		}
		os << std::endl;
	}
	return os.str();
}

template <typename T> void comparison() {
	MGRN<T> gm;
	GRN<RealCoords> gc;
	gm.params = {{12.0, 1.0}};
	gm.addProtein(ProteinType::input, "I", {{{0.69, 0.1, 0.2}}, 0.5, true, false, false});
	gm.addProtein({{{0.15, 0.8, 0.9}}, 0.5, false, false, false});
	gm.addProtein(ProteinType::output, "O", {{{0.1, 0.2, 0.3}}, 0.5, false, true, false});
	gc.setParam({{12.0, 1.0}});
	gc.addProtein(ProteinType::input, "I", {{{0.69, 0.1, 0.2}}, 0.5});
	gc.addProtein(ProteinType::regul, "R0", {{{0.15, 0.8, 0.9}}, 0.5});
	gc.addProtein(ProteinType::output, "O", {{{0.1, 0.2, 0.3}}, 0.5});
	gm.reset();
	gc.reset();

	REQUIRE(gm.getProteinConcentration("O", ProteinType::output) ==
	        gc.getProteinConcentration("O", ProteinType::output));
	gm.step(1);
	gc.step(1);
	INFO("gm s = " << printSignatures(gm));
	INFO("gc s = " << printClassicSignatures(gc));

	REQUIRE(gm.getProteinConcentration("O", ProteinType::output) ==
	        gc.getProteinConcentration("O", ProteinType::output));
}

// TEST_CASE("Scenarios", "[mgrn]") { scenario1<MGClassic>(); }
//TEST_CASE("Dynamique", "[mgrn]") { comparison<MGClassic>(); }
TEST_CASE("MGRN random declaration, init & serialization", "[mgrn]") {
	for (int i = 0; i < 100; ++i) testMGRN<MGClassic>();
}
