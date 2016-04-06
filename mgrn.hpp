#ifndef MGRN_HPP
#define MGRN_HPP
#include <assert.h>
#include <array>
#include <map>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include "common.h"

using std::array;
using std::vector;
using std::map;
using std::string;
using std::pair;
using std::ostringstream;
using std::tuple;
template <typename Implem> struct MGRN {
	friend Implem;
	using Protein = typename Implem::Protein_t;
	using json = nlohmann::json;
	using InfluenceVec = array<double, Implem::nbSignatureParams>;
	using signature_t =
	    vector<std::pair<Protein*, vector<std::pair<Protein*, InfluenceVec>>>>;
	template <typename A, typename B> using umap = std::unordered_map<A, B>;
	map<string, size_t> inputProteins, outputProteins;
	array<double, Implem::nbParams> params{};  // alpha, beta, ...
	vector<pair<Protein*, bool>> allProteinsPtr;
	size_t firstRegulPtrIndex = 0;
	size_t firstOutputPtrIndex = 0;
	vector<Protein> actualProteins;
	vector<MGRN> subNets;
	signature_t signatures;  // stores the influence of one protein onto the others (p0
	MGRN* master = this;
	MGRN* parent = nullptr;
	Implem implem;

	// mutation & crossover params
	// mutation
	double MODIF_PROT_RATE = 10;
	double ADD_PROT_RATE = 2;
	double DEL_PROT_RATE = 2;
	double ADD_GRN_RATE = 0.5;
	double DEL_GRN_RATE = 0.5;
	// crossover
	double ALIGN_TRESHOLD = 0.5;
	double APPEND_NON_ALIGNED = 0.2;
	unsigned int MAX_REGULS = 40;

	MGRN(const MGRN& grn, MGRN* m = nullptr, MGRN* p = nullptr)
	    : params(grn.params), actualProteins(grn.actualProteins), master(m), parent(p) {
		if (master == nullptr) master = this;
		for (auto& g : grn.subNets) {
			subNets.push_back(GRN(g, m, this));
		}
		updateAllProteinsPtr();
		updateSignatures();
	}

	/**************************************
	 *               GET
	 *************************************/
	size_t getNbProteins() const { return actualProteins.size(); }

	Protein& getProtein(ProteinType t, const string& name) {
		if (master != this) throw std::invalid_argument("getProtein called on a subGRN");
		if (t == ProteinType::input) return *inputProteins.at(name);
		if (t == ProteinType::output) return *outputProteins.at(name);
		throw std::invalid_argument("MGRN only has named input & output");
	}

	Protein getProtein_const(ProteinType t, const string& name) const {
		if (master != this) throw std::invalid_argument("getProtein called on a subGRN");
		if (t == ProteinType::input) return *inputProteins.at(name);
		if (t == ProteinType::output) return *outputProteins.at(name);
		throw std::invalid_argument("MGRN only has named input & output");
	}

	/**************************************
	 *               SET
	 *************************************/

	void setProteinConcentration(const string& name, ProteinType t, double c) {
		if (master != this) throw std::invalid_argument("setProteinConc on a subGRN");
		if (t == ProteinType::input) inputProteins.at(name)->setConcentration(c);
		if (t == ProteinType::output) outputProteins.at(name)->setConcentration(c);
		throw std::invalid_argument("MGRN only has named inputs & outputs");
	}

	vector<string> getProteinNames(ProteinType t) const {
		if (master != this) throw std::invalid_argument("getProteinNames on a subGRN");
		vector<string> res;
		if (t == ProteinType::output)
			for (auto& p : outputProteins) res.push_back(p.first);
		else if (t == ProteinType::input)
			for (auto& p : inputProteins) res.push_back(p.first);
		else
			throw std::invalid_argument("MGRN only has named inputs & outputs");
		return res;
	}

	/**************************************
	 *          ADDING PROTEINS
	 *************************************/
	size_t addProtein(const Protein& p) {
		actualProteins.push_back(p);
		updateAllProteinsPtr();
		return actualProteins.size() - 1;
	}

	void addProtein(ProteinType t, const string& name, const Protein& p) {
		if (master != this)
			throw std::invalid_argument("named addProtein on a subGRN");
		else if (t != ProteinType::input && t != ProteinType::output)
			throw std::invalid_argument("MGRN only has named inputs & outputs");
		else {
			auto protid = addProtein(p);
			if (t == ProteinType::input)
				inputProteins[name] = protid;
			else if (t == ProteinType::output)
				outputProteins[name] = protid;
		}
	}

	void addRandomProtein(ProteinType t, const string& name, bool m = false) {
		auto p = Protein();
		p.input = false;
		p.output = false;
		p.modifiableFlags = m;
		if (t == ProteinType::input)
			p.input = true;
		else if (t == ProteinType::output)
			p.output = true;
		addProtein(t, name, p);
	}

	void addRandomProtein(size_t n = 1) {
		for (size_t i = 0; i < n; ++n) {
			addProtein(Protein());
		}
	}

	void randomReguls(size_t n) { addRandomProtein(n); }

	template <typename T>
	static tuple<vector<tuple<size_t, size_t, double>>, vector<size_t>, vector<size_t>>
	    getAligned(const vector<T>& a, const vector<T>& b, double threshold = 1.0) {
		// returns tuple of aligned id {a,b,dist} + a's rest + b's rest
		vector<size_t> aId, bId;
		aId.reserve(a.size());
		bId.reserve(b.size());
		vector<tuple<size_t, size_t, double>> aligned;
		for (size_t i = 0; i < a.size(); ++i) aId.push_back(i);
		for (size_t i = 0; i < b.size(); ++i) bId.push_back(i);
		double lastTr = 0.0;
		while (lastTr <= threshold && aId.size() && bId.size()) {
			pair<size_t, size_t> best = {0, 0};
			double minDist = T::relativeDistance(a[aId[best.first]], b[bId[best.second]]);
			for (size_t i = 1; i < aId.size(); ++i) {
				for (size_t j = 1; j < bId.size(); ++j) {
					double dist = T::relativeDistance(a[aId[i]], b[bId[j]]);
					if (dist < minDist) {
						minDist = dist;
						best = {i, j};
					}
				}
			}
			aligned.push_back({aId[best.first], bId[best.second], minDist});
			aId.erase(std::begin(aId) + static_cast<long>(best.first));
			bId.erase(std::begin(bId) + static_cast<long>(best.second));
		}
		return {aligned, aId, bId};
	}

	static double relativeDistance(const MGRN& a, const MGRN& b) {
		const double protVsGrn = 0.4;
		// first we align proteins
		// TODO : better distance between proteins (differences in the dynamics of the whole
		// set, not just in the absolute coords of each individual protein)
		// ==> just try with all possible id offsets and find the minimal distance
		// => pass a lambda to getAligned so we can use custom distance methods
		// and define a distance btwn prots that applies an offset to the 2nd one)
		// except if it is an input / output ? or maybe use allProteinsPtr?
		double dProt = 0.0;
		auto alignedProts = getAligned(a.actualProteins, b.actualProteins);
		for (auto& al : std::get<0>(alignedProts)) dProt += std::get<2>(al);
		double nbNotAligned =
		    std::get<1>(alignedProts).size() + std::get<2>(alignedProts).size();
		dProt += nbNotAligned;
		double divisor = nbNotAligned + std::get<0>(alignedProts).size();
		if (divisor > 0) dProt = dProt / divisor;
		if (a.subNets.size() + b.subNets.size() == 0) {
			return dProt;
		} else {
			double dG = 0.0;
			auto alignedGrns = getAligned(a.subNets, b.subNets);
			for (auto& al : std::get<0>(alignedGrns)) dG += std::get<2>(al);
			double nbNotAlignedGrn =
			    std::get<1>(alignedGrns).size() + std::get<2>(alignedGrns).size();
			dG += nbNotAlignedGrn;
			double divGrn = nbNotAlignedGrn + (double)std::get<0>(alignedGrns).size();
			if (divGrn > 0) dG = dG / divGrn;
			return protVsGrn * dProt + (1.0 - protVsGrn) * dG;
		}
	}

	size_t getNbProteins() { return actualProteins.size(); }
	size_t getNbOwnProteins() { return actualProteins.size(); }
	size_t getNbOwnProteins(const ProteinType& t) {
		switch (t) {
			case (ProteinType::input):
				return firstRegulPtrIndex;
			case (ProteinType::regul):
				return firstOutputPtrIndex - firstRegulPtrIndex;
			default:
				return actualProteins.size() - firstOutputPtrIndex;
		}
	}

	/**************************************
	 *          UPDATES
	 *************************************/
	void updateSignatures() {
		/*
		When step is called, a grn needs to update the concentration of every of its actual
		Proteins. Thus, updateSignatures needs to generate, for every protein* p0 needing to
		be updated by this grn, the list of every protein* p1 (including itself) influencing
		p0 (p1 might come from a parent or a sub net), as well as the corresponding
		influenceVec (the actual signature values, e.g.  enhance, inhibit...). The signatures
		thus tell which proteins should be updated and how.
		*/
		signatures.clear();
		signatures.resize(allProteinsPtr.size());
		for (size_t i = 0; i < allProteinsPtr.size(); ++i) {
			if (allProteinsPtr[i].second) {
				auto& p = allProteinsPtr[i];
				if (parent || !p.input) {
					// same grn protein
					// its influence pool is all of the input + reguls proteins of this grn
					vector<std::pair<Protein*, InfluenceVec>> influencePool;
					influencePool.reserve(allProteinsPtr.size());
					for (size_t j = 0; j < firstOutputPtrIndex; ++j) {
						auto& p1 = allProteinsPtr[j];
						influencePool.push_back({p1, Implem::getInfluenceVec(p, p1, *this)});
					}
					if (parent && (p.input)) {
						// p is also influenced by parent's proteins
						for (size_t j = 0; j < parent->firstOutputPtrIndex; ++j) {
							auto& p1 = parent->allProteinsPtr[j];
							influencePool.push_back({p1, Implem::getInfluenceVec(p, p1, *parent)});
						}
					}
					signatures.push_back({p, influencePool});
				}
			}
		}
		for (auto& sg : subNets) sg.updateSignatures();
	}
	void updateAllProteinsPtr() {
		allProteinsPtr = vector<Protein*>();
		// we put inputs at the beginning
		for (auto& g : subNets) {
			for (auto& p : g.actualProteins) {
				// subs outputs are inputs for us
				if (p.output && !p.input) allProteinsPtr.push_back(&p);
			}
		}
		for (auto& p : actualProteins)
			if (p.input && !p.output) allProteinsPtr.push_back(&p);
		firstRegulPtrIndex = allProteinsPtr.size();
		// reguls
		for (auto& p : actualProteins)
			if ((p.input && p.output) || (!p.input && !p.output)) allProteinsPtr.push_back(&p);
		for (auto& g : subNets) {
			for (auto& p : g.actualProteins) {
				// subs outputs+inputs are reguls for us
				if (p.input && p.output) allProteinsPtr.push_back(&p);
			}
		}
		firstOutputPtrIndex = allProteinsPtr.size();
		// outputs
		for (auto& g : subNets) {
			for (auto& p : g.actualProteins) {
				// subs inputs are outputs for us
				if (p.input && !p.output) allProteinsPtr.push_back(&p);
			}
		}
	}

	signature_t getSignatures() { return signatures; }

	void step(unsigned int nbSteps = 1) { implem.step(*this, nbSteps); }

	inline void reset() {
		for (auto& p : actualProteins) p.reset();
		for (auto& g : subNets) g.reset();
	}

	void setParam(size_t i, double val) {
		if (i < params.size()) params[i] = val;
	}
	void randomParams() {
		array<pair<double, double>, Implem::nbParams> limits = Implem::paramsLimits();
		for (size_t i = 0; i < Implem::nbParams; ++i) {
			std::uniform_real_distribution<double> distrib(limits[i].first, limits[i].second);
			params[i] = distrib(grnRand);
		}
		updateSignatures();
	}

	size_t pipedRoulette(const std::vector<double>& rawCoefs) {
		auto coefs = rawCoefs;
		assert(coefs.size() > 0);
		double tot = 0.0;
		for (auto& c : coefs) tot += c;
		assert(tot > 0.0);
		for (auto& c : coefs) c /= tot;
		std::uniform_real_distribution<double> dice(0, 1);
		size_t i = 0;
		double roll = dice(grnRand);
		double sum = coefs[0];
		while (i < coefs.size() && roll > sum) {
			++i;
			sum += coefs[i];
		}
		return i;
	}

	void enforceOneInputOneOutput() {
		bool inp = false;
		bool out = false;
		for (auto& p : actualProteins) {
			if (p.input) inp = true;
			if (p.output) out = true;
		}
		std::uniform_int_distribution<size_t> dice(0, actualProteins.size() - 1);
		if (!inp) actualProteins[dice(grnRand)].input = true;
		if (!out) actualProteins[dice(grnRand)].output = true;
		parent->updateAllProteinsPtr();
	}

	std::vector<MGRN<Implem>*> getListOfAllGRNs() {
		// returns all GRNS present in the whole network
		std::unordered_set<MGRN<Implem>*> visited;
		std::vector<MGRN<Implem>*> toVisit;
		toVisit.insert(master);
		while (toVisit.size()) {
			auto g = toVisit.back();
			toVisit.pop_back();
			visited.insert(g);
			for (auto& sg : g->subNets)
				if (!visited.count(&sg)) toVisit.push_back(&sg);
		}
		std::vector<MGRN<Implem>*> result;
		result.resize(visited.size());
		for (auto& g : visited) result.push_back(g);
		return result;
	}
	std::vector<std::pair<MGRN<Implem>*, Protein*>> getListOfAllProteins() {
		// returns all of the proteins present in the whole network
		std::vector<std::pair<MGRN<Implem>*, Protein*>> result;
		auto& allGrns = getListOfAllGRNs();
		for (auto& g : allGrns)
			for (auto& p : g->actualProteins) result.push_back({g, &p});
		return result;
	}

	void deleteSubNet(MGRN<Implem>* g) {
		size_t i;
		for (i = 0; i < subNets.size(); ++i)
			if (&subNets[i] == g) break;
		deleteSubNet(i);
	}
	void deleteSubNet(size_t i) {
		assert(i < subNets.size());
		subNets.erase(subNets.begin() + i);
	}
	void deleteProtein(size_t i) {
		assert(i < actualProteins.size());
		actualProteins.erase(actualProteins.begin() + i);
	}

	/**************************************
	     *              JSON
	*************************************/
	MGRN(const string& js) {
		auto o = json::parse(js);
		assert(o.count("params"));
		json par = o.at("params");
		assert(par.size() == Implem::nbParams);
		size_t i = 0;
		for (auto& p : par) params[i++] = p.get<double>();
		assert(o.count("proteins"));
		auto& pobj = o.at("proteins");
		if (pobj.count("namedIn"))
			inputProteins = pobj.at("namedIn").get<decltype(inputProteins)>();
		if (pobj.count("namedOut"))
			inputProteins = pobj.at("namedOut").get<decltype(inputProteins)>();
		assert(pobj.count("plist"));
		actualProteins.resize(pobj.at("plist").size());
		for (auto& p : pobj.at("plist")) actualProteins.push_back(Protein(p));
		updateAllProteinsPtr();
		updateSignatures();
	}

	string toJSON() const {
		json protObj;
		if (inputProteins.size()) protObj["namedIn"] = inputProteins;
		if (outputProteins.size()) protObj["namedOut"] = outputProteins;
		protObj["plist"] = json::array();
		for (auto& p : actualProteins) protObj["plist"].push_back(p.toJSON());
		json o;
		o["params"] = params;
		o["proteins"] = protObj;
		if (subNets.size()) {
			json grnObj;
			for (auto& g : subNets) grnObj.push_back(g.toJSON());
			o["subnets"] = grnObj;
		}
		return o.dump(2);
	}

	/**************************************
	 *       MUTATION & CROSSOVER
	 *************************************/

	void mutate() {
		size_t operation = pipedRoulette(
		    {{MODIF_PROT_RATE, ADD_PROT_RATE, DEL_PROT_RATE, ADD_GRN_RATE, DEL_GRN_RATE}});
		switch (operation) {
			case 0: {
				// modif prot
				auto allP = getListOfAllProteins();
				std::uniform_int_distribution<size_t> dice(0, allP.size() - 1);
				size_t id = dice(grnRand);
				allP[id].second->mutate();
				allP[id].first->enforceOneInputOneOutput();
			} break;
			case 1: {
				// add prot
				auto allG = getListOfAllGRNs();
				std::uniform_int_distribution<size_t> dice(0, allG.size() - 1);
				allG[dice(grnRand)]->addRandomProtein();
			} break;
			case 2: {
				// del prot
				auto allG = getListOfAllGRNs();
				std::uniform_int_distribution<size_t> diceG(0, allG.size() - 1);
				auto g = allG[diceG(grnRand)];
				if (g->actualProteins.size() > 1) {
					std::uniform_int_distribution<size_t> diceP(0, g->actualProteins.size() - 1);
					auto p = g->actualProteins[diceP(grnRand)];
					if (p->modifiable) {
						deleteProtein(diceP(grnRand));
						g->enforceOneInputOneOutput();
					}
				}
			} break;
			case 3: {
				// add grn (wherever)
				auto allG = getListOfAllGRNs();
				std::uniform_int_distribution<size_t> diceG(0, allG.size() - 1);
				auto g = allG[diceG(grnRand)];
				MGRN<Implem> ng;
				ng.addRandomProtein();
				ng.enforceOneInputOneOutput();
				g.addGrn(ng);
			} break;
			case 4: {
				// delete grn (only a leaf)
				auto allG = getListOfAllGRNs();
				if (allG.size() > 1) {
					vector<MGRN<Implem>*> leaves;
					for (auto& g : allG)
						if (g->subNets.size() == 0) leaves.push_back(g);
					std::uniform_int_distribution<size_t> diceG(0, leaves.size() - 1);
					auto g = leaves[diceG(grnRand)];
					g->parent->deleteSubNet(g);
				}
			} break;
		}
	}

	MGRN<Implem> crossover(const MGRN<Implem>& other) {
		return MGRN::crossover(*this, other);
	}

	static MGRN<Implem> crossover(const MGRN<Implem>& g0, const MGRN<Implem>&) {
		return g0;
	};
};
#endif
