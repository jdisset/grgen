#ifndef MGRN_HPP
#define MGRN_HPP
#include <assert.h>
#include <array>
#include <map>
#include <unordered_map>
#include <string>
#include <sstream>
#include "common.h"

using std::array;
using std::vector;
using std::map;
using std::string;
using std::pair;
using std::ostringstream;
template <typename Implem> struct MGRN {
	struct GAConfiguration {
		// crossover
		static constexpr double ALIGN_TRESHOLD = 0.5;
		static constexpr double APPEND_NON_ALIGNED = 0.2;
		static constexpr unsigned int MAX_REGULS = 40;

		// proteins
		double MODIF_PROT_RATE = 10;
		double ADD_PROT_RATE = 2;
		double DEL_PROT_RATE = 2;
		double ADD_GRN_RATE = 0.5;
		double DEL_GRN_RATE = 0.5;
	};

	friend Implem;
	using Protein = typename Implem::Protein_t;
	using json = nlohmann::json;
	using InfluenceVec = array<double, Implem::nbSignatureParams>;
	template <typename A, typename B> using umap = std::unordered_map<A, B>;
	map<string, size_t>> inputProteins, outputProteins;
	array<double, Implem::nbParams> params{};  // alpha, beta, ...
	vector<pair<Protein*, bool>> allProteinsPtr;
	vector<Protein> actualProteins;
	vector<MGRN> subNets;
	vector<vector<InfluenceVec>>
	    signatures;  // stores the influence of one protein onto the others (p0
	MGRN* master = this;
	MGRN* parent = nullptr;
	Implem implem;

	MGRN(const MGRN& grn, MGRN* m = this, MGRN* p = nullptr)
	    : params(grn.params), actualProteins(grn.actualProteins), master(m), parent(p) {
		for (auto& g : grn.subNets) {
			subNets.push_back(GRN(g, m, this));
		}
		updateAllProteinsPtr();
		updateSignatures();
	}

	/**************************************
	 *               GET
	 *************************************/
	inline double getProteinConcentration(const string& name, const ProteinType t) const {
		try {
			return actualProteins[proteinsRefs[to_underlying(t)].at(name)].c;
		} catch (...) {
			std::cerr << "Exception raised in getProteinConcentration for name = " << name
			          << ", proteintype = " << to_underlying(t) << std::endl;
			exit(0);
		}
	}

	size_t getNbProteins() const { return internalGRN.proteinsSize(); }

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
	inline void reset() { internalGRN.reset(); }

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
			auto protid = internalGRN.addProtein(p);
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
			aId.remove(best.first);
			bId.remove(best.second);
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
		for (auto& a : get<0>(alignedProts)) dProt += get<2>(a);
		double nbNotAligned = get<1>(alignedProts).size() + get<2>(alignedProts).size();
		dProt += nbNotAligned;
		double divisor = nbNotAligned + (double)get<0>(alignedProts).size();
		if (divisor > 0) dProt = dProt / divisor;
		if (subNets.size() == 0) {
			return dProt;
		} else {
			double dG = 0.0;
			auto alignedGrns = getAligned(a.subNets, b.subNets);
			for (auto& a : get<0>(alignedGrns)) dG += get<2>(a);
			double nbNotAligned = get<1>(alignedGrns).size() + get<2>(alignedGrns).size();
			dG += nbNotAligned;
			double divGrn = nbNotAligned + (double)get<0>(alignedGrns).size();
			if (divGrn > 0) dG = dG / divGrn;
			return protVsGrn * dProt + (1.0 - protVsGrn) * dG;
		}
	}

	size_t getNbProteins() { return actualProteins.size(); }
	void randomReguls(size_t n) { internalGRN.randomReguls(n); }

	/**************************************
	 *          UPDATES
	 *************************************/
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
		firstRegulPtrIndex = allProteins.size();
		// reguls
		for (auto& p : actualProteins)
			if ((p.input && p.output) || (!p.input && !p.output)) allProteinsPtr.push_back(&p);
		for (auto& g : subNets) {
			for (auto& p : g.actualProteins) {
				// subs outputs+inputs are reguls for us
				if (p.input && p.output) allProteinsPtr.push_back(&p);
			}
		}
		firstOutputPtrIndex = allProteins.size();
		// outputs
		for (auto& g : subNets) {
			for (auto& p : g.actualProteins) {
				// subs inputs are outputs for us
				if (p.input && !p.output) allProteinsPtr.push_back(&p);
			}
		}
		for (auto& p : g.actualProteins)
			if (!p.input && p.output) allProteinsPtr.push_back(&p);
	}

	void updateSignatures() { implem.updateSignatures(*this); }

	vector<vector<InfluenceVec>> getSignatures() { return signatures; }

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

	void randomReguls(size_t n) {}

	size_t pipedRoulette(const std::vector<double>& coefs) {
		assert(coefs.size() > 0);
		double tot = 0.0;
		for (auto& c : coefs) tot += c;
		uniform_real_distribution<double> dice(0, 1);
		size_t i = 0;
		double roll = dice(grnRand);
		double sum = coefs[0];
		while (i < coefs && roll > sum) {
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
		uniform_int_distribution<size_t> dice(0, actualProteins.size() - 1);
		if (!inp) actualProteins[dice(grnRand)].input = true;
		if (!out) actualProteins[dice(grnRand)].output = true;
		parent->updateAllProteinsPtr();
	}

	void mutate() {
		size_t operation = pipedRoulette(
		    {{MODIF_PROT_RATE, ADD_PROT_RATE, DEL_PROT_RATE, ADD_GRN_RATE, DEL_GRN_RATE}});
		switch (operation) {
			case 0: {
				// modif
				auto allP = getListOfAllProteins();
				uniform_int_distribution<size_t> dice(0, allP.size() - 1);
				size_t id = dice(grnRand);
				allP.second->mutate();
				allP.first->enforceOneInputOneOutput();
			} break;
			case 1: {
				// add
				auto allG = getListOfAllGRNs();
				uniform_int_distribution<size_t> dice(0, allG.size() - 1);
				allG[dice(grnRand)]->addRandomProtein();
			} break;
			case 2: {
				// del
				auto allG = getListOfAllGRNs();
				uniform_int_distribution<size_t> diceG(0, allG.size() - 1);
				auto g = allG[diceG(grnRand)];
				uniform_int_distribution<size_t> diceP(0, g->actualProteins.size() - 1);
				deleteProtein(diceP(grnRand));
			} break;
		}
		// modify/add/delete Prot (from all grn)
		// add grn (wherever)
		// delete grn (only a leaf)
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
		for (size_t t = to_underlying(ProteinType::input);
		     t <= to_underlying(ProteinType::output); ++t) {
			json prots = o.at("proteins").at(typeToString((ProteinType)t));
			for (json::iterator it = prots.begin(); it != prots.end(); ++it) {
				addProtein((ProteinType)t, it.key(), Protein(it.value()));
			}
		}
		updateSignatures();
	}

	string toJSON() const {
		json protObj;
		for (size_t t = 0; t < proteinsRefs.size(); ++t) {
			json pr;
			for (auto p = proteinsRefs[t].begin(); p != proteinsRefs[t].end(); ++p) {
				pr[p->first] = actualProteins[p->second].toJSON();
			}
			protObj[typeToString((ProteinType)t)] = pr;
		}
		json o;
		o["proteins"] = protObj;
		o["params"] = params;
		return o.dump(2);
	}
};
/**************************************
 *       MUTATION & CROSSOVER
 *************************************/
void mutate() { internalGRN.mutate() : }

MGRN crossover(const MGRN& other) { return MGRN::crossover(*this, other); }

static MGRN crossover(const MGRN& g0, const MGRN& g1) {
	return MGRN(decltype(internalGRN)::crossover(g0.internalGRN, g1.internalGRN));
}
}
;
#endif
