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
	MGRN* master = nullptr;
	MGRN* parent = nullptr;
	Implem implem;

	// mutation & crossover params
	// mutation
	double MODIF_PROT_RATE = 10;
	double ADD_PROT_RATE = 1;
	double DEL_PROT_RATE = 1;
	double ADD_GRN_RATE = 0.1;
	double DEL_GRN_RATE = 0.05;
	// crossover
	double ALIGN_TRESHOLD = 0.5;
	double APPEND_NON_ALIGNED = 0.2;
	unsigned int MAX_REGULS = 20;

	MGRN() : master(this), parent(nullptr) {}
	MGRN& operator=(const MGRN& g) {
		inputProteins = g.inputProteins;
		outputProteins = g.outputProteins;
		params = g.params;
		actualProteins = g.actualProteins;
		master = this;
		parent = nullptr;
		subNets.clear();
		signatures.clear();
		allProteinsPtr.clear();
		for (auto& sn : g.subNets) subNets.push_back(MGRN(sn, master));
		updateSubNetsPtrsAndSignatures();
		return *this;
	}

	MGRN(const MGRN& grn, MGRN* m = nullptr)
	    : inputProteins(grn.inputProteins),
	      outputProteins(grn.outputProteins),
	      params(grn.params),
	      actualProteins(grn.actualProteins),
	      master(m ? m : this),
	      parent(nullptr) {
		for (auto& g : grn.subNets) subNets.push_back(MGRN(g, master));
		master->updateSubNetsPtrsAndSignatures();
	}

	explicit MGRN(const string& js) : MGRN(json::parse(js)) {
		master = this;
		parent = nullptr;
		master->updateSubNetsPtrsAndSignatures();
	}

	explicit MGRN(const json& o) {
		assert(o.count("params"));
		json par = o.at("params");
		assert(par.size() == Implem::nbParams);
		size_t i = 0;
		for (auto& p : par) params[i++] = p.get<double>();
		assert(o.count("proteins"));
		auto& pobj = o.at("proteins");
		if (pobj.count("namedIn") > 0)
			inputProteins = pobj.at("namedIn").get<decltype(inputProteins)>();
		if (pobj.count("namedOut") > 0)
			outputProteins = pobj.at("namedOut").get<decltype(outputProteins)>();
		assert(pobj.count("plist"));
		actualProteins.reserve(pobj.at("plist").size());
		for (auto& p : pobj.at("plist")) actualProteins.push_back(Protein(p));
		if (o.count("subnets") > 0)
			for (auto& g : o.at("subnets")) subNets.push_back(MGRN(g));
	}

	/**************************************
	 *               GET
	 *************************************/
	size_t getNbProteins() const { return actualProteins.size(); }

	Protein& getProtein(ProteinType t, const string& name) {
		if (master != this) throw std::invalid_argument("getProtein called on a subGRN");
		if (t == ProteinType::input) return actualProteins[inputProteins.at(name)];
		if (t == ProteinType::output) return actualProteins[outputProteins.at(name)];
		throw std::invalid_argument("MGRN only has named input & output");
	}

	Protein getProtein_const(ProteinType t, const string& name) const {
		if (master != this) throw std::invalid_argument("getProtein called on a subGRN");
		if (t == ProteinType::input) return actualProteins[inputProteins.at(name)];
		if (t == ProteinType::output) return actualProteins[outputProteins.at(name)];
		throw std::invalid_argument("MGRN only has named input & output");
	}

	double getProteinConcentration(const string& name, ProteinType t) const {
		if (master != this) throw std::invalid_argument("getProtein called on a subGRN");
		if (t == ProteinType::input) return actualProteins[inputProteins.at(name)].c;
		if (t == ProteinType::output) return actualProteins[outputProteins.at(name)].c;
		throw std::invalid_argument("MGRN only has named input & output");
		return 0.0;
	}

	/**************************************
	 *               SET
	 *************************************/

	void setProteinConcentration(const string& name, ProteinType t, double c) {
		if (master != this) throw std::invalid_argument("setProteinConc on a subGRN");
		if (t == ProteinType::input)
			actualProteins[inputProteins.at(name)].setConcentration(c);
		else if (t == ProteinType::output)
			actualProteins[outputProteins.at(name)].setConcentration(c);
		else
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
		if (isMaster() && (p.input || p.output)) actualProteins.back().modifiable = false;
		master->updateSubNetsPtrsAndSignatures();
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
		p.modifiable = m;
		if (t == ProteinType::input)
			p.input = true;
		else if (t == ProteinType::output)
			p.output = true;
		addProtein(t, name, p);
	}

	void addRandomProtein(size_t n = 1) {
		for (size_t i = 0; i < n; ++i) {
			auto p = Protein();
			if (master == this) {
				p.input = false;
				p.output = false;
			}
			addProtein(p);
		}
	}

	void randomReguls(size_t n) { addRandomProtein(n); }

	/**************************************
	 *          ADDING SUBNETS
	 *************************************/
	void updateSubNetsPtrsAndSignatures() {
		for (auto& sn : subNets) {
			sn.master = master;
			sn.parent = this;
			sn.updateSubNetsPtrsAndSignatures();
		}
		updateAllProteinsPtr();
		updateSignatures();
	}
	void addSubNet(const MGRN<Implem>& mg) {
		subNets.push_back(mg);
		subNets[subNets.size() - 1].inputProteins.clear();
		subNets[subNets.size() - 1].outputProteins.clear();
		master->updateSubNetsPtrsAndSignatures();
	}

	void addRandomSubNet(size_t nbP = 1) {
		MGRN<Implem> rdmNet;
		if (nbP > 0) {
			rdmNet.addRandomProtein(nbP);
			rdmNet.enforceOneInputOneOutput();
		}
		addSubNet(rdmNet);
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
			aligned.push_back(
			    tuple<size_t, size_t, double>(aId[best.first], bId[best.second], minDist));
			aId.erase(std::begin(aId) + static_cast<int64_t>(best.first));
			bId.erase(std::begin(bId) + static_cast<int64_t>(best.second));
		}
		return tuple<vector<tuple<size_t, size_t, double>>, vector<size_t>, vector<size_t>>(
		    aligned, aId, bId);
	}

	static double relativeDistance(const MGRN& a, const MGRN& b) {
		const double protVsGrn = 0.4;
		// first we align proteins
		// TODO(jean) : better distance between proteins (differences in the dynamics of the
		// whole
		// set, not just in the absolute coords of each individual protein)
		// ==> just try with all possible id offsets and find the minimal distance
		// => pass a lambda to getAligned so we can use custom distance methods
		// and define a distance btwn prots that applies an offset to the 2nd one)
		// except if it is an input / output ? or maybe use allProteinsPtr?
		double dProt = 0.0;
		auto alignedProts = getAligned(a.actualProteins, b.actualProteins);
		for (auto& al : std::get<0>(alignedProts)) dProt += std::get<2>(al);
		double nbNotAligned = static_cast<double>(std::get<1>(alignedProts).size() +
		                                          std::get<2>(alignedProts).size());
		dProt += nbNotAligned;
		double divisor = nbNotAligned + static_cast<double>(std::get<0>(alignedProts).size());
		if (divisor > 0) dProt = dProt / divisor;
		if (a.subNets.size() + b.subNets.size() == 0) return dProt;

		double dG = 0.0;
		auto alignedGrns = getAligned(a.subNets, b.subNets);
		for (auto& al : std::get<0>(alignedGrns)) dG += std::get<2>(al);
		double nbNotAlignedGrn = static_cast<double>(std::get<1>(alignedGrns).size() +
		                                             std::get<2>(alignedGrns).size());
		dG += nbNotAlignedGrn;
		double divGrn = nbNotAlignedGrn + std::get<0>(alignedGrns).size();
		if (divGrn > 0) dG = dG / divGrn;
		return protVsGrn * dProt + (1.0 - protVsGrn) * dG;
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
	bool influencesOthers(std::pair<Protein*, bool> p) {
		return (((!p.first->input && !p.first->output) ||
		         (p.first->input && p.first->output)) &&
		        p.second) ||                    // regul
		       (p.first->input && p.second) ||  // input
		       (p.first->output && !p.second);  // child's output
	}
	template <typename T, typename U>
	bool isInVectorOfPair(const std::vector<std::pair<T, U>>& v, const T& t) {
		for (auto& p : v)
			if (t == p.first) return true;
		return false;
	}
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
		for (auto& p : actualProteins) {
			if (!isMaster() || !p.input) {
				// p needs to be updated and thus needs an influence pool
				std::vector<std::pair<Protein*, InfluenceVec>> inflPool;
				// 1st case: it is influenced by the parent
				if (!isMaster() && p.input)
					for (auto& p1 : parent->allProteinsPtr) {
						if (influencesOthers(p1) && !isInVectorOfPair(inflPool, p1.first))
							inflPool.push_back(
							    {p1.first, Implem::getInfluenceVec(&p, p1.first, *this)});
					}
				// 2nd: it is influenced by the current grn
				if (p.output || (!p.input && !p.output)) {  // if output or regul
					for (auto& p1 : allProteinsPtr) {
						if (influencesOthers(p1) && !isInVectorOfPair(inflPool, p1.first)) {
							inflPool.push_back(
							    {p1.first, Implem::getInfluenceVec(&p, p1.first, *this)});
						}
					}
				}
				signatures.push_back({&p, inflPool});
			}
		}
	}

	void updateAllProteinsPtr() {
		// allProteinsPtr contains all of the actual proteins + any other protein that can
		// influence a protein at this level (using the params of this grn)
		allProteinsPtr.clear();
		// we put inputs at the beginning
		for (auto& g : subNets) {
			for (auto& p : g.actualProteins) {
				// subs outputs are inputs for us
				if (!p.input && p.output) allProteinsPtr.push_back({&p, false});
			}
		}
		// inputs
		for (auto& p : actualProteins)
			if (p.input && !p.output) allProteinsPtr.push_back({&p, true});
		firstRegulPtrIndex = allProteinsPtr.size();
		// reguls
		for (auto& p : actualProteins)
			if ((p.input && p.output) || (!p.input && !p.output))
				allProteinsPtr.push_back({&p, true});
		for (auto& g : subNets) {
			for (auto& p : g.actualProteins) {
				// subs outputs+inputs are reguls for us
				if (p.input && p.output) allProteinsPtr.push_back({&p, false});
			}
		}
		firstOutputPtrIndex = allProteinsPtr.size();
		// outputs
		for (auto& p : actualProteins)
			if (!p.input && p.output) allProteinsPtr.push_back({&p, true});
		for (auto& g : subNets) {
			for (auto& p : g.actualProteins) {
				// subs inputs are outputs for us
				if (p.input && !p.output) allProteinsPtr.push_back({&p, false});
			}
		}
		// for (auto& sn : subNets) sn.updateAllProteinsPtr();
	}

	signature_t getSignatures() { return signatures; }

	void step(unsigned int nbSteps = 1) {
		for (unsigned int i = 0; i < nbSteps; ++i) {
			Implem::step(*this);
			for (auto& g : subNets) g.step();
		}
	}

	inline void reset() {
		for (auto& p : actualProteins) p.reset();
		for (auto& g : subNets) g.reset();
	}

	void setParam(size_t i, double val) {
		if (i < params.size()) params[i] = val;
	}

	void randomParams() {
		auto limits = Implem::paramsLimits();
		for (size_t i = 0; i < Implem::nbParams; ++i) {
			std::uniform_real_distribution<double> distrib(limits[i].first, limits[i].second);
			params[i] = distrib(grnRand);
		}
		master->updateSubNetsPtrsAndSignatures();
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

	template <typename K, typename T> bool inMap(std::map<K, T> m, T elem) {
		for (auto& e : m)
			if (e.second == elem) return true;
		return false;
	}

	void enforceNamedInputsAndOutputs() {
		for (size_t i = 0; i < actualProteins.size(); ++i) {
			auto& p = actualProteins[i];
			if (!inMap(inputProteins, i)) p.input = false;
			if (!inMap(outputProteins, i)) p.output = false;
		}
		for (auto ip : inputProteins) {
			actualProteins[ip.second].modifiable = false;
		}
		for (auto op : outputProteins) {
			actualProteins[op.second].modifiable = false;
		}
	}

	void enforceOneInputOneOutput() {
		if (actualProteins.size() == 0) addRandomProtein();
		assert(actualProteins.size() > 0);
		bool inp = false;
		bool out = false;
		for (auto& p : actualProteins) {
			if (p.input) inp = true;
			if (p.output) out = true;
		}
		std::uniform_int_distribution<size_t> dice(0, actualProteins.size() - 1);
		if (!inp) actualProteins[dice(grnRand)].input = true;
		if (!out) actualProteins[dice(grnRand)].output = true;
		master->updateSubNetsPtrsAndSignatures();
	}

	std::vector<MGRN<Implem>*> getListOfAllGRNs() const {
		// returns all GRNS present in the whole network
		std::unordered_set<MGRN<Implem>*> visited;
		std::vector<MGRN<Implem>*> toVisit;
		toVisit.push_back(master);
		while (toVisit.size()) {
			auto g = toVisit.back();
			toVisit.pop_back();
			visited.insert(g);
			for (auto& sg : g->subNets)
				if (!visited.count(&sg)) toVisit.push_back(&sg);
		}
		std::vector<MGRN<Implem>*> result;
		result.reserve(visited.size());
		for (auto g : visited) result.push_back(g);
		return result;
	}

	std::vector<std::pair<MGRN<Implem>*, Protein*>> getListOfAllProteins() const {
		// returns all of the proteins present in the whole network
		std::vector<std::pair<MGRN<Implem>*, Protein*>> result;
		auto allGrns = getListOfAllGRNs();
		for (auto g : allGrns) {
			for (size_t i = 0; i < g->actualProteins.size(); ++i) {
				result.push_back({g, &g->actualProteins[i]});
			}
		}
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
		subNets.erase(subNets.begin() + static_cast<int>(i));
		master->updateSubNetsPtrsAndSignatures();
	}

	void deleteRandomProtein() {
		if (actualProteins.size()) {
			std::uniform_int_distribution<size_t> dp(0, actualProteins.size() - 1);
			deleteProtein(dp(grnRand));
		}
	}

	void deleteRandomRegul() {
		if (actualProteins.size()) {
			std::vector<size_t> regulList;
			for (size_t i = 0; i < actualProteins.size(); ++i) {
				auto& p = actualProteins[i];
				if ((p.input && p.output) || (!p.input && !p.output)) regulList.push_back(i);
			}
			if (regulList.size() > 0) {
				std::uniform_int_distribution<size_t> dp(0, regulList.size() - 1);
				deleteProtein(regulList[dp(grnRand)]);
			}
		}
	}

	void deleteProtein(size_t i) {
		assert(i < actualProteins.size());
		actualProteins.erase(actualProteins.begin() + static_cast<int>(i));
		if (parent) enforceOneInputOneOutput();
		master->updateSubNetsPtrsAndSignatures();
	}

	/**************************************
	     *              JSON
	*************************************/

	string serialize() const { return toJSON().dump(2); }

	json toJSON() const {
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
		return o;
	}

	bool isMaster() { return master == this; }

	/**************************************
	 *       MUTATION & CROSSOVER
	 *************************************/

	size_t mutate() {
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
				if (allP[id].first->isMaster()) allP[id].first->enforceNamedInputsAndOutputs();
			} break;
			case 1: {
				// add prot
				auto allG = getListOfAllGRNs();
				std::uniform_int_distribution<size_t> dice(0, allG.size() - 1);
				auto g = allG[dice(grnRand)];
				g->addRandomProtein();
			} break;
			case 2: {
				// del prot
				auto allG = getListOfAllGRNs();
				std::uniform_int_distribution<size_t> diceG(0, allG.size() - 1);
				auto g = allG[diceG(grnRand)];
				if (g->actualProteins.size() > 1) {
					std::uniform_int_distribution<size_t> diceP(0, g->actualProteins.size() - 1);
					auto pid = diceP(grnRand);
					if (g->actualProteins[pid].modifiable) {
						g->deleteProtein(pid);
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
				g->addSubNet(ng);
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
		return operation;
	}

	bool operator!=(const MGRN<Implem>& other) const { return !(*this == other); }
	bool operator==(const MGRN<Implem>& other) const {
		// checking that they have the same params
		for (size_t i = 0; i < params.size(); ++i)
			if (other.params[i] != params[i]) return false;
		// proteins and subNets orders are not importants...
		// this forces us to do more computations
		if (actualProteins.size() != other.actualProteins.size()) return false;
		if (subNets.size() != other.subNets.size()) return false;
		for (auto& p0 : actualProteins) {
			bool found = false;
			for (auto& p1 : other.actualProteins) {
				if (p1 == p0) {
					found = true;
					break;
				}
			}
			if (!found) {
				return false;
			}
		}
		for (auto& p0 : other.actualProteins) {
			bool found = false;
			for (auto& p1 : actualProteins) {
				if (p1 == p0) {
					found = true;
					break;
				}
			}
			if (!found) return false;
		}
		for (auto& s0 : subNets) {
			bool found = false;
			for (auto& s1 : other.subNets) {
				if (s1 == s0) {
					found = true;
					break;
				}
			}
			if (!found) return false;
		}
		for (auto& s0 : other.subNets) {
			bool found = false;
			for (auto& s1 : subNets) {
				if (s1 == s0) {
					found = true;
					break;
				}
			}
			if (!found) return false;
		}
		return true;
	}

	MGRN<Implem> crossover(const MGRN<Implem>& other) {
		return MGRN::crossover(*this, other);
	}

	static MGRN<Implem> crossover(const MGRN<Implem>& g0, const MGRN<Implem>& g1) {
		std::uniform_int_distribution<int> d(0, 1);
		return d(grnRand) == 1 ? g0 : g1;
	};
};
#endif
