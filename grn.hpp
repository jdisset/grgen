#ifndef GENERICGRN_HPP
#define GENERICGRN_HPP
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

template <typename Implem> class GRN {
	struct GAConfiguration {
		// crossover
		static constexpr double ALIGN_TRESHOLD = 0.4;
		static constexpr double APPEND_NON_ALIGNED = 0.4;
		static constexpr unsigned int MAX_REGULS = 40;

		// mutation
		double MODIF_RATE = 0.5;
		double ADD_RATE = 0.25;
		double DEL_RATE = 0.25;
	};

	friend Implem;

 public:
	using Protein = typename Implem::Protein_t;
	using json = nlohmann::json;
	using InfluenceVec = array<double, Implem::nbSignatureParams>;
	template <typename A, typename B> using umap = std::unordered_map<A, B>;
	GAConfiguration config;

 protected:
	array<double, Implem::nbParams> params;    // alpha, beta, ...
	array<umap<string, Protein>, 3> proteins;  // the actual proteins
	umap<Protein*, umap<Protein*, InfluenceVec>>
	    signatures;  // stores the influence of one protein onto the others (p0 influences
	                 // p1 with value [p0][p1])
	int currentStep = 0;
	Implem implem;

 public:
	GRN() { updateSignatures(); }

	GRN(const GRN& grn)
	    : params(grn.params), proteins(grn.proteins), currentStep(grn.currentStep) {
		updateSignatures();
	}
	/**************************************
	 *          UPDATES
	 *************************************/
	void updateSignatures() { implem.updateSignatures(*this); }

	void step(unsigned int nbSteps = 1) { implem.step(*this, nbSteps); }

	/**************************************
	 *               GET
	 *************************************/
	double getProteinConcentration(const string& name, const ProteinType t) const {
		return proteins.at((size_t)t).at(name).c;
	}
	array<double, Implem::nbParams> getParams() const { return params; }
	array<map<string, Protein>, 3> getProteins() const { return proteins; }
	inline size_t getProteinSize(ProteinType t) const {
		return proteins[to_underlying(t)].size();
	}
	size_t getNbProteins() const {
		size_t n = 0;
		for (const auto& i : proteins) n += i.size();
		return n;
	}
	int getCurrentStep() const { return currentStep; }
	Protein& getProtein(ProteinType t, const string& name) {
		assert(proteins[to_underlying(t)].count(name) > 0);
		return proteins[to_underlying(t)][name];
	}
	/**************************************
	 *               SET
	 *************************************/
	void reset() {
		for (size_t t0 = 0; t0 < proteins.size(); ++t0) {
			for (auto p0 = proteins[t0].begin(); p0 != proteins[t0].end(); ++p0) {
				p0->second.reset();
			}
		}
	}
	void setParam(size_t i, double val) {
		if (i < params.size()) params[i] = val;
	}
	void setProteinConcentration(const string& name, ProteinType t, double c) {
		proteins[(size_t)t][name].c = c;
	}
	vector<string> getProteinNames(ProteinType t) const {
		vector<string> res;
		for (auto& p : proteins[(size_t)t]) {
			res.push_back(p.first);
		}
		return res;
	}

	void randomParams() {
		array<pair<double, double>, Implem::nbParams> limits = Implem::paramsLimits();
		for (size_t i = 0; i < Implem::nbParams; ++i) {
			std::uniform_real_distribution<double> distrib(limits[i].first, limits[i].second);
			params[i] = distrib(grnRand);
		}
		updateSignatures();
	}
	/**************************************
	 *          ADDING PROTEINS
	 *************************************/
	void addProtein(const ProteinType t, const string& name, const Protein& p) {
		proteins[to_underlying(t)].insert(make_pair(name, Protein(p)));
		updateSignatures();
	}
	void addRandomProtein(const ProteinType t, const string& name) {
		proteins[to_underlying(t)].insert(make_pair(name, Protein()));
	}
	void addProteins(umap<string, Protein>& prots, const ProteinType t) {
		for (auto& p : prots) {
			addProtein(t, p.first, p.second);
		}
	}
	void randomReguls(size_t n) {
		ostringstream name;
		proteins[to_underlying(ProteinType::regul)].clear();
		for (size_t i = 0; i < n; ++i) {
			name.str("");
			name.clear();
			name << "r" << i;
			addProtein(ProteinType::regul, name.str(), Protein());
		}
		updateSignatures();
	}
	void updateRegulNames() {
		int id = 0;
		umap<string, Protein> newReguls;
		for (auto& i : proteins[to_underlying(ProteinType::regul)]) {
			ostringstream name;
			name << "r" << id++;
			newReguls[name.str()] = i.second;
		}
		proteins[to_underlying(ProteinType::regul)] = newReguls;
	};

	/**************************************
	 *       MUTATION & CROSSOVER
	 *************************************/
	void mutate() {
		std::uniform_real_distribution<double> dReal(0.0, 1.0);
		// mutate params
		if (dReal(grnRand) < config.MODIF_RATE / static_cast<double>(getNbProteins())) {
			array<pair<double, double>, Implem::nbParams> limits = Implem::paramsLimits();
			std::uniform_int_distribution<int> dInt(0, Implem::nbParams - 1);
			size_t mutParam = dInt(grnRand);
			std::uniform_real_distribution<double> distrib(limits[mutParam].first,
			                                               limits[mutParam].second);
			params[mutParam] = distrib(grnRand);
		}
		// mutate proteins
		else {
			vector<string> reguls = getProteinNames(ProteinType::regul);
			// mutate 1 protein
			double dval = dReal(grnRand);
			if (dval <=
			    config.MODIF_RATE / (config.MODIF_RATE + config.ADD_RATE + config.DEL_RATE)) {
				if (reguls.size() > 0) {
					std::uniform_int_distribution<int> dInt(0, reguls.size() - 1);
					int v = dInt(grnRand);
					proteins[to_underlying(ProteinType::regul)][reguls[v]].mutate();
				}
			}
			// ajout
			else if (dval <= (config.MODIF_RATE + config.ADD_RATE) /
			                     (config.MODIF_RATE + config.ADD_RATE + config.DEL_RATE)) {
				ostringstream name;
				name << "r" << reguls.size();
				addProtein(ProteinType::regul, name.str(), Protein());
			}
			// suppression
			else {
				if (reguls.size() > 0) {
					std::uniform_int_distribution<int> dInt(0, reguls.size() - 1);
					auto it = proteins[to_underlying(ProteinType::regul)].begin();
					advance(it, dInt(grnRand));
					proteins[to_underlying(ProteinType::regul)].erase(it);
					updateRegulNames();
				}
			}
		}
		updateSignatures();
	}

	GRN crossover(const GRN& other) { return GRN::crossover(*this, other); }

	static double computeDistance(const GRN& g0, const GRN& g1) {
		// a rough measure of the distance between 2 grn
		int pclones = 0;
		for (auto t = 0u; t < 3; ++t) {
			for (const auto& pr0 : g0.proteins[t]) {
				for (const auto& pr1 : g1.proteins[t]) {
					if (Protein::areSame(pr0.second, pr1.second)) {
						++pclones;
						break;
					}
				}
			}
		}
		return 1.0 - (static_cast<double>(pclones) * 2.0 /
		              (g0.getNbProteins() + g1.getNbProteins()));
	}

	static GRN crossover(const GRN& g0, const GRN& g1) {
		assert(g0.proteins.size() == g1.proteins.size());
		assert(g0.params.size() == g1.params.size());
		assert(g0.proteins[to_underlying(ProteinType::input)].size() ==
		       g1.proteins[to_underlying(ProteinType::input)].size());
		assert(g0.proteins[to_underlying(ProteinType::output)].size() ==
		       g1.proteins[to_underlying(ProteinType::output)].size());
		GRN offspring;
		assert(offspring.params.size() == g0.params.size());
		std::uniform_int_distribution<int> d5050(0, 1);
		std::uniform_real_distribution<double> dReal(0.0, 1.0);
		// 50/50 for params, inputs and outputs
		for (size_t i = 0; i < g0.params.size(); ++i) {
			offspring.params[i] = d5050(grnRand) ? g0.params[i] : g1.params[i];
		}
		offspring.proteins[to_underlying(ProteinType::input)] =
		    g0.proteins[to_underlying(ProteinType::input)];
		offspring.proteins[to_underlying(ProteinType::output)] =
		    g0.proteins[to_underlying(ProteinType::output)];
		for (auto& i : g1.proteins[to_underlying(ProteinType::input)])
			if (d5050(grnRand))
				offspring.proteins[to_underlying(ProteinType::input)][i.first] = i.second;
		for (auto& i : g1.proteins[to_underlying(ProteinType::output)])
			if (d5050(grnRand))
				offspring.proteins[to_underlying(ProteinType::output)][i.first] = i.second;
		// find closest pairs
		umap<string, Protein> r0 = g0.proteins[to_underlying(ProteinType::regul)];
		umap<string, Protein> r1 = g1.proteins[to_underlying(ProteinType::regul)];
		vector<pair<Protein, Protein>> aligned;  // first = g0's proteins, second = g1's
		double minDist = 0;
		while (minDist < GAConfiguration::ALIGN_TRESHOLD && r0.size() > 0 && r1.size() > 0 &&
		       aligned.size() < GAConfiguration::MAX_REGULS) {
			pair<string, string> closest;
			minDist = std::numeric_limits<double>::infinity();
			for (auto i = r0.begin(); i != r0.end(); ++i) {
				for (auto j = r1.begin(); j != r1.end(); ++j) {
					double dist = i->second.getDistanceWith(j->second);
					if (dist < minDist) {
						closest = make_pair(i->first, j->first);
						minDist = dist;
					}
				}
			}
			if (minDist < GAConfiguration::ALIGN_TRESHOLD) {
				aligned.push_back(
				    pair<Protein, Protein>(r0.at(closest.first), r1.at(closest.second)));
				r0.erase(closest.first);
				r1.erase(closest.second);
			}
		}
		// ProteinType::regul : 50/50 with aligned
		int id = offspring.proteins[to_underlying(ProteinType::regul)].size();
		for (auto& i : aligned) {
			ostringstream name;
			name << "r" << id++;
			if (d5050(grnRand))
				offspring.proteins[to_underlying(ProteinType::regul)][name.str()] = i.first;
			else
				offspring.proteins[to_underlying(ProteinType::regul)][name.str()] = i.second;
		}
		// append the rest (about 1/2 chance)
		for (auto& i : r0) {
			if (offspring.proteins[to_underlying(ProteinType::regul)].size() <
			    GAConfiguration::MAX_REGULS) {
				if (dReal(grnRand) < GAConfiguration::APPEND_NON_ALIGNED) {
					ostringstream name;
					name << "r" << id++;
					offspring.proteins[to_underlying(ProteinType::regul)][name.str()] = i.second;
				}
			}
		}
		for (auto& i : r1) {
			if (offspring.proteins[to_underlying(ProteinType::regul)].size() <
			    GAConfiguration::MAX_REGULS) {
				if (dReal(grnRand) < GAConfiguration::APPEND_NON_ALIGNED) {
					ostringstream name;
					name << "r" << id++;
					offspring.proteins[to_underlying(ProteinType::regul)][name.str()] = i.second;
				}
			}
		}
		offspring.updateSignatures();
		return offspring;
	}

	/**************************************
	 *              JSON
	 *************************************/
	GRN(const string& js) {
		auto o = json::parse(js);
		assert(o.count("params"));
		json par = o.at("params");
		assert(par.size() == Implem::nbParams);
		size_t i = 0;
		for (auto& p : par) {
			double v;
			sscanf(p.get<string>().c_str(), "%lf", &v);
			params[i++] = v;
		}
		assert(o.count("proteins"));
		for (size_t t = to_underlying(ProteinType::input);
		     t <= to_underlying(ProteinType::output); ++t) {
			json prots = o.at("proteins").at(typeToString((ProteinType)t));
			for (json::iterator it = prots.begin(); it != prots.end(); ++it) {
				proteins[t][it.key()] = Protein(it.value());
			}
		}
		updateSignatures();
	}

	string toJSON() const {
		json parArray;
		for (auto& p : params) {
			char buf[50];
			snprintf(buf, sizeof(buf), "%a", p);
			parArray.push_back(buf);
		}
		json protObj;
		for (size_t t = 0; t < proteins.size(); ++t) {
			json pr;
			for (auto p = proteins[t].begin(); p != proteins[t].end(); ++p) {
				pr[p->first] = p->second.toJSON();
			}
			protObj[typeToString((ProteinType)t)] = pr;
		}
		json o;
		o["proteins"] = protObj;
		o["params"] = parArray;
		return o.dump();
	}

	std::string typeToString(ProteinType t) const {
		switch (t) {
			case ProteinType::input:
				return "input";
			case ProteinType::regul:
				return "regul";
			case ProteinType::output:
				return "output";
		}
		return "unknown_type";
	}
};
#endif
