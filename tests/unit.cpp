#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "../common.h"
#include "../grn.hpp"
#include "../classic.hpp"
#include "../json/json.hpp"

TEST_CASE("Proteins are ok", "[proteins]") {
	Protein<1> p0({{21.1}}, 0.5);
	Protein<1> p1({{0.1234567890}}, 0.7);
	Protein<2, int, 100, 110> p2({{12, 107}}, 2.1234567891234567891234567891234);
	REQUIRE(decltype(p0)::getMaxDistance() == 1.0);
	REQUIRE(decltype(p2)::getMaxDistance() == sqrt(200.0));

	SECTION("Proteins are correctly explicitly constructed") {
		REQUIRE(p0.coords.size() == 1);
		REQUIRE(p1.coords.size() == 1);
		REQUIRE(p2.coords.size() == 2);
		REQUIRE(p0.coords[0] == 1.0);
		REQUIRE(p0.c == 0.5);
		REQUIRE(p1.coords[0] == 0.1234567890);
		REQUIRE(p1.c == 0.7);
		REQUIRE(p2.coords[0] == 100);
		REQUIRE(p2.coords[1] == 107);
		REQUIRE(p2.c == 2.1234567891234567891234567891234);
	}

	SECTION("Proteins are correctly randomly constructed") {
		Protein<3> pr;
		REQUIRE(decltype(pr)::getMaxDistance() == sqrt(3.0));
		Protein<2, int, 100, 110> prbase;
		bool atLeast1Different0 = false;
		bool atLeast1Different1 = false;
		for (int i = 0; i < 500; ++i) {
			Protein<2, int, 100, 110> pr2;
			REQUIRE(pr2.c == pr.c);
			REQUIRE(prbase.c == pr.c);
			REQUIRE(pr2.coords[0] <= 110);
			REQUIRE(pr2.coords[0] >= 100);
			REQUIRE(pr2.coords[1] <= 110);
			REQUIRE(pr2.coords[1] >= 100);
			if (prbase.coords[0] != pr2.coords[0]) atLeast1Different0 = true;
			if (prbase.coords[1] != pr2.coords[1]) atLeast1Different1 = true;
		}
		REQUIRE(atLeast1Different0);
		REQUIRE(atLeast1Different1);
	}

	SECTION("Proteins are correctly copy constructed") {
		auto pc0 = p0;
		auto pc1 = p1;
		auto pc2 = p2;
		decltype(pc2) pcc2(pc2);
		REQUIRE(!decltype(pc0)::areSame(pc0, pc1));
		REQUIRE(decltype(pc0)::areSame(pc0, p0));
		REQUIRE(decltype(pc0)::areSame(pc1, p1));
		REQUIRE(decltype(pc2)::areSame(pc2, p2));
		REQUIRE(decltype(pc2)::areSame(pcc2, pc2));
		REQUIRE(pc0.coords[0] == 1.0);
		REQUIRE(pc0.c == 0.5);
		REQUIRE(pc1.coords[0] == 0.1234567890);
		REQUIRE(pc1.c == 0.7);
		REQUIRE(pc2.coords[0] == 100);
		REQUIRE(pc2.coords[1] == 107);
		REQUIRE(pc2.c == 2.1234567891234567891234567891234);
		REQUIRE(pcc2.coords[0] == 100);
		REQUIRE(pcc2.coords[1] == 107);
		REQUIRE(pcc2.c == 2.1234567891234567891234567891234);
	}

	SECTION("Proteins are correctly saved & read") {
		auto js0 = p0.toJSON().dump();
		auto pc0 = decltype(p0)(nlohmann::json::parse(js0));
		REQUIRE(decltype(pc0)::areSame(pc0, p0));
		REQUIRE(pc0.c == p0.c);
		REQUIRE(pc0.coords[0] == p0.coords[0]);

		p2.c = 12435894651.4646545812345678901;
		auto js2 = p2.toJSON();
		auto pc2 = decltype(p2)(js2);
		REQUIRE(decltype(pc2)::areSame(pc2, p2));
		REQUIRE(pc2.c == p2.c);
		REQUIRE(pc2.c == 12435894651.4646545812345678901);

		auto js3 = p2.toJSON().dump(2);
		auto pc3 = decltype(p2)(nlohmann::json::parse(js3));
		REQUIRE(decltype(pc3)::areSame(pc3, p2));
		REQUIRE(pc3.c == p2.c);
		REQUIRE(pc3.c == 12435894651.4646545812345678901);
	}

	SECTION("distance") {
		REQUIRE(p0.getDistanceWith(p0) == 0.0);
		REQUIRE(p1.getDistanceWith(p0) > 0.0);
	}
}
TEST_CASE("GRNs are ok", "[grn]") {
	GRN<Classic> grn;
	auto params = grn.getParams();
	REQUIRE(params.size() == 2);
	REQUIRE(params[0] == 0);
	REQUIRE(params[1] == 0);
	using Prot = typename Classic::Protein_t;
	Prot p0({{1, 2, 1}}, 0.5);
	Prot p1({{5, 15, 3}}, 0.5);
	Prot p2({{10, 9, 14}}, 0.5);
	Prot p3({{8, 12, 7}}, 0.5);
	grn.addProtein(ProteinType::input, "i0", p0);
	grn.addProtein(ProteinType::input, "i1", p1);
	grn.addProtein(ProteinType::regul, "r0", p2);
	grn.addProtein(ProteinType::output, "o0", p3);
	REQUIRE(grn.getNbProteins() == 4);
	REQUIRE(grn.getProteinSize(ProteinType::input) == 2);
	REQUIRE(grn.getProteinSize(ProteinType::regul) == 1);
	REQUIRE(grn.getProteinSize(ProteinType::output) == 1);
	REQUIRE(Prot::areSame(grn.getProtein(ProteinType::input, "i0"), p0));
	REQUIRE(Prot::areSame(grn.getProtein(ProteinType::input, "i1"), p1));
	REQUIRE(Prot::areSame(grn.getProtein(ProteinType::regul, "r0"), p2));
	REQUIRE(Prot::areSame(grn.getProtein(ProteinType::output, "o0"), p3));
	auto sig = grn.getSignatures();
	REQUIRE(sig.size() == grn.getNbProteins());
	REQUIRE(sig[0].size() == grn.getNbProteins());

	/////////////////////////////////////

	grn.step(100);
	params = grn.getParams();
	REQUIRE(params.size() == 2);
	REQUIRE(params[0] == 0);
	REQUIRE(params[1] == 0);
	REQUIRE(sig.size() == grn.getNbProteins());
	REQUIRE(sig[0].size() == grn.getNbProteins());
	REQUIRE(grn.getNbProteins() == 4);
	REQUIRE(grn.getProteinSize(ProteinType::input) == 2);
	REQUIRE(grn.getProteinSize(ProteinType::regul) == 1);
	REQUIRE(grn.getProteinSize(ProteinType::output) == 1);
	REQUIRE(Prot::areSame(grn.getProtein(ProteinType::input, "i0"), p0));
	REQUIRE(Prot::areSame(grn.getProtein(ProteinType::input, "i1"), p1));
	REQUIRE(Prot::areSame(grn.getProtein(ProteinType::regul, "r0"), p2));
	REQUIRE(Prot::areSame(grn.getProtein(ProteinType::output, "o0"), p3));
	REQUIRE(grn.getProtein(ProteinType::input, "i0").c == 0.5);
	REQUIRE(grn.getProtein(ProteinType::input, "i1").c == 0.5);
	REQUIRE(grn.getProtein(ProteinType::regul, "r0").c == 0.5);
	REQUIRE(grn.getProtein(ProteinType::output, "o0").c == 0.5);
}
