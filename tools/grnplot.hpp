#ifndef GRNPLOT_HPP
#define GRNPLOT_HPP
#include <stdio.h>
#include <map>
#include <vector>
#include <string>
#include <functional>
#include <utility>
#include "../common.h"

template <typename GRN>
void plot(const std::function<void(GRN &g)> &updtFunc, GRN &grn, const int nbSteps,
          const std::string &name) {
	std::vector<std::string> outNames = grn.getProteinNames(ProteinType::output);
	std::map<std::string, std::vector<std::pair<double, double>>> outputs;
	for (double i = 0.0; i < static_cast<double>(nbSteps); ++i) {
		for (auto &pname : outNames) {
			outputs[pname].push_back(std::make_pair(i, grn.getProteinConcentration(pname, ProteinType::output)));
		}
		updtFunc(grn);
	}
	FILE *gnuplotPipe = popen("gnuplot -persistent", "w");
	fprintf(
	    gnuplotPipe,
	    "set terminal pngcairo transparent enhanced font \"arial,10\" fontscale 1.0 size "
	    "1200,800\nset output '%s.png'\nset yrange [*<-0.05:1.05<*]\nset ylabel "
	    "\"concenntration\"\nset xlabel \"step\"\nset size ratio 0.5\nplot ",
	    name.c_str());
	size_t i = 0;
	for (auto &o : outputs) {
		fprintf(gnuplotPipe, "'-' using 1:2 title \"%s\" with lines", o.first.c_str());
		if (i++ == outputs.size() - 1)
			fprintf(gnuplotPipe, "\n");
		else
			fprintf(gnuplotPipe, ", \\\n");
	}
	for (auto &o : outputs) {
		for (auto &p : o.second) {
			fprintf(gnuplotPipe, "%lf %lf\n", p.first, p.second);
		}
		fprintf(gnuplotPipe, "EOF");
	}
	fprintf(gnuplotPipe, "e");
}
#endif
