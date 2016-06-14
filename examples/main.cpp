#include <thread>
#include "../mgclassic.hpp"
#include "../mgrn.hpp"

int main(int, char**) {
	const int nbSteps = 500;
	MGRN<MGClassic> grn;
	grn.addRandomProtein(ProteinType::input, "Input 0");
	grn.addRandomProtein(ProteinType::input, "Input 1");
	grn.addRandomProtein(ProteinType::output, "Output 0");
	grn.addRandomProtein(ProteinType::output, "Output 1");
	grn.randomReguls(20);
	grn.randomParams();
	for (int i = 0; i < nbSteps; ++i) {
		grn.step();
		std::cout << "{ \"GRN\":" << grn.serialize() << "}" << std::endl;
		std::this_thread::sleep_for(0.5s);
	}
	return 0;
}

