// Prints SNaQ.jl's tcgidentifiable(N) for each network of a file (one extended
// Newick per line), after semidirectnetwork!(N), as tools/julia/tcgcheck.jl does.
//
// usage: snaqpp_tcgcheck <networks file>
// output: one line per network: "true", "false" or "error: <message>"
#include <fstream>
#include <iostream>
#include <string>

#include "identifiability.hpp"
#include "newick.hpp"
#include "snaq_network.hpp"

using namespace snaqpp;

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "usage: snaqpp_tcgcheck <networks file>\n";
        return 2;
    }
    std::ifstream in(argv[1]);
    if (!in) {
        std::cerr << "cannot open " << argv[1] << "\n";
        return 2;
    }
    std::string line;
    while (std::getline(in, line)) {
        if (line.find_first_not_of(" \t\r\n") == std::string::npos) continue;
        try {
            Network net = readnewick_string(line);
            semidirectnetwork(net);
            std::cout << (tcgidentifiable(net) ? "true" : "false") << "\n";
        } catch (const std::exception& e) {
            std::cout << "error: " << e.what() << "\n";
        }
    }
    return 0;
}
