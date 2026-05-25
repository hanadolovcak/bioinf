// Author: Hana Dolovčak
#include "fasta_parser.h"

#include <cctype>
#include <fstream>
#include <stdexcept>

namespace bamboo {

namespace {

bool is_acgt(char c) {
    return c == 'A' || c == 'C' || c == 'G' || c == 'T';
}

}

// Loads a FASTA file and returns its concatenated sequence (header lines
// dropped). When strict_acgt is true (default), only A/C/G/T characters
// are kept; otherwise whitespace is the only thing dropped.
std::string load_fasta(const std::string& path, bool strict_acgt) {
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("load_fasta: cannot open " + path);
    }
    std::string line;
    std::string out;
    out.reserve(1u << 20);
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') continue;
        for (char c : line) {
            char up = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
            if (strict_acgt) {
                if (is_acgt(up)) out.push_back(up);
            } else {
                if (!std::isspace(static_cast<unsigned char>(up))) {
                    out.push_back(up);
                }
            }
        }
    }
    return out;
}

}
