// Author: Hana Dolovčak
#include <cassert>
#include <cstdio>
#include <fstream>
#include <string>

#include "fasta_parser.h"

int main() {
    using namespace bamboo;

    const char* path = "/tmp/bamboo_test_fasta.fa";
    {
        std::ofstream f(path);
        f << ">seq1 some header\n";
        f << "ACGTacgt\n";
        f << "NNNAA\n";
        f << ">seq2\n";
        f << "  GG\tCC\n";
    }
    {
        std::string s = load_fasta(path, /*strict_acgt=*/true);
        assert(s == "ACGTACGTAAGGCC");
    }
    {
        std::ofstream f(path);
        f << ">x\n";
        f << "ACG\n";
    }
    {
        std::string s = load_fasta(path);
        assert(s == "ACG");
    }
    std::remove(path);
    std::puts("test_fasta OK");
    return 0;
}
