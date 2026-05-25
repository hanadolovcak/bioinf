// Author: Hana Dolovčak
#ifndef BAMBOO_FASTA_PARSER_H
#define BAMBOO_FASTA_PARSER_H

#include <string>

namespace bamboo {

// Loads a FASTA file and returns its concatenated sequence (header lines
// dropped). When strict_acgt is true (default), only A/C/G/T characters are
// kept and everything else is filtered out. The ntHash function we use is
// only defined on the ACGT alphabet, so callers feeding the filter should
// keep this on.
std::string load_fasta(const std::string& path, bool strict_acgt = true);

}

#endif
