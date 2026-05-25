// Author: Hana Dolovčak
#ifndef BAMBOO_MEMORY_UTIL_H
#define BAMBOO_MEMORY_UTIL_H

#include <cstddef>

namespace bamboo {

// Peak resident-set memory the process has used so far, in kilobytes.
size_t peak_rss_kb();

// Current resident-set memory of the process, in kilobytes.
size_t current_rss_kb();

}

#endif
