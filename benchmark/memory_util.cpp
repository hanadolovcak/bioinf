// Author: Hana Dolovčak
#include "memory_util.h"

#include <sys/resource.h>

#include <cstdio>
#include <cstring>
#include <fstream>
#include <string>

#if defined(__APPLE__)
#include <mach/mach.h>
#endif

namespace bamboo {

// Peak resident-set memory the process has used so far, in kilobytes.
size_t peak_rss_kb() {
    struct rusage ru;
    if (getrusage(RUSAGE_SELF, &ru) != 0) return 0;
#if defined(__APPLE__)
    // macOS reports ru_maxrss in bytes; Linux reports kilobytes.
    return static_cast<size_t>(ru.ru_maxrss) / 1024u;
#else
    return static_cast<size_t>(ru.ru_maxrss);
#endif
}

// Current resident-set memory of the process, in kilobytes.
size_t current_rss_kb() {
#if defined(__APPLE__)
    mach_task_basic_info info;
    mach_msg_type_number_t count = MACH_TASK_BASIC_INFO_COUNT;
    if (task_info(mach_task_self(), MACH_TASK_BASIC_INFO,
                  reinterpret_cast<task_info_t>(&info), &count) != KERN_SUCCESS) {
        return 0;
    }
    return static_cast<size_t>(info.resident_size) / 1024u;
#else
    std::ifstream f("/proc/self/status");
    if (!f) return 0;
    std::string line;
    while (std::getline(f, line)) {
        if (line.rfind("VmRSS:", 0) == 0) {
            size_t kb = 0;
            std::sscanf(line.c_str(), "VmRSS: %zu kB", &kb);
            return kb;
        }
    }
    return 0;
#endif
}

}
