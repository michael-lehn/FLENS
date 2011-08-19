#include <cassert>
#include <cctype>

#define LSAMEN  lsamen_

extern "C" {

int
LSAMEN(const int *n, const char *str1, const char *str2)
{
    for (int i=0; i<*n; ++i) {
        if (toupper(str1[i])!=toupper(str2[i])) {
            return 0;
        }
    }
    return 1;
}

} // extern "C"