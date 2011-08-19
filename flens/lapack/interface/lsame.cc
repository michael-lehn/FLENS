#include <cassert>
#include <cctype>


#define LSAME   lsame_

extern "C" {

int
LSAME(const char *c1, const char *c2)
{
    return (toupper(*c1)==toupper(*c2)) ? 1 : 0;
}

} // extern "C"