#define ILAVER  ilaver_

extern "C" {

void
ILAVER(int *major, int *minor, int *patch)
{
    *major = 3;
    *minor = 3;
    *patch = 1;
}

} // extern "C"