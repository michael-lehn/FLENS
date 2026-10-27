#include <flens/lapack/interface/include/config.h>

extern "C" {

//-- lsame ---------------------------------------------------------------------
LOGICAL
LAPACK_DECL(lsame)(const char       *CA,
                   const char       *CB)
{
//    DEBUG_FLENS_LAPACK("lsame");
//    std::cerr << "*CA = " << *CA << ", *CB = " << *CB << std::endl;
//
//  Compare regardless of case (ASCII) like the reference implementation
//
    char a = *CA;
    char b = *CB;
    if (a>='a' && a<='z') {
        a = a - 'a' + 'A';
    }
    if (b>='a' && b<='z') {
        b = b - 'a' + 'A';
    }
    return (a==b) ? LOGICAL(1) : LOGICAL(0);
}

} // extern "C"
