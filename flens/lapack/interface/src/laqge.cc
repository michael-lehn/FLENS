#define STR(x)      #x
#define STRING(x)   STR(x)

#include <flens/lapack/interface/include/config.h>


namespace flens { namespace lapack {

extern "C" {

//-- dlaqge --------------------------------------------------------------------
void
LAPACK_DECL(dlaqge)(const INTEGER    *M,
                    const INTEGER    *N,
                    DOUBLE           *A,
                    const INTEGER    *LDA,
                    const DOUBLE     *R,
                    const DOUBLE     *C,
                    const DOUBLE     *ROWCND,
                    const DOUBLE     *COLCND,
                    const DOUBLE     *AMAX,
                    char             *EQUED)
{
//
//  Call FLENS implementation
//
    DGeMatrixView          A_  = DFSView(*M, *N, A, *LDA);
    DConstDenseVectorView  R_  = DConstArrayView(*M, R, 1);
    DConstDenseVectorView  C_  = DConstArrayView(*N, C, 1);

    const auto equed = laq(A_, R_, C_, *ROWCND, *COLCND, *AMAX);
    *EQUED = char(equed);
}

} // extern "C"

} } // namespace lapack, flens
