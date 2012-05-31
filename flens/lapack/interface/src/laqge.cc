#define STR(x)      #x
#define STRING(x)   STR(x)

#define FLENS_DEFAULT_INDEXTYPE int

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
    DGeMatrixView          _A  = DFSView(*M, *N, A, *LDA);
    DConstDenseVectorView  _R  = DConstArrayView(*M, R, 1);
    DConstDenseVectorView  _C  = DConstArrayView(*N, C, 1);

    const auto equed = laq(_A, _R, _C, *ROWCND, *COLCND, *AMAX);
    *EQUED = char(equed);
}

} // extern "C"

} } // namespace lapack, flens
