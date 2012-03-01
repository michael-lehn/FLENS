#include <flens/lapack/interface/include/config.h>


namespace flens { namespace lapack {

#if defined CHECK_CXXLAPACK

extern "C" {

//-- ilaenv --------------------------------------------------------------------
INTEGER
LAPACK_DECL(ilaenv)(const INTEGER   *SPEC,
                    const char      *NAME,
                    const char      *OPTS,
                    const INTEGER   *N1,
                    const INTEGER   *N2,
                    const INTEGER   *N3,
                    const INTEGER   *N4,
                    int             ,
                    int             )
{
    //std::cerr << "LAPACK_DECL(ilaenv) called!" << std::endl;
    return ilaenv_generic<double>(*SPEC, NAME, OPTS, *N1, *N2, *N3, *N4);
}

} // extern "C"

#endif

} } // namespace lapack, flens
