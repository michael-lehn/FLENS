#include <flens/blas/interface/blas/config.h>


using namespace flens;

extern "C" {

void
BLAS(snrm2_sub)(const INTEGER   *N,
                const float     *X,
                const INTEGER   *INCX,
                float           *RES)
{
#   ifdef TEST_DIRECT_CBLAS

        *RES = cblas_snrm2(*N, X, *INCX);

#   else

        using std::abs;

        SDenseVectorConstView  x(SConstArrayView(*N, X, abs(*INCX)), *INCX<0);

        *RES = blas::nrm2(x);

#   endif

}

void
BLAS(dnrm2_sub)(const INTEGER   *N,
                const double    *X,
                const INTEGER   *INCX,
                double          *RES)
{
#   ifdef TEST_DIRECT_CBLAS

    *RES = cblas_dnrm2(*N, X, *INCX);

#   else

        using std::abs;

        DDenseVectorConstView  x(DConstArrayView(*N, X, abs(*INCX)), *INCX<0);

        *RES = blas::nrm2(x);

#   endif
}

void
BLAS(scnrm2_sub)(const INTEGER   *N,
                 const cfloat    *X,
                 const INTEGER   *INCX,
                 float           *RES)
{
#   ifdef TEST_DIRECT_CBLAS

        *RES = cblas_scnrm2(*N, reinterpret_cast<const float *>(X), *INCX);

#   else

        using std::abs;

        CDenseVectorConstView  x(CConstArrayView(*N, X, abs(*INCX)), *INCX<0);

        *RES = blas::nrm2(x);

#   endif
}

void
BLAS(dznrm2_sub)(const INTEGER   *N,
                 const cdouble   *X,
                 const INTEGER   *INCX,
                 double          *RES)
{
#   ifdef TEST_DIRECT_CBLAS

        *RES = cblas_dznrm2(*N, reinterpret_cast<const double *>(X), *INCX);

#   else

        using std::abs;

        ZDenseVectorConstView  x(ZConstArrayView(*N, X, abs(*INCX)), *INCX<0);

        *RES = blas::nrm2(x);
#   endif
}

} // extern "C"
