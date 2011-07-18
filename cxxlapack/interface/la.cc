#include <cxxlapack/cxxlapack.cxx>
#include <cxxblas/cxxblas.cxx>
#include <cxxlapack/interface/aux.h>

#include <cmath>


#ifdef SINGLE
#   define LABAD        slabad_
#   define LACPY        slacpy_
#   define LANGE        slange_
#   define LANTR        slantr_
#   define LARNV        slarnv_
#   define LARTG        slartg_
#   define LASET        slaset_
#   define LASWP        slaswp_
#   define LAQGE        slaqge_
#endif


#ifdef SINGLE
#   define LABAD_REF    slabad
#   define LACPY_REF    slacpy
#   define LANGE_REF    slange
#   define LANTR_REF    slantr
#   define LARNV_REF    slarnv
#   define LARTG_REF    slartg
#   define LASET_REF    slaset
#   define LASWP_REF    slaswp
#   define LAQGE_REF    slaqge
#endif


using namespace cxxblas;
using namespace cxxlapack;

using std::log10;
using std::sqrt;

extern "C" {

//- labad ----------------------------------------------------------------------

#ifdef USE_LABAD_REF
void
LABAD_REF(FLOAT *, FLOAT *);
#endif

void
LABAD(FLOAT *SMALL, FLOAT *LARGE)
{
#   ifdef DEBUG_INTERFACE
    fprintf(stderr, "LABAD\n");
#   endif


#   ifdef USE_LABAD_REF
    FLOAT SMALL_ = *SMALL;
    FLOAT LARGE_ = *LARGE;
    LABAD_REF(&SMALL_, &LARGE_);
#   else
    if (log10(*LARGE)>2000) {
        *SMALL = sqrt(*SMALL);
        *LARGE = sqrt(*LARGE);
    }
#   endif


#   ifdef USE_LABAD_REF
    if ((*SMALL!=SMALL_) || (*LARGE!=LARGE_)) {
        fprintf(stderr, "SMALL  = %f, LARGE  = %f\n", *SMALL, *LARGE);
        fprintf(stderr, "SMALL_ = %f, LARGE_ = %f\n", SMALL_, LARGE_);
    }
#   endif
}

//- lacpy ----------------------------------------------------------------------

void
LACPY(const char *UPLO, const INT *M, const INT *N,
      const FLOAT *A, const INT *LDA,
      FLOAT *B, const INT *LDB)
{
#   ifdef DEBUG_INTERFACE
    fprintf(stderr, "LACPY\n");
#   endif

#   ifdef USE_LACPY_REF
#       ifdef DEBUG_INTERFACE
        fprintf(stderr, "WARNING: USE_LACPY_REF not supported\n");
#       endif
#   endif

    if ((*UPLO=='u') || (*UPLO=='U')) {
        fprintf(stderr, "LACPY: (U) not implemented\n");
        return;
    }
    if ((*UPLO=='l') || (*UPLO=='L')) {
        fprintf(stderr, "LACPY: (L) not implemented\n");
        return;
    }
    cxxblas::gecopy(ColMajor, NoTrans, *M, *N, A, *LDA, B, *LDB);
}

//- lange ----------------------------------------------------------------------

#ifdef USE_LANGE_REF
FLOAT
LANGE_REF(const char *NORM, const INT *M, const INT *N,
          const FLOAT *A, const INT *LDA, FLOAT *WORK);
#endif


FLOAT
LANGE(const char *NORM, const INT *M, const INT *N,
      const FLOAT *A, const INT *LDA, FLOAT *WORK)
{
#   ifdef DEBUG_INTERFACE
    fprintf(stderr, "LANGE\n");
#   endif

#   ifdef USE_LANGE_REF
    FLOAT value_ = LANGE_REF(NORM, M, N, A, LDA, WORK);
#   endif

    Norm norm;

    if ((*NORM=='M') || (*NORM=='m')) {
        norm = MaximumNorm;
    } else if ((*NORM=='1') || (*NORM=='O') || (*NORM=='o')) {
        norm = OneNorm;
    } else if ((*NORM=='I') || (*NORM=='i')) {
        norm = InfinityNorm;
    } else if ((*NORM=='F') || (*NORM=='f') || (*NORM=='E') || (*NORM=='e')) {
        norm = FrobeniusNorm;
    } else {
        fprintf(stderr, "LANGE: illegal value NORM=%c\n", *NORM);
        assert(0);
    }
    FLOAT value;
    lange(ColMajor, norm, *M, *N, A, *LDA, WORK, value);

#   ifdef USE_LANGE_REF
    if (value_!=value) {
        fprintf(stderr, "LANGE:     value = %f\n", value);
        fprintf(stderr, "LANGE_REF: value = %f\n", value_);
        fprintf(stderr, "diff             = %f\n", value-value_);
        fprintf(stderr, "NORM=%c, M=%d, N=%d, LDA=%d", *NORM, *M, *N, *LDA);
        assert(0);
    }
#   endif

    return value;
}

//- lantr ----------------------------------------------------------------------

#ifdef USE_LANTR_REF
FLOAT
LANTR_REF(const char *NORM, const char *UPLO, const char *DIAG,
          const int *M, const int *N, FLOAT *A, const int *LDA,
          FLOAT *WORK);
#endif

FLOAT
LANTR(const char *NORM, const char *UPLO, const char *DIAG,
      const int *M, const int *N, FLOAT *A, const int *LDA,
      FLOAT *WORK)
{
#   ifdef DEBUG_INTERFACE
    fprintf(stderr, "LANTR\n");
#   endif

#   ifdef USE_LANTR_REF
    FLOAT value_ = LANTR_REF(NORM, UPLO, DIAG, M, N, A, LDA, WORK);
#   endif

    Norm            normA;
    StorageUpLo     upLoA;
    Diag            diagA;

    if ((*NORM=='M') || (*NORM=='m')) {
        normA = MaximumNorm;
    } else if ((*NORM=='1') || (*NORM=='O') || (*NORM=='o')) {
        normA = OneNorm;
    } else if ((*NORM=='I') || (*NORM=='i')) {
        normA = InfinityNorm;
    } else if ((*NORM=='F') || (*NORM=='f') || (*NORM=='E') || (*NORM=='e')) {
        normA = FrobeniusNorm;
    } else {
        fprintf(stderr, "LANTR: illegal value NORM=%c\n", *NORM);
        assert(0);
    }

    if ((*UPLO=='U') || (*UPLO=='u')) {
        upLoA = Upper;
    } else if ((*UPLO=='L') || (*UPLO=='l')) {
        upLoA = Lower;
    } else {
        fprintf(stderr, "LANTR: illegal value UPLO=%c\n", *UPLO);
        assert(0);
    }

    if ((*DIAG=='U') || (*DIAG=='u')) {
        diagA = Unit;
    } else if ((*DIAG=='N') || (*DIAG=='n')) {
        diagA = NonUnit;
    } else {
        fprintf(stderr, "LANTR: illegal value DIAG=%c\n", *DIAG);
        assert(0);
    }

    FLOAT value;
    lantr(ColMajor, normA, upLoA, diagA, *M, *N, A, *LDA, WORK, value);

#   ifdef USE_LANTR_REF
    if (value_!=value) {
        fprintf(stderr, "LANTR:     value = %f\n", value);
        fprintf(stderr, "LANTR_REF: value = %f\n", value_);
        fprintf(stderr, "diff             = %f\n", value-value_);
        fprintf(stderr, "NORM=%c, UPLO=%c, DIAG=%c, M=%d, N=%d, LDA=%d",
                *NORM, *UPLO, *DIAG, *M, *N, *LDA);
    }
#   endif

    return value;
}

//- larnv ----------------------------------------------------------------------

#ifdef USE_LARNV_REF
void
LARNV_REF(const INT *IDIST, INT *ISEED, const INT *N, FLOAT *X);
#endif

void
LARNV(const INT *IDIST, INT *ISEED, const INT *N, FLOAT *X)
{
#   ifdef DEBUG_INTERFACE
    fprintf(stderr, "LARNV\n");
#   endif

    ProbabilityDistribution probDist;
    switch (*IDIST) {
        case 1: probDist = Uniform01;
                break;
        case 2: probDist = Uniform_11;
                break;
        case 3: probDist = StandardNormal;
                break;
        default:
            fprintf(stderr, "LARNV: illegal value IDIST=%d\n", *IDIST);
            assert(0);
    }
    // larnv(probDist, ISEED, *N, X);
#   ifdef USE_LARNV_REF
    LARNV_REF(IDIST, ISEED, N, X);
#   endif
}

//- laset ----------------------------------------------------------------------

#ifdef USE_LASET_REF
void
LASET_REF(const char *UPLO, const INT *M, const INT *N,
          const FLOAT *ALPHA, const FLOAT *BETA,
          FLOAT *A, const INT *LDA);
#endif

void
LASET(const char *UPLO, const INT *M, const INT *N,
      const FLOAT *ALPHA, const FLOAT *BETA,
      FLOAT *A, const INT *LDA)
{
#   ifdef DEBUG_INTERFACE
    fprintf(stderr, "LASET\n");
#   endif

#   ifdef USE_LASET_REF
    LASET_REF(UPLO, M, N, ALPHA, BETA, A, LDA);
#   else
    bool setFull = false;
    StorageUpLo upLoA;
    if ((*UPLO=='u') || (*UPLO=='U')) {
        upLoA = Upper;
    } else if ((*UPLO=='l') || (*UPLO=='L')) {
        upLoA = Lower;
    } else if ((*UPLO=='f') || (*UPLO=='F')) {
        setFull = true;
    } else {
        fprintf(stderr, "LASET: illegal value UPLO=%c\n", *UPLO);
        assert(0);
    }
    if (setFull) {
        laset(ColMajor, *M, *N, *ALPHA, *BETA, A, *LDA);
    } else {
        laset(ColMajor, upLoA, *M, *N, *ALPHA, *BETA, A, *LDA);
    }
#   endif
}

//- laswp ----------------------------------------------------------------------

#ifdef USE_LASWP_REF
void
LASWP_REF(const INT *N, FLOAT *A, const INT *LDA,
          const INT *K1, const INT *K2,
          const INT *IPIV, const INT *INCX);
#endif

void
LASWP(const INT *N, FLOAT *A, const INT *LDA,
      const INT *K1, const INT *K2,
      const INT *IPIV, const INT *INCX)
{
#   ifdef DEBUG_INTERFACE
    fprintf(stderr, "LASWP\n");
#   endif

#   ifdef USE_LASWP_REF
    LASWP_REF(N, A, LDA, K1, K2, IPIV, INCX);
#   else
    laswp(ColMajor, *N, A, *LDA, *K1-1, *K2-1, IPIV, *INCX, 1);
#   endif
}

//- lartg ----------------------------------------------------------------------

#ifdef USE_LARTG_REF
void
LARTG_REF(const FLOAT *F, const FLOAT *G, FLOAT *CS, FLOAT *SN, FLOAT *R);
#endif

void
LARTG(const FLOAT *F, const FLOAT *G, FLOAT *CS, FLOAT *SN, FLOAT *R)
{
#   ifdef DEBUG_INTERFACE
    fprintf(stderr, "LARTG\n");
#   endif

#   ifdef USE_LARTG_REF
    LARTG_REF(F, G, CS, SN, R);
#   else
    lartg(*F, *G, *CS, *SN, *R);
#endif
}

//- laqge ----------------------------------------------------------------------

#ifdef USE_LAQGE_REF
void
LAQGE_REF(const INT *M, const INT *N,
          FLOAT *A, const INT *LDA,
          const FLOAT *R, const FLOAT *C,
          const FLOAT *ROWCND, const FLOAT *COLCND, const FLOAT *AMAX,
          char *EQUED);
#endif

void
LAQGE(const INT *M, const INT *N,
      FLOAT *A, const INT *LDA,
      const FLOAT *R, const FLOAT *C,
      const FLOAT *ROWCND, const FLOAT *COLCND, const FLOAT *AMAX,
      char *EQUED)
{
#   ifdef DEBUG_INTERFACE
    fprintf(stderr, "LAQGE\n");
#   endif
    
#ifdef USE_LAQGE_REF
    LAQGE_REF(M, N,
              A, LDA,
              R, C,
              ROWCND, COLCND, AMAX,
              EQUED);
#else
    Equilibration equilibration;
    laqge(ColMajor, *M, *N, A, *LDA,
          R, C, *ROWCND, *COLCND, *AMAX,
          equilibration);

    switch (equilibration) {
        case NoEquilibration:
            *EQUED = 'N';
            break;
        case RowEquilibration:
            *EQUED = 'R';
            break;
        case ColEquilibration:
            *EQUED = 'C';
            break;
        case RowColEquilibration:
            *EQUED = 'B';
            break;
        default:
            assert(0);
    }
#endif
}

} // extern "C"