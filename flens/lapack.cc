/*
 *   Copyright (c) 2007, Michael Lehn
 *
 *   All rights reserved.
 *
 *   Redistribution and use in source and binary forms, with or without
 *   modification, are permitted provided that the following conditions
 *   are met:
 *
 *   1) Redistributions of source code must retain the above copyright
 *      notice, this list of conditions and the following disclaimer.
 *   2) Redistributions in binary form must reproduce the above copyright
 *      notice, this list of conditions and the following disclaimer in
 *      the documentation and/or other materials provided with the
 *      distribution.
 *   3) Neither the name of the FLENS development group nor the names of
 *      its contributors may be used to endorse or promote products derived
 *      from this software without specific prior written permission.
 *
 *   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 *   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 *   OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 *   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 *   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 *   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 *   THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 *   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 *   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <cassert>
#include <complex>
#include <flens/lapack.h>

#ifdef VECLIB
#    include <Accelerate/Accelerate.h>
#elif defined MKL
#    ifdef MAC
#        include <Intel_MKL/mkl_lapack.h>
#    else
#        include <mkl_lapack.h>
#    endif
#elif defined GSL_CBLAS
#    include <gsl/gsl_cblas.h>
#else
#    include <cblas.h>
#endif

namespace flens {

extern "C" {

    void
    dpotrf_(char *uplo, int *n, double *ab, int *ldab, int *info);

    void
    sgetrf_(int *m, int *n, float *a, int *lda, int *ipiv, int *info);

    void
    dgetrf_(int *m, int *n, double *a, int *lda, int *ipiv, int *info);

    void
    cgetrf_(int *m, int *n, complex<float> *a, int *lda, int *ipiv, int *info);

    void
    zgetrf_(int *m, int *n, complex<double> *a, int *lda, int *ipiv, int *info);

    void
    sgbtrf_(int *m, int *n, int *kl, int *ku,
            float *ab, int *ldab, int *ipiv, int *info);

    void
    dgbtrf_(int *m, int *n, int *kl, int *ku,
            double *ab, int *ldab, int *ipiv, int *info);

    void
    sgetri_(int *n, float *a, int *lda, const int *ipiv, float *work,
            int *lwork, int *info);

    void
    dgetri_(int *n, double *a, int *lda, const int *ipiv, double *work,
            int *lwork, int *info);

    void
    cgetri_(int *n, complex<float> *a, int *lda, const int *ipiv,
            complex<float> *work, int *lwork, int *info);

    void
    zgetri_(int *n, complex<double> *a, int *lda, const int *ipiv,
            complex<double> *work, int *lwork, int *info);

    void
    sgetrs_(char *trans, int *n, int *nrhs, const float *a, int *lda,
            const int *ipiv, float *b, int *ldb, int *info);

    void
    dgetrs_(char *trans, int *n, int *nrhs, const double *a, int *lda,
            const int *ipiv, double *b, int *ldb, int *info);

    void
    sgbtrs_(char *trans, int *n, int *kl, int *ku, int *nrhs, const float *ab,
            int *ldab, const int *ipiv, float *b, int *ldb, int *info);

    void
    dgbtrs_(char *trans, int *n, int *kl, int *ku, int *nrhs, const double *ab,
            int *ldab, const int *ipiv, double *b, int *ldb, int *info);

    void
    sgesv_(int *n, int *nrhs, float *a, int *lda,
           int *ipiv, float *b, int *ldb, int *info);

    void
    dgesv_(int *n, int *nrhs, double *a, int *lda,
           int *ipiv, double *b, int *ldb, int *info);

    void
    cgesv_(int *n, int *nrhs, complex<float> *a, int *lda,
           int *ipiv, complex<float> *b, int *ldb, int *info);

    void
    zgesv_(int *n, int *nrhs, complex<double> *a, int *lda,
           int *ipiv, complex<double> *b, int *ldb, int *info);

    void
    sgbsv_(int *n, int *kl, int *ku, int *nrhs, float *ab, int *ldab, int *ipiv,
           float *b, int *ldb, int *info);

    void
    dgbsv_(int *n, int *kl, int *ku, int *nrhs, double *ab, int *ldab,
           int *ipiv, double *b, int *ldb, int *info);

    void
    cgbsv_(int *n, int *kl, int *ku, int *nrhs, complex<float> *ab,
           int *ldab, int *ipiv, complex<float> *b, int *ldb, int *info);

    void
    zgbsv_(int *n, int *kl, int *ku, int *nrhs, complex<double> *ab,
           int *ldab, int *ipiv, complex<double> *b, int *ldb, int *info);

    void
    strtrs_(char *uplo, char *trans, char *diag, int *n, int *nrhs,
            const float *a, int *lda, float *b, int *ldb, int *info);

    void
    dtrtrs_(char *uplo, char *trans, char *diag, int *n, int *nrhs,
            const double *a, int *lda, double *b, int *ldb, int *info);

    void
    sgeqrf_(int *m, int *n, float *a, int *lda, float *tau,
            float *work, int *lwork, int *info);

    void
    dgeqrf_(int *m, int *n, double *a, int *lda, double *tau,
            double *work, int *lwork, int *info);

    void
    sorgqr_(int *m, int *n, int *k, float *a, int *lda, const float *tau,
            float *work, int *lwork, int *info);

    void
    dorgqr_(int *m, int *n, int *k, double *a, int *lda, const double *tau,
            double *work, int *lwork, int *info);

    void
    sormqr_(char *side, char *trans, int *m, int *n, int *k,
            const float *a, int *lda, const float *tau, float *c, int *ldc,
            float *work, int *lwork, int *info);

    void
    dormqr_(char *side, char *trans, int *m, int *n, int *k,
            const double *a, int *lda, const double *tau, double *c, int *ldc,
            double *work, int *lwork, int *info);

    void
    sgels_(char *trans, int *m, int *n, int *nrhs, float *a, int *lda,
           float *b, int *ldb, float *work, int *lwork, int *info);

    void
    dgels_(char *trans, int *m, int *n, int *nrhs, double *a, int *lda,
           double *b, int *ldb, double *work, int *lwork, int *info);

    void
    cgels_(char *trans, int *m, int *n, int *nrhs, complex<float> *a,
           int *lda, complex<float> *b, int *ldb,
           complex<float> *work, int *lwork, int *info);

    void
    zgels_(char *trans, int *m, int *n, int *nrhs, complex<double> *a,
           int *lda, complex<double> *b, int *ldb,
           complex<double> *work, int *lwork, int *info);

    void
    sgelss_(int *m, int *n, int *nrhs,
            float *a, int *lda, float *b, int *ldb,
            float *s, float *rcond, int *rank,
            float *work, int *lwork, int *info);

    void
    dgelss_(int *m, int *n, int *nrhs,
            double *a, int *lda, double *b, int *ldb,
            double *s, double *rcond, int *rank,
            double *work, int *lwork, int *info);

    void
    cgelss_(int *m, int *n, int *nrhs,
            complex<float> *a, int *lda,
            complex<float> *b, int *ldb,
            complex<float> *s, complex<float> *rcond, int *rank,
            complex<float> *work, int *lwork, int *info);

    void
    zgelss_(int *m, int *n, int *nrhs,
            complex<double> *a, int *lda,
            complex<double> *b, int *ldb,
            complex<double> *s, complex<double> *rcond, int *rank,
            complex<double> *work, int *lwork, int *info);

    void
    sgees_(char * jobvs, char * sort, sgees_select *select,
           int *n, float *a, int *lda, int *sdim, float *wr, float *wi,
           float *vs, int *ldvs, float *work, int *lwork, int *bwork,
           int *info);

    void
    dgees_(char * jobvs, char * sort, dgees_select *select,
           int *n, double *a, int *lda, int *sdim, double *wr, double *wi,
           double *vs, int *ldvs, double *work, int *lwork, int *bwork,
           int *info);

    void
    cgees_(char * jobvs, char * sort, cgees_select *select,
           int *n, complex<float> *a, int *lda, int *sdim, complex<float> *w,
           complex<float> *vs, int *ldvs,
           complex<float> *work, int *lwork, float *rwork, int *bwork,
           int *info);

    void
    zgees_(char * jobvs, char * sort, zgees_select *select,
           int *n, complex<double> *a, int *lda, int *sdim, complex<double> *w,
           complex<double> *vs, int *ldvs,
           complex<double> *work, int *lwork, double *rwork, int *bwork,
           int *info);

    void
    sgeev_(char *jobvl, char *jobvr, int *n, float *a, int *lda,
           float *wr, float *wi,
           float *vl, int *ldvl,
           float *vr, int *ldvr,
           float *work, int *lwork, int *info);

    void
    dgeev_(char *jobvl, char *jobvr, int *n, double *a, int *lda,
           double *wr, double *wi,
           double *vl, int *ldvl,
           double *vr, int *ldvr,
           double *work, int *lwork, int *info);

    void
    cgeev_(char *jobvl, char *jobvr, int *n, complex<float> *a, int *lda,
           complex<float> *w,
           complex<float> *vl, int *ldvl,
           complex<float> *vr, int *ldvr,
           complex<float> *work,int *lwork,float *rwork,int *info);

    void
    zgeev_(char *jobvl, char *jobvr, int *n, complex<double> *a, int *lda,
           complex<double> *w,
           complex<double> *vl, int *ldvl,
           complex<double> *vr, int *ldvr,
           complex<double> *work,int *lwork,double *rwork,int *info);

    void
    ssyev_(char *jobz, char *uplo, int *n, float *a, int *lda,
           float *w, float *work, int *lwork, int *info);

    void
    dsyev_(char *jobz, char *uplo, int *n, double *a, int *lda,
           double *w, double *work, int *lwork, int *info);

    void
    ssbev_(char *jobz, char *upLo, int *n, int *kd, float *ab, int *ldab,
           float *w, float *z, int *ldz, float *work, int *info);

    void
    dsbev_(char *jobz, char *upLo, int *n, int *kd, double *ab, int *ldab,
           double *w, double *z, int *ldz, double *work, int *info);

    void
    sspev_(char *jobz, char *upLo, int *n, float *ap, float *w,
           float *z, int *ldz, float *work, int *info);

    void
    dspev_(char *jobz, char *upLo, int *n, double *ap, double *w,
           double *z, int *ldz, double *work, int *info);

    void
    cheev_(char *jobz, char *uplo, int *n, complex<float> *a, int *lda,
           float *w, complex<float> *work, int *lwork, float *rwork,
           int *info);

    void
    zheev_(char *jobz, char *uplo, int *n, complex<double> *a, int *lda,
           double *w, complex<double> *work, int *lwork, double *rwork,
           int *info);

    int
    chbev_(char *jobz, char *upLo, int *n, int *kd, complex<float> *ab,
           int *ldab, float *w, complex<float> *Z, int *ldz,
           complex<float> *work, float *rwork, int *info);

    int
    zhbev_(char *jobz, char *upLo, int *n, int *kd, complex<double> *ab,
           int *ldab, double *w, complex<double> *Z, int *ldz,
           complex<double> *work, double *rwork, int *info);

    int
    chpev_(char *jobz, char *upLo, int *n, complex<float> *ap,
           float *w, complex<float> *Z, int *ldz,
           complex<float> *work, float *rwork, int *info);

    int
    zhpev_(char *jobz, char *upLo, int *n, complex<double> *ap,
           double *w, complex<double> *Z, int *ldz,
           complex<double> *work, double *rwork, int *info);

    void
    sgesvd_(char *jobu, char *jobvt, int *m, int *n, float *a, int *lda,
            float *s, float *u, int *ldu, float *vt, int *ldvt,
            float *work, int *lwork, int *info);

    void
    dgesvd_(char *jobu, char *jobvt, int *m, int *n, double *a, int *lda,
            double *s, double *u, int *ldu, double *vt, int *ldvt,
            double *work, int *lwork, int *info);

    void
    cgesvd_(char *jobu, char *jobvt,
            int *m, int *n, complex<float> *a, int *lda, float *s,
            complex<float> *u, int *ldu,
            complex<float> *vt, int *ldvt,
            complex<float> *work, int *lwork, float *rwork, int *info);

    void
    zgesvd_(char *jobu, char *jobvt,
            int *m, int *n, complex<double> *a, int *lda, double *s,
            complex<double> *u, int *ldu,
            complex<double> *vt, int *ldvt,
            complex<double> *work, int *lwork, double *rwork, int *info);

    void
    sgesdd_(char *jobz, int *m, int *n, float *a, int *lda,
            float *s, float *u, int *ldu, float *vt, int *ldvt,
            float *work, int *lwork, int *iwork, int *info);

    void
    dgesdd_(char *jobz, int *m, int *n, double *a, int *lda,
            double *s, double *u, int *ldu, double *vt, int *ldvt,
            double *work, int *lwork, int *iwork, int *info);

    void
    cgesdd_(char *jobz,
            int *m, int *n, complex<float> *a, int *lda, float *s,
            complex<float> *u, int *ldu,
            complex<float> *vt, int *ldvt,
            complex<float> *work, int *lwork, float *rwork, int *iwork,
            int *info);

    void
    zgesdd_(char *jobzu,
            int *m, int *n, complex<double> *a, int *lda, double *s,
            complex<double> *u, int *ldu,
            complex<double> *vt, int *ldvt,
            complex<double> *work, int *lwork, double *rwork, int *iwork,
            int *info);

    void
    dgecon_(char *norm, int *n, double *a, int *lda, double *anorm,
            double *rcond, double *work, int *iwork, int *info);
}

int
potrf(StorageUpLo upLo, int n, double *a, int lda)
{
    int info;
    char _upLo = (upLo==Upper) ? 'U' : 'L';
    dpotrf_(&_upLo,&n,a,&lda,&info);
    return info;
}

int
getrf(int m, int n, float *a, int lda, int *ipiv)
{
    int info;
    sgetrf_(&m, &n, a, &lda, ipiv, &info);
    return info;
}

int
getrf(int m, int n, double *a, int lda, int *ipiv)
{
    int info;
    dgetrf_(&m, &n, a, &lda, ipiv, &info);
    return info;
}

int
getrf(int m, int n, complex<float> *a, int lda, int *ipiv)
{
    int info;
    cgetrf_(&m, &n, a, &lda, ipiv, &info);
    return info;
}

int
getrf(int m, int n, complex<double> *a, int lda, int *ipiv)
{
    int info;
    zgetrf_(&m, &n, a, &lda, ipiv, &info);
    return info;
}

int
gbtrf(int m, int n, int kl, int ku, float *ab, int ldab, int *ipiv)
{
    int info;
    sgbtrf_(&m, &n, &kl, &ku, ab, &ldab, ipiv, &info);
    return info;
}

int
gbtrf(int m, int n, int kl, int ku, double *ab, int ldab, int *ipiv)
{
    int info;
    dgbtrf_(&m, &n, &kl, &ku, ab, &ldab, ipiv, &info);
    return info;
}

int
getri(int n, float *a, int lda, const int *ipiv,
      float *work, int lwork)
{
    int info;
    sgetri_(&n, a, &lda, ipiv, work, &lwork, &info);
    return info;
}

int
getri(int n, double *a, int lda, const int *ipiv,
      double *work, int lwork)
{
    int info;
    dgetri_(&n, a, &lda, ipiv, work, &lwork, &info);
    return info;
}

int
getri(int n, complex<float> *a, int lda, const int *ipiv,
      complex<float> *work, int lwork)
{
    int info;
    cgetri_(&n, a, &lda, ipiv, work, &lwork, &info);
    return info;
}


int
getri(int n, complex<double> *a, int lda, const int *ipiv,
      complex<double> *work, int lwork)
{
    int info;
    zgetri_(&n, a, &lda, ipiv, work, &lwork, &info);
    return info;
}

int
getrs(Transpose trans, int n, int nrhs, const float *a, int lda,
      const int *ipiv, float *b, int ldb)
{
    int info;
    char _trans = (trans==NoTrans) ? 'N' : 'T';

    sgetrs_(&_trans, &n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
    return info;
}

int
getrs(Transpose trans, int n, int nrhs, const double *a, int lda,
      const int *ipiv, double *b, int ldb)
{
    int info;
    char _trans = (trans==NoTrans) ? 'N' : 'T';

    dgetrs_(&_trans, &n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
    return info;
}

int
gbtrs(Transpose trans, int n, int kl, int ku, int nrhs,
      const float *ab, int ldab, const int *ipiv, float *b, int ldb)
{
    int info;
    char _trans = (trans==NoTrans) ? 'N' : 'T';
    sgbtrs_(&_trans, &n, &kl, &ku, &nrhs, ab, &ldab, ipiv, b, &ldb, &info);
    return info;
}

int
gbtrs(Transpose trans, int n, int kl, int ku, int nrhs,
      const double *ab, int ldab, const int *ipiv, double *b, int ldb)
{
    int info;
    char _trans = (trans==NoTrans) ? 'N' : 'T';
    dgbtrs_(&_trans, &n, &kl, &ku, &nrhs, ab, &ldab, ipiv, b, &ldb, &info);
    return info;
}

int
gesv(int n, int nrhs, float *a, int lda, int *ipiv, float *b, int ldb)
{
    int info;
    sgesv_(&n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
    return info;
}

int
gesv(int n, int nrhs, double *a, int lda, int *ipiv, double *b, int ldb)
{
    int info;
    dgesv_(&n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
    return info;
}

int
gesv(int n, int nrhs, complex<float> *a, int lda, int *ipiv,
     complex<float> *b, int ldb)
{
    int info;
    cgesv_(&n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
    return info;
}

int
gesv(int n, int nrhs, complex<double> *a, int lda, int *ipiv,
     complex<double> *b, int ldb)
{
    int info;
    zgesv_(&n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
    return info;
}

int
gbsv(int n, int kl, int ku, int nrhs, float *ab, int ldab,
     int *ipiv, float *b, int ldb)
{
    int info;
    sgbsv_(&n, &kl, &ku, &nrhs, ab, &ldab, ipiv, b, &ldb, &info);
    return info;
}

int
gbsv(int n, int kl, int ku, int nrhs, double *ab, int ldab,
     int *ipiv, double *b, int ldb)
{
    int info;
    dgbsv_(&n, &kl, &ku, &nrhs, ab, &ldab, ipiv, b, &ldb, &info);
    return info;
}

int
gbsv(int n, int kl, int ku, int nrhs, complex<float> *ab, int ldab,
     int *ipiv, complex<float> *b, int ldb)
{
    int info;
    cgbsv_(&n, &kl, &ku, &nrhs, ab, &ldab, ipiv, b, &ldb, &info);
    return info;
}

int
gbsv(int n, int kl, int ku, int nrhs, complex<double> *ab, int ldab,
     int *ipiv, complex<double> *b, int ldb)
{
    int info;
    zgbsv_(&n, &kl, &ku, &nrhs, ab, &ldab, ipiv, b, &ldb, &info);
    return info;
}

int
trtrs(StorageUpLo upLo, Transpose trans, UnitDiag diag, int n, int nrhs,
      const float *a, int lda, float *b, int ldb)
{
    int info;
    char _upLo = (upLo==Upper) ? 'U' : 'L';
    char _trans = (trans==NoTrans) ? 'N' : 'T';
    char _diag = (diag==Unit) ? 'U' : 'N';

    strtrs_(&_upLo, &_trans, &_diag, &n, &nrhs, a, &lda, b, &ldb, &info);
    return info;
}

int
trtrs(StorageUpLo upLo, Transpose trans, UnitDiag diag, int n, int nrhs,
      const double *a, int lda, double *b, int ldb)
{
    int info;
    char _upLo = (upLo==Upper) ? 'U' : 'L';
    char _trans = (trans==NoTrans) ? 'N' : 'T';
    char _diag = (diag==Unit) ? 'U' : 'N';

    dtrtrs_(&_upLo, &_trans, &_diag, &n, &nrhs, a, &lda, b, &ldb, &info);
    return info;
}

int
geqrf(int m, int n, float *a, int lda, float *tau, float *work, int lwork)
{
    int info;
    sgeqrf_(&m, &n, a, &lda, tau, work, &lwork, &info);
    return info;
}

int
geqrf(int m, int n, double *a, int lda, double *tau, double *work, int lwork)
{
    int info;
    dgeqrf_(&m, &n, a, &lda, tau, work, &lwork, &info);
    return info;
}

int
orgqr(int m, int n, int k, float *a, int lda, const float *tau,
      float *work, int lwork)
{
    int info;
    sorgqr_(&m, &n, &k, a, &lda, tau, work, &lwork, &info);
    return info;
}

int
orgqr(int m, int n, int k, double *a, int lda, const double *tau,
      double *work, int lwork)
{
    int info;
    dorgqr_(&m, &n, &k, a, &lda, tau, work, &lwork, &info);
    return info;
}

int
ormqr(BlasSide side, Transpose trans, int m, int n, int k,
      const float *a, int lda, const float *tau, float *c, int ldc,
      float *work, int lwork)
{
    int info;
    char _side = (side==Left) ? 'L' : 'R';
    char _trans = (trans==NoTrans) ? 'N' : 'T';

    sormqr_(&_side, &_trans, &m, &n, &k, a, &lda, tau,
            c, &ldc, work, &lwork, &info);
    return info;
}

int
ormqr(BlasSide side, Transpose trans, int m, int n, int k,
      const double *a, int lda, const double *tau, double *c, int ldc,
      double *work, int lwork)
{
    int info;
    char _side = (side==Left) ? 'L' : 'R';
    char _trans = (trans==NoTrans) ? 'N' : 'T';

    dormqr_(&_side, &_trans, &m, &n, &k, a, &lda, tau,
            c, &ldc, work, &lwork, &info);
    return info;
}

int
gels(Transpose trans, int m, int n, int nrhs, float *a, int lda,
     float *b, int ldb, float *work, int lwork)
{
    int info;
    char _trans = (trans==NoTrans) ? 'N' : 'T';

    sgels_(&_trans, &m, &n, &nrhs, a, &lda, b, &ldb, work, &lwork, &info);
    return info;
}

int
gels(Transpose trans, int m, int n, int nrhs, double *a, int lda,
     double *b, int ldb, double *work, int lwork)
{
    int info;
    char _trans = (trans==NoTrans) ? 'N' : 'T';

    dgels_(&_trans, &m, &n, &nrhs, a, &lda, b, &ldb, work, &lwork, &info);
    return info;
}

int
gels(Transpose trans, int m, int n, int nrhs, complex<float> *a, int lda,
     complex<float> *b, int ldb, complex<float> *work, int lwork)
{
    int info;
    char _trans = (trans==NoTrans) ? 'N' : 'T';

    cgels_(&_trans, &m, &n, &nrhs, a, &lda, b, &ldb, work, &lwork, &info);
    return info;
}

int
gels(Transpose trans, int m, int n, int nrhs, complex<double> *a, int lda,
     complex<double> *b, int ldb, complex<double> *work, int lwork)
{
    int info;
    char _trans = (trans==NoTrans) ? 'N' : 'T';

    zgels_(&_trans, &m, &n, &nrhs, a, &lda, b, &ldb, work, &lwork, &info);
    return info;
}

int
gelss(int m, int n, int nrhs, float *a, int lda, float *b, int ldb,
     float *s, float rcond, int rank, float *work, int lwork)
{
    int info;
    sgelss_(&m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, &rank,
            work, &lwork, &info);
    return info;
}

int
gelss(int m, int n, int nrhs, double *a, int lda, double *b, int ldb,
     double *s, double rcond, int rank, double *work, int lwork)
{
    int info;
    dgelss_(&m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, &rank,
            work, &lwork, &info);
    return info;
}

int
gelss(int m, int n, int nrhs, complex<float> *a, int lda,
      complex<float> *b, int ldb, complex<float> *s,
      complex<float> rcond, int rank, complex<float> *work,
      int lwork)
{
    int info;
    cgelss_(&m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, &rank,
            work, &lwork, &info);
    return info;
}

int
gelss(int m, int n, int nrhs, complex<double> *a, int lda,
      complex<double> *b, int ldb, complex<double> *s,
      complex<double> rcond, int rank, complex<double> *work,
      int lwork)
{
    int info;
    zgelss_(&m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, &rank,
            work, &lwork, &info);
    return info;
}

int
geev(bool jobvl, bool jobvr, int n, float *a, int lda,
     float *wr, float *wi,
     float *vl, int ldvl,
     float *vr, int ldvr,
     float *work, int lwork)
{
    int info;
    char _jobvl = (jobvl==true) ? 'V' : 'N';
    char _jobvr = (jobvr==true) ? 'V' : 'N';

    sgeev_(&_jobvl, &_jobvr, &n, a, &lda, wr, wi, vl, &ldvl, vr, &ldvr,
           work, &lwork, &info);
    return info;
}

//- syev
int
syev(bool jobz, StorageUpLo upLo, int n, float *a, int lda,
     float *w, float *work, int lwork)
{
    int info;
    char _jobz = (jobz==true) ? 'V' : 'N';
    char _upLo = (upLo==Upper) ? 'U' : 'L';

    ssyev_(&_jobz, &_upLo, &n, a, &lda, w, work, &lwork, &info);
    return info;
}

int
syev(bool jobz, StorageUpLo upLo, int n, double *a, int lda,
     double *w, double *work, int lwork)
{
    int info;
    char _jobz = (jobz==true) ? 'V' : 'N';
    char _upLo = (upLo==Upper) ? 'U' : 'L';

    dsyev_(&_jobz, &_upLo, &n, a, &lda, w, work, &lwork, &info);
    return info;
}

//- sbev
int
sbev(bool jobz, StorageUpLo upLo, int n, int kd, float *ab, int ldab,
     float *w, float *z, int ldz, float *work)
{
    int info;
    char _jobz = (jobz==true) ? 'V' : 'N';
    char _upLo = (upLo==Upper) ? 'U' : 'L';

    ssbev_(&_jobz, &_upLo, &n, &kd, ab, &ldab, w, z, &ldz, work, &info);
    return info;
}

int
sbev(bool jobz, StorageUpLo upLo, int n, int kd, double *ab, int ldab,
     double *w, double *z, int ldz, double *work)
{
    int info;
    char _jobz = (jobz==true) ? 'V' : 'N';
    char _upLo = (upLo==Upper) ? 'U' : 'L';

    dsbev_(&_jobz, &_upLo, &n, &kd, ab, &ldab, w, z, &ldz, work, &info);
    return info;
}

//- spev
int
spev(bool jobz, StorageUpLo upLo, int n, float *ap, float *w,
     float *z, int ldz, float *work)
{
    int info;
    char _jobz = (jobz==true) ? 'V' : 'N';
    char _upLo = (upLo==Upper) ? 'U' : 'L';

    sspev_(&_jobz, &_upLo, &n, ap, w, z, &ldz, work, &info);
    return info;
}

int
spev(bool jobz, StorageUpLo upLo, int n, double *ap, double *w,
     double *z, int ldz, double *work)
{
    int info;
    char _jobz = (jobz==true) ? 'V' : 'N';
    char _upLo = (upLo==Upper) ? 'U' : 'L';

    dspev_(&_jobz, &_upLo, &n, ap, w, z, &ldz, work, &info);
    return info;
}

//- heev
int
heev(bool jobz, StorageUpLo upLo, int n, complex<float> *a, int lda, float *w,
     complex<float> *work, int lwork, float *rwork )
{
    int info;
    char _jobz = (jobz==true) ? 'V' : 'N';
    char _upLo = (upLo==Upper) ? 'U' : 'L';

    cheev_(&_jobz, &_upLo, &n, a, &lda, w, work, &lwork, rwork, &info);

    return info;
}

int
heev(bool jobz, StorageUpLo upLo, int n, complex<double> *a, int lda, double *w,
     complex<double> *work, int lwork, double *rwork )
{
    int info;
    char _jobz = (jobz==true) ? 'V' : 'N';
    char _upLo = (upLo==Upper) ? 'U' : 'L';

    zheev_(&_jobz, &_upLo, &n, a, &lda, w, work, &lwork, rwork, &info);

    return info;
}

//- hbev
int
hbev(bool jobz, StorageUpLo upLo, int n, int kd, complex<float> *ab, int ldab,
     float *w, complex<float> *Z, int ldz,
     complex<float> *work, float *rwork)
{
    int info;
    char _jobz = (jobz==true) ? 'V' : 'N';
    char _upLo = (upLo==Upper) ? 'U' : 'L';

    chbev_(&_jobz, &_upLo, &n, &kd, ab, &ldab, w, Z, &ldz, work, rwork, &info);
    return info;
}

int
hbev(bool jobz, StorageUpLo upLo, int n, int kd, complex<double> *ab, int ldab,
     double *w, complex<double> *Z, int ldz,
     complex<double> *work, double *rwork)
{
    int info;
    char _jobz = (jobz==true) ? 'V' : 'N';
    char _upLo = (upLo==Upper) ? 'U' : 'L';

    zhbev_(&_jobz, &_upLo, &n, &kd, ab, &ldab, w, Z, &ldz, work, rwork, &info);
    return info;
}

//- hpev
int
hpev(bool jobz, StorageUpLo upLo, int n, complex<float> *ap, float *w,
     complex<float> *Z, int ldz, complex<float> *work, float *rwork)
{
    int info;
    char _jobz = (jobz==true) ? 'V' : 'N';
    char _upLo = (upLo==Upper) ? 'U' : 'L';

    chpev_(&_jobz, &_upLo, &n, ap, w, Z, &ldz, work, rwork, &info);
    return info;
}

int
hpev(bool jobz, StorageUpLo upLo, int n, complex<double> *ap, double *w,
     complex<double> *Z, int ldz, complex<double> *work, double *rwork)
{
    int info;
    char _jobz = (jobz==true) ? 'V' : 'N';
    char _upLo = (upLo==Upper) ? 'U' : 'L';

    zhpev_(&_jobz, &_upLo, &n, ap, w, Z, &ldz, work, rwork, &info);
    return info;
}

int
gees(bool jobvs, bool sort, sgees_select *select,
     int n, float *a, int lda, int &sdim, float *wr, float *wi,
     float *vs, int ldvs, float *work, int lwork, int *bwork)
{
    int info;
    char _jobvs = (jobvs==true) ? 'V' : 'N';
    char _sort  = (sort==true) ? 'S' : 'N';

    sgees_(&_jobvs, &_sort, select, &n, a, &lda, &sdim,
           wr, wi, vs, &ldvs,
           work, &lwork, bwork, &info);
    return info;
}

int
gees(bool jobvs, bool sort, dgees_select *select,
     int n, double *a, int lda, int &sdim, double *wr, double *wi,
     double *vs, int ldvs, double *work, int lwork, int *bwork)
{
    int info;
    char _jobvs = (jobvs==true) ? 'V' : 'N';
    char _sort  = (sort==true) ? 'S' : 'N';

    dgees_(&_jobvs, &_sort, select, &n, a, &lda, &sdim,
           wr, wi, vs, &ldvs,
           work, &lwork, bwork, &info);
    return info;
}

int
gees(bool jobvs, bool sort, cgees_select *select,
     int n, complex<float> *a, int lda, int &sdim, complex<float> *w,
     complex<float> *vs, int ldvs,
     complex<float> *work, int lwork, float *rwork, int *bwork)
{
    int info;
    char _jobvs = (jobvs==true) ? 'V' : 'N';
    char _sort  = (sort==true) ? 'S' : 'N';

    cgees_(&_jobvs, &_sort, select, &n, a, &lda, &sdim,
           w, vs, &ldvs,
           work, &lwork, rwork, bwork, &info);
    return info;
}

int
gees(bool jobvs, bool sort, zgees_select *select,
     int n, complex<double> *a, int lda, int &sdim, complex<double> *w,
     complex<double> *vs, int ldvs,
     complex<double> *work, int lwork, double *rwork, int *bwork)
{
    int info;
    char _jobvs = (jobvs==true) ? 'V' : 'N';
    char _sort  = (sort==true) ? 'S' : 'N';

    zgees_(&_jobvs, &_sort, select, &n, a, &lda, &sdim,
           w, vs, &ldvs,
           work, &lwork, rwork, bwork, &info);
    return info;
}

int
geev(bool jobvl, bool jobvr, int n, double *a, int lda,
     double *wr, double *wi,
     double *vl, int ldvl,
     double *vr, int ldvr,
     double *work, int lwork)
{
    int info;
    char _jobvl = (jobvl==true) ? 'V' : 'N';
    char _jobvr = (jobvr==true) ? 'V' : 'N';

    dgeev_(&_jobvl, &_jobvr, &n, a, &lda, wr, wi, vl, &ldvl, vr, &ldvr,
           work, &lwork, &info);
    return info;
}

int
geev(bool jobvl, bool jobvr, int n, complex<float> *a, int lda,
     complex<float> *w,
     complex<float> *vl, int ldvl,
     complex<float> *vr, int ldvr,
     complex<float> *work, int lwork, float *rwork)
{
    int info;
    char _jobvl = (jobvl==true) ? 'V' : 'N';
    char _jobvr = (jobvr==true) ? 'V' : 'N';

    cgeev_(&_jobvl, &_jobvr, &n, a, &lda, w, vl, &ldvl, vr, &ldvr,
           work, &lwork, rwork, &info);
    return info;
}

int
geev(bool jobvl, bool jobvr, int n, complex<double> *a, int lda,
     complex<double> *w,
     complex<double> *vl, int ldvl,
     complex<double> *vr, int ldvr,
     complex<double> *work, int lwork, double *rwork)
{
    int info;
    char _jobvl = (jobvl==true) ? 'V' : 'N';
    char _jobvr = (jobvr==true) ? 'V' : 'N';

    zgeev_(&_jobvl, &_jobvr, &n, a, &lda, w, vl, &ldvl, vr, &ldvr,
           work, &lwork, rwork, &info);
    return info;
}

static char jobchar[4] = {'A','S','O','N'};

int
gesvd(SVectorsJob jobu, SVectorsJob jobvt,
      int m, int n, float *a, int lda,                // A
      float *s,                                      // singular values
      float *u, int ldu,                         // left singular vectors
      float *vt, int ldvt,                       // right singular vectors
      float *work, int lwork)
{
    int info;
    sgesvd_(&jobchar[jobu], &jobchar[jobvt],
            &m, &n, a, &lda, s, u, &ldu, vt, &ldvt,
            work, &lwork, &info);
    return info;
}

int
gesvd(SVectorsJob jobu, SVectorsJob jobvt,
      int m, int n, double *a, int lda,                // A
      double *s,                                      // singular values
      double *u, int ldu,                         // left singular vectors
      double *vt, int ldvt,                       // right singular vectors
      double *work, int lwork)
{
    int info;
    dgesvd_(&jobchar[jobu], &jobchar[jobvt],
            &m, &n, a, &lda, s, u, &ldu, vt, &ldvt,
            work, &lwork, &info);
    return info;
}

int
gesvd(SVectorsJob jobu, SVectorsJob jobvt,
      int m, int n, complex<float> *a, int lda,       // A
      float *s,                                      // singular values
      complex<float> *u, int ldu,                // left singular vectors
      complex<float> *vt, int ldvt,              // right singular vectors
      complex<float> *work, int lwork, float *rwork)
{
    int info;
    cgesvd_(&jobchar[jobu], &jobchar[jobvt],
            &m, &n, a, &lda, s, u, &ldu, vt, &ldvt,
            work, &lwork, rwork, &info);
    return info;
}

int
gesvd(SVectorsJob jobu, SVectorsJob jobvt,
      int m, int n, complex<double> *a, int lda,       // A
      double *s,                                      // singular values
      complex<double> *u, int ldu,                // left singular vectors
      complex<double> *vt, int ldvt,              // right singular vectors
      complex<double> *work, int lwork, double *rwork)
{
    int info;
    zgesvd_(&jobchar[jobu], &jobchar[jobvt],
            &m, &n, a, &lda, s, u, &ldu, vt, &ldvt,
            work, &lwork, rwork, &info);
    return info;
}

int
gesdd(SVectorsJob jobz,
      int m, int n, float *a, int lda,                // A
      float *s,                                      // singular values
      float *u, int ldu,                         // left singular vectors
      float *vt, int ldvt,                       // right singular vectors
      float *work, int lwork, int *iwork)
{
    int info;
    sgesdd_(&jobchar[jobz],
            &m, &n, a, &lda, s, u, &ldu, vt, &ldvt,
            work, &lwork, iwork, &info);
    return info;
}

int
gesdd(SVectorsJob jobz,
      int m, int n, double *a, int lda,                // A
      double *s,                                      // singular values
      double *u, int ldu,                         // left singular vectors
      double *vt, int ldvt,                       // right singular vectors
      double *work, int lwork, int *iwork)
{
    int info;
    dgesdd_(&jobchar[jobz],
            &m, &n, a, &lda, s, u, &ldu, vt, &ldvt,
            work, &lwork, iwork, &info);
    return info;
}

int
gesdd(SVectorsJob jobz,
      int m, int n, complex<float> *a, int lda,       // A
      float *s,                                      // singular values
      complex<float> *u, int ldu,                // left singular vectors
      complex<float> *vt, int ldvt,              // right singular vectors
      complex<float> *work, int lwork, float *rwork, int *iwork)
{
    int info;
    cgesdd_(&jobchar[jobz],
            &m, &n, a, &lda, s, u, &ldu, vt, &ldvt,
            work, &lwork, rwork, iwork, &info);
    return info;
}

int
gesdd(SVectorsJob jobz,
      int m, int n, complex<double> *a, int lda,       // A
      double *s,                                      // singular values
      complex<double> *u, int ldu,                // left singular vectors
      complex<double> *vt, int ldvt,              // right singular vectors
      complex<double> *work, int lwork, double *rwork, int *iwork)
{
    int info;
    zgesdd_(&jobchar[jobz],
            &m, &n, a, &lda, s, u, &ldu, vt, &ldvt,
            work, &lwork, rwork, iwork, &info);
    return info;
}

int
gecon(char norm,
      int n, double *a, int lda,                  // A
      double *anorm,                              // the norm of A
      double *rcond,                              // reciprocal condition number
      double *work, int *iwork)
{
    int info;
    dgecon_(&norm, &n, a, &lda, anorm, rcond, work, iwork, &info);
    return info;
}


} // namespace flens
