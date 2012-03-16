/*
 *   Copyright (c) 2011, Michael Lehn
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

#ifndef FLENS_LAPACK_LAPACK_TCC
#define FLENS_LAPACK_LAPACK_TCC 1

#include <flens/lapack/aux/con.tcc>
#include <flens/lapack/aux/convert.tcc>
#include <flens/lapack/aux/equ.tcc>
#include <flens/lapack/aux/ilaenv.tcc>
#include <flens/lapack/aux/ilalc.tcc>
#include <flens/lapack/aux/ilalr.tcc>
#include <flens/lapack/aux/iparmq.tcc>
#include <flens/lapack/aux/labad.tcc>
#include <flens/lapack/aux/lacn2.tcc>
#include <flens/lapack/aux/ladiv.tcc>
#include <flens/lapack/aux/laln2.tcc>
#include <flens/lapack/aux/lamch.tcc>
#include <flens/lapack/aux/lan.tcc>
#include <flens/lapack/aux/laq.tcc>
#include <flens/lapack/aux/larf.tcc>
#include <flens/lapack/aux/larfb.tcc>
#include <flens/lapack/aux/larfg.tcc>
#include <flens/lapack/aux/larft.tcc>
#include <flens/lapack/aux/larfx.tcc>
#include <flens/lapack/aux/lartg.tcc>
#include <flens/lapack/aux/lascl.tcc>
#include <flens/lapack/aux/latrs.tcc>
#include <flens/lapack/aux/lapy2.tcc>
#include <flens/lapack/aux/lassq.tcc>
#include <flens/lapack/aux/laswp.tcc>
#include <flens/lapack/aux/lasy2.tcc>
#include <flens/lapack/aux/nint.tcc>
#include <flens/lapack/aux/pocon.tcc>
#include <flens/lapack/aux/rscl.tcc>
#include <flens/lapack/aux/sign.tcc>

#include <flens/lapack/debug/isidentical.tcc>

#include <flens/lapack/eig/bak.tcc>
#include <flens/lapack/eig/bal.tcc>
#include <flens/lapack/eig/es.tcc>
#include <flens/lapack/eig/esx.tcc>
#include <flens/lapack/eig/ev.tcc>
#include <flens/lapack/eig/evx.tcc>
#include <flens/lapack/eig/hd2.tcc>
#include <flens/lapack/eig/hrd.tcc>
#include <flens/lapack/eig/hseqr.tcc>
#include <flens/lapack/eig/lahr2.tcc>
#include <flens/lapack/eig/laexc.tcc>
#include <flens/lapack/eig/lahqr.tcc>
#include <flens/lapack/eig/lanv2.tcc>
#include <flens/lapack/eig/laqr0.tcc>
#include <flens/lapack/eig/laqr1.tcc>
#include <flens/lapack/eig/laqr2.tcc>
#include <flens/lapack/eig/laqr3.tcc>
#include <flens/lapack/eig/laqr4.tcc>
#include <flens/lapack/eig/laqr5.tcc>
#include <flens/lapack/eig/laqtr.tcc>
#include <flens/lapack/eig/orghr.tcc>
#include <flens/lapack/eig/ormhr.tcc>
#include <flens/lapack/eig/trevc.tcc>
#include <flens/lapack/eig/trexc.tcc>
#include <flens/lapack/eig/trsen.tcc>
#include <flens/lapack/eig/trsna.tcc>
#include <flens/lapack/eig/trsyl.tcc>

#include <flens/lapack/gesv/lauu2.tcc>
#include <flens/lapack/gesv/lauum.tcc>
#include <flens/lapack/gesv/posv.tcc>
#include <flens/lapack/gesv/potf2.tcc>
#include <flens/lapack/gesv/potrf.tcc>
#include <flens/lapack/gesv/potri.tcc>
#include <flens/lapack/gesv/potrs.tcc>
#include <flens/lapack/gesv/rfs.tcc>
#include <flens/lapack/gesv/sv.tcc>
#include <flens/lapack/gesv/svx.tcc>
#include <flens/lapack/gesv/tf2.tcc>
#include <flens/lapack/gesv/trf.tcc>
#include <flens/lapack/gesv/ti2.tcc>
#include <flens/lapack/gesv/tri.tcc>
#include <flens/lapack/gesv/trs.tcc>

#include <flens/lapack/lq/lq2.tcc>
#include <flens/lapack/lq/lqf.tcc>
#include <flens/lapack/lq/orgl2.tcc>
#include <flens/lapack/lq/orglq.tcc>
#include <flens/lapack/lq/orml2.tcc>
#include <flens/lapack/lq/ormlq.tcc>

#include <flens/lapack/ls/ls.tcc>

#include <flens/lapack/qr/laqp2.tcc>
#include <flens/lapack/qr/laqps.tcc>
#include <flens/lapack/qr/org2r.tcc>
#include <flens/lapack/qr/orgqr.tcc>
#include <flens/lapack/qr/orm2r.tcc>
#include <flens/lapack/qr/ormqr.tcc>
#include <flens/lapack/qr/qp3.tcc>
#include <flens/lapack/qr/qr2.tcc>
#include <flens/lapack/qr/qrf.tcc>
#include <flens/lapack/qr/qrs.tcc>

#include <flens/lapack/svd/jsv.tcc>
#include <flens/lapack/svd/svj.tcc>
#include <flens/lapack/svd/svj0.tcc>
#include <flens/lapack/svd/svj1.tcc>

#if defined CHECK_CXXLAPACK || defined USE_NATIVE_ILAENV
#   include <flens/lapack/interface/include/cxxlapack.tcc>
#endif

#endif // FLENS_LAPACK_LAPACK_TCC
