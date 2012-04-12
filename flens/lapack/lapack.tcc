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

#ifdef USE_CXXLAPACK
#   include <cxxlapack/cxxlapack.tcc>
#endif

#include <flens/lapack/debug/isidentical.tcc>
#include <flens/lapack/debug/hex.tcc>

#include <flens/lapack/auxiliary/getf77char.tcc>
#include <flens/lapack/auxiliary/nint.tcc>
#include <flens/lapack/auxiliary/pow.tcc>
#include <flens/lapack/auxiliary/sign.tcc>

#include <flens/lapack/impl/bak.tcc>
#include <flens/lapack/impl/bal.tcc>
#include <flens/lapack/impl/con.tcc>
#include <flens/lapack/impl/equ.tcc>
#include <flens/lapack/impl/es.tcc>
#include <flens/lapack/impl/esx.tcc>
#include <flens/lapack/impl/ev.tcc>
#include <flens/lapack/impl/evx.tcc>
#include <flens/lapack/impl/hd2.tcc>
#include <flens/lapack/impl/hrd.tcc>
#include <flens/lapack/impl/hseqr.tcc>
#include <flens/lapack/impl/ilaenv.tcc>
#include <flens/lapack/impl/ilalc.tcc>
#include <flens/lapack/impl/ilalr.tcc>
#include <flens/lapack/impl/iparmq.tcc>
#include <flens/lapack/impl/jsv.tcc>
#include <flens/lapack/impl/labad.tcc>
#include <flens/lapack/impl/lacn2.tcc>
#include <flens/lapack/impl/ladiv.tcc>
#include <flens/lapack/impl/laexc.tcc>
#include <flens/lapack/impl/lahqr.tcc>
#include <flens/lapack/impl/lahr2.tcc>
#include <flens/lapack/impl/laln2.tcc>
#include <flens/lapack/impl/lamch.tcc>
#include <flens/lapack/impl/lan.tcc>
#include <flens/lapack/impl/lanv2.tcc>
#include <flens/lapack/impl/lapy2.tcc>
#include <flens/lapack/impl/laq.tcc>
#include <flens/lapack/impl/laqp2.tcc>
#include <flens/lapack/impl/laqps.tcc>
#include <flens/lapack/impl/laqr0.tcc>
#include <flens/lapack/impl/laqr1.tcc>
#include <flens/lapack/impl/laqr2.tcc>
#include <flens/lapack/impl/laqr3.tcc>
#include <flens/lapack/impl/laqr4.tcc>
#include <flens/lapack/impl/laqr5.tcc>
#include <flens/lapack/impl/laqtr.tcc>
#include <flens/lapack/impl/larf.tcc>
#include <flens/lapack/impl/larfb.tcc>
#include <flens/lapack/impl/larfg.tcc>
#include <flens/lapack/impl/larft.tcc>
#include <flens/lapack/impl/larfx.tcc>
#include <flens/lapack/impl/lartg.tcc>
#include <flens/lapack/impl/lascl.tcc>
#include <flens/lapack/impl/lassq.tcc>
#include <flens/lapack/impl/laswp.tcc>
#include <flens/lapack/impl/lasy2.tcc>
#include <flens/lapack/impl/latrs.tcc>
#include <flens/lapack/impl/lauu2.tcc>
#include <flens/lapack/impl/lauum.tcc>
#include <flens/lapack/impl/lq2.tcc>
#include <flens/lapack/impl/lqf.tcc>
#include <flens/lapack/impl/ls.tcc>
#include <flens/lapack/impl/org2r.tcc>
#include <flens/lapack/impl/orghr.tcc>
#include <flens/lapack/impl/orgl2.tcc>
#include <flens/lapack/impl/orglq.tcc>
#include <flens/lapack/impl/orgqr.tcc>
#include <flens/lapack/impl/orm2r.tcc>
#include <flens/lapack/impl/ormhr.tcc>
#include <flens/lapack/impl/orml2.tcc>
#include <flens/lapack/impl/ormlq.tcc>
#include <flens/lapack/impl/ormqr.tcc>
#include <flens/lapack/impl/pocon.tcc>
#include <flens/lapack/impl/posv.tcc>
#include <flens/lapack/impl/potf2.tcc>
#include <flens/lapack/impl/potrf.tcc>
#include <flens/lapack/impl/potri.tcc>
#include <flens/lapack/impl/potrs.tcc>
#include <flens/lapack/impl/qp3.tcc>
#include <flens/lapack/impl/qr2.tcc>
#include <flens/lapack/impl/qrf.tcc>
#include <flens/lapack/impl/qrs.tcc>
#include <flens/lapack/impl/rfs.tcc>
#include <flens/lapack/impl/rscl.tcc>
#include <flens/lapack/impl/sv.tcc>
#include <flens/lapack/impl/svj.tcc>
#include <flens/lapack/impl/svj0.tcc>
#include <flens/lapack/impl/svj1.tcc>
#include <flens/lapack/impl/svx.tcc>
#include <flens/lapack/impl/tf2.tcc>
#include <flens/lapack/impl/ti2.tcc>
#include <flens/lapack/impl/trevc.tcc>
#include <flens/lapack/impl/trexc.tcc>
#include <flens/lapack/impl/trf.tcc>
#include <flens/lapack/impl/tri.tcc>
#include <flens/lapack/impl/trs.tcc>
#include <flens/lapack/impl/trsen.tcc>
#include <flens/lapack/impl/trsna.tcc>
#include <flens/lapack/impl/trsyl.tcc>

#endif // FLENS_LAPACK_LAPACK_TCC
