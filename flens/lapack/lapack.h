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

#ifndef FLENS_LAPACK_LAPACK_H
#define FLENS_LAPACK_LAPACK_H 1

//
//  Control debugging of LAPACK functions
//
#ifdef LAPACK_DEBUG
#   define LAPACK_DEBUG_OUT(msg)    std::cerr << "debug: " << msg << std::endl
#else
#   define LAPACK_DEBUG_OUT(msg)
#endif

//
//  Select LAPACK preferences generic vs external implementation.
//
#ifdef ALWAYS_USE_CXXLAPACK
#   ifndef USE_CXXLAPACK
#       define USE_CXXLAPACK
#   endif
#   define LAPACK_SELECT  external
#else
#   define LAPACK_SELECT  generic
#endif

//
//  If an external LAPACK implementation is available include headers
//
#ifdef USE_CXXLAPACK
#   include <cxxlapack/cxxlapack.h>
#endif

#include <cmath>
#include <flens/lapack/typedefs.h>

#include <flens/lapack/debug/isidentical.h>
#include <flens/lapack/debug/hex.h>

#include <flens/lapack/auxiliary/getf77char.h>
#include <flens/lapack/auxiliary/nint.h>
#include <flens/lapack/auxiliary/pow.h>
#include <flens/lapack/auxiliary/sign.h>

#include <flens/lapack/impl/bak.h>
#include <flens/lapack/impl/bal.h>
#include <flens/lapack/impl/con.h>
#include <flens/lapack/impl/equ.h>
#include <flens/lapack/impl/es.h>
#include <flens/lapack/impl/esx.h>
#include <flens/lapack/impl/ev.h>
#include <flens/lapack/impl/evx.h>
#include <flens/lapack/impl/hd2.h>
#include <flens/lapack/impl/hrd.h>
#include <flens/lapack/impl/hseqr.h>
#include <flens/lapack/impl/ilaenv.h>
#include <flens/lapack/impl/ilalc.h>
#include <flens/lapack/impl/ilalr.h>
#include <flens/lapack/impl/iparmq.h>
#include <flens/lapack/impl/jsv.h>
#include <flens/lapack/impl/labad.h>
#include <flens/lapack/impl/lacn2.h>
#include <flens/lapack/impl/ladiv.h>
#include <flens/lapack/impl/laexc.h>
#include <flens/lapack/impl/lahqr.h>
#include <flens/lapack/impl/lahr2.h>
#include <flens/lapack/impl/laln2.h>
#include <flens/lapack/impl/lamch.h>
#include <flens/lapack/impl/lan.h>
#include <flens/lapack/impl/lanv2.h>
#include <flens/lapack/impl/lapy2.h>
#include <flens/lapack/impl/laq.h>
#include <flens/lapack/impl/laqp2.h>
#include <flens/lapack/impl/laqps.h>
#include <flens/lapack/impl/laqr0.h>
#include <flens/lapack/impl/laqr1.h>
#include <flens/lapack/impl/laqr2.h>
#include <flens/lapack/impl/laqr3.h>
#include <flens/lapack/impl/laqr4.h>
#include <flens/lapack/impl/laqr5.h>
#include <flens/lapack/impl/laqtr.h>
#include <flens/lapack/impl/larf.h>
#include <flens/lapack/impl/larfb.h>
#include <flens/lapack/impl/larfg.h>
#include <flens/lapack/impl/larft.h>
#include <flens/lapack/impl/larfx.h>
#include <flens/lapack/impl/lartg.h>
#include <flens/lapack/impl/lascl.h>
#include <flens/lapack/impl/lassq.h>
#include <flens/lapack/impl/laswp.h>
#include <flens/lapack/impl/lasy2.h>
#include <flens/lapack/impl/latrs.h>
#include <flens/lapack/impl/lauu2.h>
#include <flens/lapack/impl/lauum.h>
#include <flens/lapack/impl/lq2.h>
#include <flens/lapack/impl/lqf.h>
#include <flens/lapack/impl/ls.h>
#include <flens/lapack/impl/org2r.h>
#include <flens/lapack/impl/orghr.h>
#include <flens/lapack/impl/orgl2.h>
#include <flens/lapack/impl/orglq.h>
#include <flens/lapack/impl/orgqr.h>
#include <flens/lapack/impl/orm2r.h>
#include <flens/lapack/impl/ormhr.h>
#include <flens/lapack/impl/orml2.h>
#include <flens/lapack/impl/ormlq.h>
#include <flens/lapack/impl/ormqr.h>
#include <flens/lapack/impl/pocon.h>
#include <flens/lapack/impl/posv.h>
#include <flens/lapack/impl/potf2.h>
#include <flens/lapack/impl/potrf.h>
#include <flens/lapack/impl/potri.h>
#include <flens/lapack/impl/potrs.h>
#include <flens/lapack/impl/qp3.h>
#include <flens/lapack/impl/qr2.h>
#include <flens/lapack/impl/qrf.h>
#include <flens/lapack/impl/qrs.h>
#include <flens/lapack/impl/rfs.h>
#include <flens/lapack/impl/rscl.h>
#include <flens/lapack/impl/sv.h>
#include <flens/lapack/impl/svd.h>
#include <flens/lapack/impl/svj.h>
#include <flens/lapack/impl/svj0.h>
#include <flens/lapack/impl/svj1.h>
#include <flens/lapack/impl/svx.h>
#include <flens/lapack/impl/tf2.h>
#include <flens/lapack/impl/ti2.h>
#include <flens/lapack/impl/trevc.h>
#include <flens/lapack/impl/trexc.h>
#include <flens/lapack/impl/trf.h>
#include <flens/lapack/impl/tri.h>
#include <flens/lapack/impl/trs.h>
#include <flens/lapack/impl/trsen.h>
#include <flens/lapack/impl/trsna.h>
#include <flens/lapack/impl/trsyl.h>
#include <flens/lapack/impl/unglq.h>
#include <flens/lapack/impl/ungqr.h>
#include <flens/lapack/impl/unmlq.h>
#include <flens/lapack/impl/unmqr.h>

#endif // FLENS_LAPACK_LAPACK_H
