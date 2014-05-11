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

#include <flens/lapack/auxiliary/getf77char.h>
#include <flens/lapack/auxiliary/nint.h>
#include <flens/lapack/auxiliary/sign.h>

#include <flens/lapack/debug/hex.h>
#include <flens/lapack/debug/isidentical.h>

#include <flens/lapack/gb/sv.h>
#include <flens/lapack/gb/tf2.h>
#include <flens/lapack/gb/trf.h>
#include <flens/lapack/gb/trs.h>

#include <flens/lapack/ge/bak.h>
#include <flens/lapack/ge/bal.h>
#include <flens/lapack/ge/con.h>
#include <flens/lapack/ge/equ.h>
#include <flens/lapack/ge/es.h>
#include <flens/lapack/ge/esx.h>
#include <flens/lapack/ge/ev.h>
#include <flens/lapack/ge/evx.h>
#include <flens/lapack/ge/hd2.h>
#include <flens/lapack/ge/hrd.h>
#include <flens/lapack/ge/jsv.h>
#include <flens/lapack/ge/lq2.h>
#include <flens/lapack/ge/lqf.h>
#include <flens/lapack/ge/ls.h>
#include <flens/lapack/ge/lsy.h>
#include <flens/lapack/ge/qp3.h>
#include <flens/lapack/ge/qr2.h>
#include <flens/lapack/ge/qrf.h>
#include <flens/lapack/ge/qrs.h>
#include <flens/lapack/ge/rfs.h>
#include <flens/lapack/ge/rscl.h>
#include <flens/lapack/ge/sv.h>
#include <flens/lapack/ge/svd.h>
#include <flens/lapack/ge/svj.h>
#include <flens/lapack/ge/svj0.h>
#include <flens/lapack/ge/svj1.h>
#include <flens/lapack/ge/svx.h>
#include <flens/lapack/ge/tf2.h>
#include <flens/lapack/ge/trf.h>
#include <flens/lapack/ge/tri.h>
#include <flens/lapack/ge/trs.h>

#include <flens/lapack/hb/ev.h>

#include <flens/lapack/he/ev.h>
#include <flens/lapack/he/sv.h>
#include <flens/lapack/he/td2.h>
#include <flens/lapack/he/trd.h>
#include <flens/lapack/he/trf.h>
#include <flens/lapack/he/tri.h>
#include <flens/lapack/he/trs.h>

#include <flens/lapack/hp/ev.h>
#include <flens/lapack/hp/sv.h>
#include <flens/lapack/hp/trf.h>
#include <flens/lapack/hp/tri.h>
#include <flens/lapack/hp/trs.h>

#include <flens/lapack/impl/hseqr.h>
#include <flens/lapack/impl/iparmq.h>
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
#include <flens/lapack/impl/ormr3.h>
#include <flens/lapack/impl/ormrz.h>
#include <flens/lapack/impl/steqr.h>
#include <flens/lapack/impl/sterf.h>
#include <flens/lapack/impl/trevc.h>
#include <flens/lapack/impl/trexc.h>
#include <flens/lapack/impl/trsen.h>
#include <flens/lapack/impl/trsna.h>
#include <flens/lapack/impl/trsyl.h>
#include <flens/lapack/impl/tzrzf.h>
#include <flens/lapack/impl/ung2l.h>
#include <flens/lapack/impl/ung2r.h>
#include <flens/lapack/impl/unghr.h>
#include <flens/lapack/impl/ungl2.h>
#include <flens/lapack/impl/unglq.h>
#include <flens/lapack/impl/ungql.h>
#include <flens/lapack/impl/ungqr.h>
#include <flens/lapack/impl/ungtr.h>
#include <flens/lapack/impl/unm2r.h>
#include <flens/lapack/impl/unmhr.h>
#include <flens/lapack/impl/unml2.h>
#include <flens/lapack/impl/unmlq.h>
#include <flens/lapack/impl/unmqr.h>
#include <flens/lapack/impl/unmr3.h>
#include <flens/lapack/impl/unmrz.h>

#include <flens/lapack/la/ilaenv.h>
#include <flens/lapack/la/ilalc.h>
#include <flens/lapack/la/ilalr.h>
#include <flens/lapack/la/labad.h>
#include <flens/lapack/la/lacn2.h>
#include <flens/lapack/la/ladiv.h>
#include <flens/lapack/la/lae2.h>
#include <flens/lapack/la/laev2.h>
#include <flens/lapack/la/laexc.h>
#include <flens/lapack/la/lahqr.h>
#include <flens/lapack/la/lahr2.h>
#include <flens/lapack/la/laic1.h>
#include <flens/lapack/la/laln2.h>
#include <flens/lapack/la/lamch.h>
#include <flens/lapack/la/lan.h>
#include <flens/lapack/la/lanst.h>
#include <flens/lapack/la/lanv2.h>
#include <flens/lapack/la/lapy2.h>
#include <flens/lapack/la/lapy3.h>
#include <flens/lapack/la/laq.h>
#include <flens/lapack/la/laqp2.h>
#include <flens/lapack/la/laqps.h>
#include <flens/lapack/la/laqr0.h>
#include <flens/lapack/la/laqr1.h>
#include <flens/lapack/la/laqr2.h>
#include <flens/lapack/la/laqr3.h>
#include <flens/lapack/la/laqr4.h>
#include <flens/lapack/la/laqr5.h>
#include <flens/lapack/la/laqtr.h>
#include <flens/lapack/la/larf.h>
#include <flens/lapack/la/larfb.h>
#include <flens/lapack/la/larfg.h>
#include <flens/lapack/la/larft.h>
#include <flens/lapack/la/larfx.h>
#include <flens/lapack/la/lartg.h>
#include <flens/lapack/la/larz.h>
#include <flens/lapack/la/larzb.h>
#include <flens/lapack/la/larzt.h>
#include <flens/lapack/la/lascl.h>
#include <flens/lapack/la/lasr.h>
#include <flens/lapack/la/lasrt.h>
#include <flens/lapack/la/lassq.h>
#include <flens/lapack/la/laswp.h>
#include <flens/lapack/la/lasy2.h>
#include <flens/lapack/la/latrd.h>
#include <flens/lapack/la/latrs.h>
#include <flens/lapack/la/latrz.h>
#include <flens/lapack/la/lauu2.h>
#include <flens/lapack/la/lauum.h>

#include <flens/lapack/pb/pbsv.h>
#include <flens/lapack/pb/pbtrf.h>
#include <flens/lapack/pb/pbtrs.h>

#include <flens/lapack/po/pocon.h>
#include <flens/lapack/po/posv.h>
#include <flens/lapack/po/potf2.h>
#include <flens/lapack/po/potrf.h>
#include <flens/lapack/po/potri.h>
#include <flens/lapack/po/potrs.h>

#include <flens/lapack/sb/ev.h>
#include <flens/lapack/sp/ev.h>
#include <flens/lapack/sp/sv.h>
#include <flens/lapack/sp/trf.h>
#include <flens/lapack/sp/tri.h>
#include <flens/lapack/sp/trs.h>

#include <flens/lapack/sy/ev.h>
#include <flens/lapack/sy/sv.h>
#include <flens/lapack/sy/trf.h>
#include <flens/lapack/sy/tri.h>
#include <flens/lapack/sy/trs.h>

#include <flens/lapack/tb/trs.h>

#include <flens/lapack/tp/tri.h>
#include <flens/lapack/tp/trs.h>

#include <flens/lapack/tr/ti2.h>
#include <flens/lapack/tr/tri.h>
#include <flens/lapack/tr/trs.h>

#include <flens/lapack/typedefs.h>

#endif // FLENS_LAPACK_LAPACK_H
