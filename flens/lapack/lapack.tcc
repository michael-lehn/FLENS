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

#include <flens/lapack/auxiliary/getf77char.tcc>
#include <flens/lapack/auxiliary/nint.tcc>
#include <flens/lapack/auxiliary/sign.tcc>

#include <flens/lapack/debug/hex.tcc>
#include <flens/lapack/debug/isidentical.tcc>

#include <flens/lapack/gb/sv.tcc>
#include <flens/lapack/gb/tf2.tcc>
#include <flens/lapack/gb/trf.tcc>
#include <flens/lapack/gb/trs.tcc>

#include <flens/lapack/ge/bak.tcc>
#include <flens/lapack/ge/bal.tcc>
#include <flens/lapack/ge/con.tcc>
#include <flens/lapack/ge/equ.tcc>
#include <flens/lapack/ge/es.tcc>
#include <flens/lapack/ge/esx.tcc>
#include <flens/lapack/ge/ev.tcc>
#include <flens/lapack/ge/evx.tcc>
#include <flens/lapack/ge/hd2.tcc>
#include <flens/lapack/ge/hrd.tcc>
#include <flens/lapack/ge/jsv.tcc>
#include <flens/lapack/ge/lq2.tcc>
#include <flens/lapack/ge/lqf.tcc>
#include <flens/lapack/ge/ls.tcc>
#include <flens/lapack/ge/lsy.tcc>
#include <flens/lapack/ge/qp3.tcc>
#include <flens/lapack/ge/qr2.tcc>
#include <flens/lapack/ge/qrf.tcc>
#include <flens/lapack/ge/qrs.tcc>
#include <flens/lapack/ge/rfs.tcc>
#include <flens/lapack/ge/rscl.tcc>
#include <flens/lapack/ge/sv.tcc>
#include <flens/lapack/ge/svd.tcc>
#include <flens/lapack/ge/svj.tcc>
#include <flens/lapack/ge/svj0.tcc>
#include <flens/lapack/ge/svj1.tcc>
#include <flens/lapack/ge/svx.tcc>
#include <flens/lapack/ge/tf2.tcc>
#include <flens/lapack/ge/trf.tcc>
#include <flens/lapack/ge/tri.tcc>
#include <flens/lapack/ge/trs.tcc>

#include <flens/lapack/hb/ev.tcc>

#include <flens/lapack/he/ev.tcc>
#include <flens/lapack/he/sv.tcc>
#include <flens/lapack/he/trf.tcc>
#include <flens/lapack/he/tri.tcc>
#include <flens/lapack/he/trs.tcc>

#include <flens/lapack/hp/ev.tcc>
#include <flens/lapack/hp/sv.tcc>
#include <flens/lapack/hp/trf.tcc>
#include <flens/lapack/hp/tri.tcc>
#include <flens/lapack/hp/trs.tcc>

#include <flens/lapack/impl/hseqr.tcc>
#include <flens/lapack/impl/iparmq.tcc>
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
#include <flens/lapack/impl/ormr3.tcc>
#include <flens/lapack/impl/ormrz.tcc>
#include <flens/lapack/impl/trevc.tcc>
#include <flens/lapack/impl/trexc.tcc>
#include <flens/lapack/impl/trsen.tcc>
#include <flens/lapack/impl/trsna.tcc>
#include <flens/lapack/impl/trsyl.tcc>
#include <flens/lapack/impl/tzrzf.tcc>
#include <flens/lapack/impl/ung2r.tcc>
#include <flens/lapack/impl/unghr.tcc>
#include <flens/lapack/impl/ungl2.tcc>
#include <flens/lapack/impl/unglq.tcc>
#include <flens/lapack/impl/ungqr.tcc>
#include <flens/lapack/impl/unm2r.tcc>
#include <flens/lapack/impl/unmhr.tcc>
#include <flens/lapack/impl/unml2.tcc>
#include <flens/lapack/impl/unmlq.tcc>
#include <flens/lapack/impl/unmqr.tcc>
#include <flens/lapack/impl/unmr3.tcc>
#include <flens/lapack/impl/unmrz.tcc>

#include <flens/lapack/la/ilaenv.tcc>
#include <flens/lapack/la/ilalc.tcc>
#include <flens/lapack/la/ilalr.tcc>
#include <flens/lapack/la/labad.tcc>
#include <flens/lapack/la/lacn2.tcc>
#include <flens/lapack/la/ladiv.tcc>
#include <flens/lapack/la/laexc.tcc>
#include <flens/lapack/la/lahqr.tcc>
#include <flens/lapack/la/lahr2.tcc>
#include <flens/lapack/la/laic1.tcc>
#include <flens/lapack/la/laln2.tcc>
#include <flens/lapack/la/lamch.tcc>
#include <flens/lapack/la/lan.tcc>
#include <flens/lapack/la/lanv2.tcc>
#include <flens/lapack/la/lapy2.tcc>
#include <flens/lapack/la/lapy3.tcc>
#include <flens/lapack/la/laq.tcc>
#include <flens/lapack/la/laqp2.tcc>
#include <flens/lapack/la/laqps.tcc>
#include <flens/lapack/la/laqr0.tcc>
#include <flens/lapack/la/laqr1.tcc>
#include <flens/lapack/la/laqr2.tcc>
#include <flens/lapack/la/laqr3.tcc>
#include <flens/lapack/la/laqr4.tcc>
#include <flens/lapack/la/laqr5.tcc>
#include <flens/lapack/la/laqtr.tcc>
#include <flens/lapack/la/larf.tcc>
#include <flens/lapack/la/larfb.tcc>
#include <flens/lapack/la/larfg.tcc>
#include <flens/lapack/la/larft.tcc>
#include <flens/lapack/la/larfx.tcc>
#include <flens/lapack/la/lartg.tcc>
#include <flens/lapack/la/larz.tcc>
#include <flens/lapack/la/larzb.tcc>
#include <flens/lapack/la/larzt.tcc>
#include <flens/lapack/la/lascl.tcc>
#include <flens/lapack/la/lassq.tcc>
#include <flens/lapack/la/laswp.tcc>
#include <flens/lapack/la/lasy2.tcc>
#include <flens/lapack/la/latrs.tcc>
#include <flens/lapack/la/latrz.tcc>
#include <flens/lapack/la/lauu2.tcc>
#include <flens/lapack/la/lauum.tcc>

#include <flens/lapack/pb/pbsv.tcc>
#include <flens/lapack/pb/pbtrf.tcc>
#include <flens/lapack/pb/pbtrs.tcc>

#include <flens/lapack/po/pocon.tcc>
#include <flens/lapack/po/posv.tcc>
#include <flens/lapack/po/potf2.tcc>
#include <flens/lapack/po/potrf.tcc>
#include <flens/lapack/po/potri.tcc>
#include <flens/lapack/po/potrs.tcc>

#include <flens/lapack/sb/ev.tcc>

#include <flens/lapack/sp/ev.tcc>
#include <flens/lapack/sp/sv.tcc>
#include <flens/lapack/sp/trf.tcc>
#include <flens/lapack/sp/tri.tcc>
#include <flens/lapack/sp/trs.tcc>

#include <flens/lapack/sy/ev.tcc>
#include <flens/lapack/sy/sv.tcc>
#include <flens/lapack/sy/trf.tcc>
#include <flens/lapack/sy/tri.tcc>
#include <flens/lapack/sy/trs.tcc>

#include <flens/lapack/tb/trs.tcc>

#include <flens/lapack/tp/tri.tcc>
#include <flens/lapack/tp/trs.tcc>

#include <flens/lapack/tr/ti2.tcc>
#include <flens/lapack/tr/tri.tcc>
#include <flens/lapack/tr/trs.tcc>

#endif // FLENS_LAPACK_LAPACK_TCC
