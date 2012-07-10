/*
 *   Copyright (c) 2012, Michael Lehn
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

#ifndef CXXLAPACK_INTERFACE_INTERFACE_TCC
#define CXXLAPACK_INTERFACE_INTERFACE_TCC 1

#include <cxxlapack/interface/gbsv.tcc>
#include <cxxlapack/interface/gbtrf.tcc>
#include <cxxlapack/interface/gbtrs.tcc>
#include <cxxlapack/interface/gebak.tcc>
#include <cxxlapack/interface/gebal.tcc>
#include <cxxlapack/interface/gecon.tcc>
#include <cxxlapack/interface/geequ.tcc>
#include <cxxlapack/interface/gees.tcc>
#include <cxxlapack/interface/geesx.tcc>
#include <cxxlapack/interface/geev.tcc>
#include <cxxlapack/interface/geevx.tcc>
#include <cxxlapack/interface/gehd2.tcc>
#include <cxxlapack/interface/gehrd.tcc>
#include <cxxlapack/interface/gejsv.tcc>
#include <cxxlapack/interface/gelq2.tcc>
#include <cxxlapack/interface/gelqf.tcc>
#include <cxxlapack/interface/gels.tcc>
#include <cxxlapack/interface/gelss.tcc>
#include <cxxlapack/interface/gelsy.tcc>
#include <cxxlapack/interface/geqp3.tcc>
#include <cxxlapack/interface/geqr2.tcc>
#include <cxxlapack/interface/geqrf.tcc>
#include <cxxlapack/interface/gerfs.tcc>
#include <cxxlapack/interface/gesdd.tcc>
#include <cxxlapack/interface/gesv.tcc>
#include <cxxlapack/interface/gesvd.tcc>
#include <cxxlapack/interface/gesvj.tcc>
#include <cxxlapack/interface/gesvx.tcc>
#include <cxxlapack/interface/getf2.tcc>
#include <cxxlapack/interface/getrf.tcc>
#include <cxxlapack/interface/getri.tcc>
#include <cxxlapack/interface/getrs.tcc>
#include <cxxlapack/interface/gsvj0.tcc>
#include <cxxlapack/interface/gsvj1.tcc>
#include <cxxlapack/interface/hbev.tcc>
#include <cxxlapack/interface/heev.tcc>
#include <cxxlapack/interface/hpev.tcc>
#include <cxxlapack/interface/hseqr.tcc>
#include <cxxlapack/interface/ilalc.tcc>
#include <cxxlapack/interface/ilalr.tcc>
#include <cxxlapack/interface/labad.tcc>
#include <cxxlapack/interface/lacn2.tcc>
#include <cxxlapack/interface/laexc.tcc>
#include <cxxlapack/interface/lahqr.tcc>
#include <cxxlapack/interface/lahr2.tcc>
#include <cxxlapack/interface/laic1.tcc>
#include <cxxlapack/interface/laln2.tcc>
#include <cxxlapack/interface/lamch.tcc>
#include <cxxlapack/interface/lange.tcc>
#include <cxxlapack/interface/lantr.tcc>
#include <cxxlapack/interface/lanv2.tcc>
#include <cxxlapack/interface/lapy2.tcc>
#include <cxxlapack/interface/lapy3.tcc>
#include <cxxlapack/interface/laqge.tcc>
#include <cxxlapack/interface/laqp2.tcc>
#include <cxxlapack/interface/laqps.tcc>
#include <cxxlapack/interface/laqr0.tcc>
#include <cxxlapack/interface/laqr1.tcc>
#include <cxxlapack/interface/laqr2.tcc>
#include <cxxlapack/interface/laqr3.tcc>
#include <cxxlapack/interface/laqr4.tcc>
#include <cxxlapack/interface/laqr5.tcc>
#include <cxxlapack/interface/laqtr.tcc>
#include <cxxlapack/interface/larf.tcc>
#include <cxxlapack/interface/larfb.tcc>
#include <cxxlapack/interface/larfg.tcc>
#include <cxxlapack/interface/larft.tcc>
#include <cxxlapack/interface/larfx.tcc>
#include <cxxlapack/interface/lartg.tcc>
#include <cxxlapack/interface/larz.tcc>
#include <cxxlapack/interface/larzb.tcc>
#include <cxxlapack/interface/larzt.tcc>
#include <cxxlapack/interface/lascl.tcc>
#include <cxxlapack/interface/lassq.tcc>
#include <cxxlapack/interface/laswp.tcc>
#include <cxxlapack/interface/lasy2.tcc>
#include <cxxlapack/interface/latrs.tcc>
#include <cxxlapack/interface/latrz.tcc>
#include <cxxlapack/interface/lauu2.tcc>
#include <cxxlapack/interface/lauum.tcc>
#include <cxxlapack/interface/org2r.tcc>
#include <cxxlapack/interface/orghr.tcc>
#include <cxxlapack/interface/orgl2.tcc>
#include <cxxlapack/interface/orglq.tcc>
#include <cxxlapack/interface/orgqr.tcc>
#include <cxxlapack/interface/orm2r.tcc>
#include <cxxlapack/interface/ormhr.tcc>
#include <cxxlapack/interface/orml2.tcc>
#include <cxxlapack/interface/ormlq.tcc>
#include <cxxlapack/interface/ormqr.tcc>
#include <cxxlapack/interface/ormr3.tcc>
#include <cxxlapack/interface/ormrz.tcc>
#include <cxxlapack/interface/pocon.tcc>
#include <cxxlapack/interface/posv.tcc>
#include <cxxlapack/interface/potf2.tcc>
#include <cxxlapack/interface/potrf.tcc>
#include <cxxlapack/interface/potri.tcc>
#include <cxxlapack/interface/potrs.tcc>
#include <cxxlapack/interface/rscl.tcc>
#include <cxxlapack/interface/sbev.tcc>
#include <cxxlapack/interface/syev.tcc>
#include <cxxlapack/interface/trevc.tcc>
#include <cxxlapack/interface/trexc.tcc>
#include <cxxlapack/interface/trsen.tcc>
#include <cxxlapack/interface/trsna.tcc>
#include <cxxlapack/interface/trsyl.tcc>
#include <cxxlapack/interface/trti2.tcc>
#include <cxxlapack/interface/trtri.tcc>
#include <cxxlapack/interface/trtrs.tcc>
#include <cxxlapack/interface/tzrzf.tcc>
#include <cxxlapack/interface/unglq.tcc>
#include <cxxlapack/interface/ungqr.tcc>
#include <cxxlapack/interface/unmlq.tcc>
#include <cxxlapack/interface/unmqr.tcc>
#include <cxxlapack/interface/unmr3.tcc>
#include <cxxlapack/interface/unmrz.tcc>

#endif // CXXLAPACK_INTERFACE_INTERFACE_TCC
