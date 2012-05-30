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

#ifndef CXXLAPACK_INTERFACE_INTERFACE_H
#define CXXLAPACK_INTERFACE_INTERFACE_H 1

#include <cxxlapack/interface/gbsv.h>
#include <cxxlapack/interface/gbtrf.h>
#include <cxxlapack/interface/gbtrs.h>
#include <cxxlapack/interface/gebak.h>
#include <cxxlapack/interface/gebal.h>
#include <cxxlapack/interface/gecon.h>
#include <cxxlapack/interface/geequ.h>
#include <cxxlapack/interface/gees.h>
#include <cxxlapack/interface/geesx.h>
#include <cxxlapack/interface/geev.h>
#include <cxxlapack/interface/geevx.h>
#include <cxxlapack/interface/gehd2.h>
#include <cxxlapack/interface/gehrd.h>
#include <cxxlapack/interface/gejsv.h>
#include <cxxlapack/interface/gelq2.h>
#include <cxxlapack/interface/gelqf.h>
#include <cxxlapack/interface/gels.h>
#include <cxxlapack/interface/gelss.h>
#include <cxxlapack/interface/geqp3.h>
#include <cxxlapack/interface/geqr2.h>
#include <cxxlapack/interface/geqrf.h>
#include <cxxlapack/interface/gerfs.h>
#include <cxxlapack/interface/gesdd.h>
#include <cxxlapack/interface/gesv.h>
#include <cxxlapack/interface/gesvd.h>
#include <cxxlapack/interface/gesvj.h>
#include <cxxlapack/interface/gesvx.h>
#include <cxxlapack/interface/getf2.h>
#include <cxxlapack/interface/getrf.h>
#include <cxxlapack/interface/getri.h>
#include <cxxlapack/interface/getrs.h>
#include <cxxlapack/interface/gsvj0.h>
#include <cxxlapack/interface/gsvj1.h>
#include <cxxlapack/interface/hbev.h>
#include <cxxlapack/interface/heev.h>
#include <cxxlapack/interface/hpev.h>
#include <cxxlapack/interface/hseqr.h>
#include <cxxlapack/interface/ilalc.h>
#include <cxxlapack/interface/ilalr.h>
#include <cxxlapack/interface/labad.h>
#include <cxxlapack/interface/lacn2.h>
#include <cxxlapack/interface/laexc.h>
#include <cxxlapack/interface/lahqr.h>
#include <cxxlapack/interface/lahr2.h>
#include <cxxlapack/interface/laln2.h>
#include <cxxlapack/interface/lamch.h>
#include <cxxlapack/interface/lange.h>
#include <cxxlapack/interface/lantr.h>
#include <cxxlapack/interface/lanv2.h>
#include <cxxlapack/interface/lapy2.h>
#include <cxxlapack/interface/lapy3.h>
#include <cxxlapack/interface/laqge.h>
#include <cxxlapack/interface/laqp2.h>
#include <cxxlapack/interface/laqps.h>
#include <cxxlapack/interface/laqr0.h>
#include <cxxlapack/interface/laqr1.h>
#include <cxxlapack/interface/laqr2.h>
#include <cxxlapack/interface/laqr3.h>
#include <cxxlapack/interface/laqr4.h>
#include <cxxlapack/interface/laqr5.h>
#include <cxxlapack/interface/laqtr.h>
#include <cxxlapack/interface/larf.h>
#include <cxxlapack/interface/larfb.h>
#include <cxxlapack/interface/larfg.h>
#include <cxxlapack/interface/larft.h>
#include <cxxlapack/interface/larfx.h>
#include <cxxlapack/interface/lartg.h>
#include <cxxlapack/interface/lascl.h>
#include <cxxlapack/interface/lassq.h>
#include <cxxlapack/interface/laswp.h>
#include <cxxlapack/interface/lasy2.h>
#include <cxxlapack/interface/latrs.h>
#include <cxxlapack/interface/lauu2.h>
#include <cxxlapack/interface/lauum.h>
#include <cxxlapack/interface/org2r.h>
#include <cxxlapack/interface/orghr.h>
#include <cxxlapack/interface/orgl2.h>
#include <cxxlapack/interface/orglq.h>
#include <cxxlapack/interface/orgqr.h>
#include <cxxlapack/interface/orm2r.h>
#include <cxxlapack/interface/ormhr.h>
#include <cxxlapack/interface/orml2.h>
#include <cxxlapack/interface/ormlq.h>
#include <cxxlapack/interface/ormqr.h>
#include <cxxlapack/interface/pocon.h>
#include <cxxlapack/interface/posv.h>
#include <cxxlapack/interface/potf2.h>
#include <cxxlapack/interface/potrf.h>
#include <cxxlapack/interface/potri.h>
#include <cxxlapack/interface/potrs.h>
#include <cxxlapack/interface/rscl.h>
#include <cxxlapack/interface/sbev.h>
#include <cxxlapack/interface/syev.h>
#include <cxxlapack/interface/trevc.h>
#include <cxxlapack/interface/trexc.h>
#include <cxxlapack/interface/trsen.h>
#include <cxxlapack/interface/trsna.h>
#include <cxxlapack/interface/trsyl.h>
#include <cxxlapack/interface/trti2.h>
#include <cxxlapack/interface/trtri.h>
#include <cxxlapack/interface/trtrs.h>
#include <cxxlapack/interface/unglq.h>
#include <cxxlapack/interface/ungqr.h>
#include <cxxlapack/interface/unmlq.h>
#include <cxxlapack/interface/unmqr.h>

#endif // CXXLAPACK_INTERFACE_INTERFACE_H
