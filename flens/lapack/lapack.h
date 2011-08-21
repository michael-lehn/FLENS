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

#include <flens/lapack/aux/ilaenv.h>
#include <flens/lapack/aux/ilalc.h>
#include <flens/lapack/aux/ilalr.h>
#include <flens/lapack/aux/labad.h>
#include <flens/lapack/aux/lamch.h>
#include <flens/lapack/aux/lanv2.h>
#include <flens/lapack/aux/lapy2.h>
#include <flens/lapack/aux/larf.h>
#include <flens/lapack/aux/larfb.h>
#include <flens/lapack/aux/larfg.h>
#include <flens/lapack/aux/larft.h>
#include <flens/lapack/aux/laswp.h>
#include <flens/lapack/aux/sign.h>

#include <flens/lapack/eig/hd2.h>

#include <flens/lapack/gesv/sv.h>
#include <flens/lapack/gesv/tf2.h>
#include <flens/lapack/gesv/trf.h>
#include <flens/lapack/gesv/trs.h>

#include <flens/lapack/qr/org2r.h>
#include <flens/lapack/qr/orgqr.h>
#include <flens/lapack/qr/orm2r.h>
#include <flens/lapack/qr/ormqr.h>
#include <flens/lapack/qr/qr2.h>
#include <flens/lapack/qr/qrf.h>
#include <flens/lapack/qr/qrs.h>

#include <flens/lapack/typedefs.h>

#endif // FLENS_LAPACK_LAPACK_H
