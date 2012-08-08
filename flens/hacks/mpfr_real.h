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

#if defined(MPFR_REAL_HPP) && !defined(FLENS_HACKS_MPFR_REAL_H)
#define FLENS_HACKS_MPFR_REAL_H 1

#include <limits>
#include <flens/auxiliary/explicit_cast.h>

namespace mpfr {

//
// auxiliary-functions used in numeric_limits
//

template <mpfr::real_prec_t _prec, mpfr::real_rnd_t _rnd>
    const mpfr::real<_prec,_rnd>
    nextabove(const mpfr::real<_prec,_rnd> &x);

template <mpfr::real_prec_t _prec, mpfr::real_rnd_t _rnd>
    const mpfr::real<_prec,_rnd>
    get_max();

} // namespace mpfr


namespace std {

//
// numeric_limits
// TODO: This just defines what gets used in FLENS-LAPACK get just copy
//       the rest from numeric_limits<double>
//

template <mpfr::real_prec_t _prec, mpfr::real_rnd_t _rnd>
class numeric_limits<mpfr::real<_prec,_rnd> >
    : public numeric_limits<double>
{
    public:
        typedef mpfr::real<_prec,_rnd>  T;

        static const T _max;
        static const T _min;
        static const T _eps;

        static const T
        epsilon();

        static const T
        max();

        static const T
        min();
};

//
// import isnan to namespace std
//
template <mpfr::real_prec_t _prec, mpfr::real_rnd_t _rnd>
    bool
    isnan(const mpfr::real<_prec,_rnd> &x);

//
// import mpfr_real to namespace std
//
template <mpfr::real_prec_t _prec, mpfr::real_rnd_t _rnd>
    int
    signbit(const mpfr::real<_prec,_rnd> &x);

} // namespace std

#endif // defined(MPFR_REAL_HPP) && !defined(FLENS_HACKS_MPFR_REAL_H)
