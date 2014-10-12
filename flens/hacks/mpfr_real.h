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

#include <cxxstd/limits.h>
#include <flens/auxiliary/explicit_cast.h>

namespace mpfr {

//
// auxiliary-functions used in numeric_limits
//

template <mpfr::real_prec_t prec_, mpfr::real_rnd_t rnd_>
    const mpfr::real<prec_,rnd_>
    nextabove(const mpfr::real<prec_,rnd_> &x);

template <mpfr::real_prec_t prec_, mpfr::real_rnd_t rnd_>
    const mpfr::real<prec_,rnd_>
    get_max();

} // namespace mpfr


namespace std {

//
// numeric_limits
// TODO: This just defines what gets used in FLENS-LAPACK get just copy
//       the rest from numeric_limits<double>
//

template <mpfr::real_prec_t prec_, mpfr::real_rnd_t rnd_>
class numeric_limits<mpfr::real<prec_,rnd_> >
    : public numeric_limits<double>
{
    public:
        typedef mpfr::real<prec_,rnd_>  T;

        static const T max_;
        static const T min_;
        static const T eps_;

        static const T
        epsilon();

        static const T
        max();

        static const T
        min();
};

//
// import is_floating_point to namespace std
//
template <mpfr::real_prec_t prec_, mpfr::real_rnd_t rnd_>
struct is_floating_point<mpfr::real<prec_,rnd_> >
    : public is_floating_point<double>
{
    public:

        static const bool value = true;
};

//
// import is_arithmetic to namespace std
//
template <mpfr::real_prec_t prec_, mpfr::real_rnd_t rnd_>
struct is_arithmetic<mpfr::real<prec_,rnd_> >
    : public is_arithmetic<double>
{
    public:

        static const bool value = true;
};

//
// import isnan to namespace std
//
template <mpfr::real_prec_t prec_, mpfr::real_rnd_t rnd_>
    bool
    isnan(const mpfr::real<prec_,rnd_> &x);

//
// import mpfr_real to namespace std
//
template <mpfr::real_prec_t prec_, mpfr::real_rnd_t rnd_>
    int
    signbit(const mpfr::real<prec_,rnd_> &x);

} // namespace std

#endif // defined(MPFR_REAL_HPP) && !defined(FLENS_HACKS_MPFR_REAL_H)
