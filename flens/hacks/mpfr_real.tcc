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

#if defined(MPFR_REAL_HPP) && !defined(FLENS_HACKS_MPFR_REAL_TCC)
#define FLENS_HACKS_MPFR_REAL_TCC 1

#include <external/real.hpp>
#include <flens/hacks/mpfr_real.h>

/*
 * NOTE: This hack requires that in the mpfr::real class the 'x_' attribute
 *       is made public!
 */

namespace mpfr {

//
// aux-functions used in numeric_limits
//

template <mpfr::real_prec_t prec_, mpfr::real_rnd_t rnd_>
const mpfr::real<prec_,rnd_>
nextabove(const mpfr::real<prec_,rnd_> &x)
{
    mpfr::real<prec_,rnd_> tmp = x;
    mpfr_nextabove(tmp.x_);
    return tmp;
}

template <mpfr::real_prec_t prec_, mpfr::real_rnd_t rnd_>
const mpfr::real<prec_,rnd_>
get_max()
{
    const unsigned long emax = mpfr_get_emax();
    mpfr::real<prec_,rnd_> tmp(1);

    mpfr_mul_2ui(tmp.x_, tmp.x_, emax-1, rnd_);
    return tmp;
}

} // namespace mpfr

//
// numeric_limits
//

namespace std {

template <mpfr::real_prec_t prec_, mpfr::real_rnd_t rnd_>
    const mpfr::real<prec_,rnd_>
    numeric_limits<mpfr::real<prec_,rnd_> >::max_
        = mpfr::get_max<prec_,rnd_>();

template <mpfr::real_prec_t prec_, mpfr::real_rnd_t rnd_>
    const mpfr::real<prec_,rnd_>
    numeric_limits<mpfr::real<prec_,rnd_> >::min_
        = mpfr::nextabove(T(0));

template <mpfr::real_prec_t prec_, mpfr::real_rnd_t rnd_>
    const mpfr::real<prec_,rnd_>
    numeric_limits<mpfr::real<prec_,rnd_> >::eps_
        = mpfr::nextabove(T(1)) - T(1);

template <mpfr::real_prec_t prec_, mpfr::real_rnd_t rnd_>
const mpfr::real<prec_,rnd_>
numeric_limits<mpfr::real<prec_,rnd_> >::epsilon()
{
     return eps_;
}

template <mpfr::real_prec_t prec_, mpfr::real_rnd_t rnd_>
const mpfr::real<prec_,rnd_>
numeric_limits<mpfr::real<prec_,rnd_> >::max()
{
    return max_;
}

template <mpfr::real_prec_t prec_, mpfr::real_rnd_t rnd_>
const mpfr::real<prec_,rnd_>
numeric_limits<mpfr::real<prec_,rnd_> >::min()
{
    return min_;
}

//
// import isnan to namespace std
//
template <mpfr::real_prec_t prec_, mpfr::real_rnd_t rnd_>
bool
isnan(const mpfr::real<prec_,rnd_> &x)
{
    return mpfr_nan_p(x.x_);
}

//
// import mpfr_real to namespace std
//
template <mpfr::real_prec_t prec_, mpfr::real_rnd_t rnd_>
int
signbit(const mpfr::real<prec_,rnd_> &x)
{
    return mpfr_signbit(x.x_);
}

} // namespace std

#endif // defined(MPFR_REAL_HPP) && !defined(FLENS_HACKS_MPFR_REAL_TCC)
