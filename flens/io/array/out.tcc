/*
 *   Copyright (c) 2010, Michael Lehn
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

#ifndef FLENS_IO_ARRAY_OUT_TCC
#define FLENS_IO_ARRAY_OUT_TCC 1

#include <flens/auxiliary/iscomplex.h>
#include <flens/io/array/out.h>

namespace flens {

template <typename A>
std::ostream &
operator<<(std::ostream &out, const DenseVector<A> &x)
{
    typedef typename DenseVector<A>::IndexType IndexType;
    typedef typename DenseVector<A>::ElementType ElementType;

#   ifdef FLENS_IO_WITH_RANGES
    IndexType defaultIndexBase = A::defaultIndexBase;

    out << std::endl << "[";
    if ((x.firstIndex()==defaultIndexBase) && (x.inc()>0)) {
        out << x.length();
    } else {
        out << x.firstIndex()
            << ".."
            << x.lastIndex();
    }
    out << "] ";
#   endif // FLENS_IO_WITH_RANGES

    out << std::endl;

    for (IndexType i=x.firstIndex(); i!=x.endIndex(); i+=x.inc()) {
        if (IsNotComplex<ElementType>::value)
                out.width(13);
            else
                out.width(28);

        out << x(i) << " ";
        if (i!=x.lastIndex()) {
            out << " ";
        }
    }
    out << std::endl;
    return out;
}

template <typename MA>
std::ostream &
operator<<(std::ostream &out, const DiagMatrix<MA> &A)
{
    typedef typename DiagMatrix<MA>::IndexType IndexType;
    typedef typename DiagMatrix<MA>::ElementType ElementType;
    
    const auto x = A.diag();

#   ifdef FLENS_IO_WITH_RANGES
    IndexType defaultIndexBase = A::defaultIndexBase;

    out << std::endl << "[";
    if ((x.firstIndex()==defaultIndexBase) && (x.inc()>0)) {
        out << x.length();
    } else {
        out << x.firstIndex()
            << ".."
            << x.lastIndex();
    }
    out << "] ";
#   endif // FLENS_IO_WITH_RANGES

    out << std::endl;

    for (IndexType i=x.firstIndex(); i!=x.endIndex(); i+=x.inc()) {
        if (IsNotComplex<ElementType>::value)
                out.width(13);
            else
                out.width(28);

        out << x(i) << " ";
        if (i!=x.lastIndex()) {
            out << " ";
        }
    }
    out << std::endl;
    return out;
}

} // namespace flens

#endif // FLENS_IO_ARRAY_OUT_TCC
