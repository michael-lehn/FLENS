/*
 *   Copyright (c) 2014, Michael Lehn
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

/* Based on
 *
       SUBROUTINE DLASRT( ID, N, D, INFO )
 *
 *  -- LAPACK routine (version 3.2) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *     November 2006
 */

#ifndef FLENS_LAPACK_LA_LASRT_TCC
#define FLENS_LAPACK_LA_LASRT_TCC 1

#include <flens/auxiliary/auxiliary.h>
#include <flens/lapack/typedefs.h>
#include <flens/matrixtypes/matrixtypes.h>
#include <flens/vectortypes/vectortypes.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

namespace generic {

template <typename VD>
void
lasrt_impl(bool              increasing,
           DenseVector<VD>   &d)
{
    using std::swap;

    typedef typename VD::ElementType  T;
    typedef typename VD::IndexType    IndexType;

    const IndexType   select = 20;
    const IndexType   n      = d.length();

    IndexType   stack_data[2*32];
    IndexType   *stack1 = &stack_data[0] - 1;
    IndexType   *stack2 = &stack_data[32] - 1;

    IndexType   stackPointer;

    stackPointer = 1;
    stack1[1] = 1;
    stack2[1] = n;

    do {
        const IndexType start = stack1[stackPointer];
        const IndexType end   = stack2[stackPointer];
        --stackPointer;

        if (end-start<=select && end-start>0) {
//
//          Do Insertion sort on D( START:ENDD )
//
            if (! increasing) {
//
//              Sort into decreasing order
//
                for (IndexType i=start+1; i<=end; ++i) {
                    for (IndexType j=i; j>=start+1; --j) {
                        if (d(j)>d(j-1)) {
                            swap(d(j-1), d(j));
                        } else {
                            break;
                        }
                    }
                }

            } else {
//
//              Sort into increasing order
//
                for (IndexType i=start+1; i<=end; ++i) {
                    for (IndexType j=i; j>=start+1; --j) {
                        if (d(j)<d(j-1)) {
                            swap(d(j-1),d(j));
                        } else {
                            break;
                        }
                    }
                }

            }
//
        } else if (end-start>select) {
//
//          Partition D( START:ENDD ) and stack parts, largest one first
//
//          Choose partition entry as median of 3
//
            IndexType i  = (start+end)/2;

            T  d1 = d(start);
            T  d2 = d(end);
            T  d3 = d(i);

            T dx;
            if (d1<d2) {
                if (d3<d1) {
                    dx = d1;
                } else if (d3<d2) {
                    dx = d3;
                } else {
                    dx = d2;
                }
            } else {
                if (d3<d2) {
                    dx = d2;
                } else if (d3<d1) {
                    dx = d3;
                } else {
                    dx = d1;
                }
            }

            if (! increasing) {
//
//              Sort into decreasing order
//
                IndexType j;

                i = start - 1;
                j = end + 1;

                while (true) {

                    do {
                        --j;
                    } while (d(j)<dx);

                    do {
                        ++i;
                    } while (d(i)>dx);

                    if (i<j) {
                        swap(d(i), d(j));
                    } else {
                        break;
                    }
                }

                if (j-start>end-j-1) {
                    ++stackPointer;
                    ASSERT(stackPointer<=32);
                    stack1[stackPointer] = start;
                    stack2[stackPointer] = j;

                    ++stackPointer;
                    ASSERT(stackPointer<=32);
                    stack1[stackPointer] = j+1;
                    stack2[stackPointer] = end;
                } else {
                    ++stackPointer;
                    ASSERT(stackPointer<=32);
                    stack1[stackPointer] = j+1;
                    stack2[stackPointer] = end;

                    ++stackPointer;
                    ASSERT(stackPointer<=32);
                    stack1[stackPointer] = start;
                    stack2[stackPointer] = j;
                }
            } else {
//
//              Sort into increasing order
//
                IndexType j;

                i = start - 1;
                j = end + 1;

                while (true) {

                    do {
                        --j;
                    } while (d(j)>dx);

                    do {
                        ++i;
                    } while (d(i)<dx);

                    if (i<j) {
                        swap(d(i), d(j));
                    } else {
                        break;
                    }
                }
                if (j-start>end-j-1) {
                    ++stackPointer;
                    ASSERT(stackPointer<=32);
                    stack1[stackPointer] = start;
                    stack2[stackPointer] = j;

                    ++stackPointer;
                    ASSERT(stackPointer<=32);
                    stack1[stackPointer] = j+1;
                    stack2[stackPointer] = end;
                } else {
                    ++stackPointer;
                    ASSERT(stackPointer<=32);
                    stack1[stackPointer] = j+1;
                    stack2[stackPointer] = end;

                    ++stackPointer;
                    ASSERT(stackPointer<=32);
                    stack1[stackPointer] = start;
                    stack2[stackPointer] = j;
                }
            }
        }
    } while (stackPointer>0);
}

} // namespace generic

//== interface for native lapack ===============================================

#ifdef USE_CXXLAPACK

namespace external {

template <typename VD>
void
lasrt_impl(bool             increasing,
           DenseVector<VD>  &d)
{
    cxxlapack::lasrt(increasing ? 'I' : 'D',
                     d.length(),
                     d.data());
}

} // namespace external

#endif

//== public interface ==========================================================

template <typename VD>
typename RestrictTo<IsRealDenseVector<VD>::value,
         void>::Type
lasrt(bool increasing,
      VD   &&d)
{
    LAPACK_DEBUG_OUT("lasrt");
//
//  Remove references from rvalue types
//
#   ifdef CHECK_CXXLAPACK
    typedef typename RemoveRef<VD>::Type    VectorD;
#   endif

//
//  Test the input parameters
//
#   ifndef NDEBUG
    ASSERT(d.firstIndex()==1);
#   endif

//
//  Make copies of output arguments
//
#   ifdef CHECK_CXXLAPACK
    typename VectorD::NoView      d_org    = d;
#   endif

//
//  Call implementation
//
    LAPACK_SELECT::lasrt_impl(increasing, d);

#   ifdef CHECK_CXXLAPACK
//
//  Compare results
//
    typename VectorD::NoView      d_generic    = d;

    d = d_org;

    external::lasrt_impl(increasing, d);

    if (! isIdentical(d_generic, d, "d_generic", "d")) {
        std::cerr << "CXXLAPACK: d_generic = " << d_generic << std::endl;
        std::cerr << "F77LAPACK: d = " << d << std::endl;
        ASSERT(0);
    }
#   endif
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_LA_LASRT_TCC
