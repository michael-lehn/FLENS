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

/* Based on
 *
       SUBROUTINE DTREXC( COMPQ, N, T, LDT, Q, LDQ, IFST, ILST, WORK,
      $                   INFO )
 *
 *  -- LAPACK routine (version 3.2) --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *     November 2006
 */

#ifndef FLENS_LAPACK_EIG_TREXC_TCC
#define FLENS_LAPACK_EIG_TREXC_TCC 1

#include <flens/blas/blas.h>
#include <flens/lapack/lapack.h>

namespace flens { namespace lapack {

//== generic lapack implementation =============================================

template <typename MT, typename MQ, typename IndexType, typename VWORK>
IndexType
trexc_generic(bool                          computeQ,
              GeMatrix<MT>                  &T,
              GeMatrix<MQ>                  &Q,
              IndexType                     &iFirst,
              IndexType                     &iLast,
              DenseVector<VWORK>            &work)
{
    typedef typename GeMatrix<MT>::ElementType  ElementType;

    const IndexType     n = T.numRows();
    const ElementType   Zero(0);

    IndexType nBlockFirst, nBlockLast, nBlockNext, here;
    IndexType info = 0;
//
//  Quick return if possible
//
    if (n<=1) {
        return info;
    }
//
//  Determine the first row of specified block
//  and find out if it is 1 by 1 or 2 by 2.
//
    if (iFirst>1) {
        if (T(iFirst,iFirst-1)!=Zero) {
            --iFirst;
        }
    }
    nBlockFirst = 1;
    if (iFirst<n) {
        if (T(iFirst+1,iFirst)!=Zero) {
            nBlockFirst = 2;
        }
    }
//
//  Determine the first row of the final block
//  and find out if it is 1 by 1 or 2 by 2.
//
    if (iLast>1) {
        if (T(iLast,iLast-1)!=Zero) {
            --iLast;
        }
    }
    nBlockLast = 1;
    if (iLast<n) {
        if (T(iLast+1,iLast)!=Zero) {
            nBlockLast = 2;
        }
    }

    if (iFirst==iLast) {
        return info;
    }

    if (iFirst<iLast) {
//
//      Update ILST
//
        if (nBlockFirst==2 && nBlockLast==1) {
            --iLast;
        }
        if (nBlockFirst==1 && nBlockLast==2) {
            ++iLast;
        }

        here = iFirst;

        do {
//
//          Swap block with next one below
//
            if (nBlockFirst==1 || nBlockFirst==2) {
//
//              Current block either 1 by 1 or 2 by 2
//
                nBlockNext = 1;
                if (here+nBlockFirst+1<=n) {
                    if (T(here+nBlockFirst+1,here+nBlockFirst)!=Zero) {
                        nBlockNext = 2;
                    }
                }
                info = laexc(computeQ, T, Q,
                             here, nBlockFirst, nBlockNext,
                             work);
                if (info!=0) {
                    iLast = here;
                    return info;
                }
                here += nBlockNext;
//
//              Test if 2 by 2 block breaks into two 1 by 1 blocks
//
                if (nBlockFirst==2) {
                    if (T(here+1,here)==Zero) {
                        nBlockFirst = 3;
                    }
                }

            } else {
//
//              Current block consists of two 1 by 1 blocks each of which
//              must be swapped individually
//
                nBlockNext = 1;
                if (here+3<=n) {
                    if (T(here+3,here+2)!=Zero) {
                        nBlockNext = 2;
                    }
                }
                info = laexc(computeQ, T, Q, here+1,
                             IndexType(1), nBlockNext,
                             work);
                if (info!=0) {
                    iLast = here;
                    return info;
                }
                if (nBlockNext==1) {
//
//                  Swap two 1 by 1 blocks, no problems possible
//
                    laexc(computeQ, T, Q, here, IndexType(1), nBlockNext, work);
                    ++here;
                } else {
//
//                  Recompute NBNEXT in case 2 by 2 split
//
                    if (T(here+2,here+1)==Zero) {
                        nBlockNext = 1;
                    }
                    if (nBlockNext==2) {
//
//                      2 by 2 Block did not split
//
                        info = laexc(computeQ, T, Q, here,
                                     IndexType(1), nBlockNext,
                                     work);
                        if (info!=0) {
                            iLast = here;
                            return info;
                        }
                        here += 2;
                    } else {
//
//                      2 by 2 Block did split
//
                        laexc(computeQ, T, Q,
                              here, IndexType(1), IndexType(1),
                              work);
                        laexc(computeQ, T, Q,
                              here+1, IndexType(1), IndexType(1),
                              work);
                        here += 2;
                    }
                }
            }
        } while (here<iLast);
    } else {

        here = iFirst;
        do {
//
//          Swap block with next one above
//
            if (nBlockFirst==1 || nBlockFirst==2) {
//
//              Current block either 1 by 1 or 2 by 2
//
                nBlockNext = 1;
                if (here>=3) {
                    if (T(here-1,here-2)!=Zero) {
                        nBlockNext = 2;
                    }
                }
                info = laexc(computeQ, T, Q,
                             here-nBlockNext, nBlockNext, nBlockFirst,
                             work);
                if (info!=0) {
                    iLast = here;
                    return info;
                }
                here -= nBlockNext;
//
//              Test if 2 by 2 block breaks into two 1 by 1 blocks
//
                if (nBlockFirst==2) {
                    if (T(here+1,here)==Zero) {
                        nBlockFirst = 3;
                    }
                }

            } else {
//
//              Current block consists of two 1 by 1 blocks each of which
//              must be swapped individually
//
                nBlockNext = 1;
                if (here>=3) {
                    if (T(here-1,here-2)!=Zero) {
                        nBlockNext = 2;
                    }
                }
                info = laexc(computeQ, T, Q,
                             here-nBlockNext, nBlockNext, IndexType(1),
                             work);
                if (info!=0) {
                    iLast = here;
                    return info;
                } 
                if (nBlockNext==1) {
//
//                  Swap two 1 by 1 blocks, no problems possible
//
                    laexc(computeQ, T, Q,
                          here, nBlockNext, IndexType(1),
                          work);
                    --here;
                } else {
//
//                  Recompute NBNEXT in case 2 by 2 split
//
                    if (T(here,here-1)==Zero) {
                        nBlockNext = 1;
                    }
                    if (nBlockNext==2) {
//
//                      2 by 2 Block did not split
//
                        info = laexc(computeQ, T, Q,
                                     here-1, IndexType(2), IndexType(1),
                                     work);
                        if (info!=0) {
                            iLast = here;
                            return info;
                        }
                        here -= 2;
                    } else {
//
//                      2 by 2 Block did split
//
                        laexc(computeQ, T, Q,
                              here, IndexType(1), IndexType(1),
                              work);
                        laexc(computeQ, T, Q,
                              here-1, IndexType(1), IndexType(1),
                              work);
                        here -= 2;
                    }
                }
            }
        } while (here>iLast);
    }
    iLast = here;
    return info;
}

//== interface for native lapack ===============================================

#ifdef TODO_CHECK_CXXLAPACK

template <typename MT, typename MQ, typename IndexType, typename VWORK>
IndexType
trexc_native(bool                          computeQ,
             GeMatrix<MT>                  &T,
             GeMatrix<MQ>                  &Q,
             IndexType                     &iFirst,
             IndexType                     &iLast,
             DenseVector<VWORK>            &work)
{
    typedef typename GeMatrix<MT>::ElementType ElementType;

    const char       COMPQ  = (computeQ) ? 'V' : 'N';
    const INTEGER    N      = T.numRows();
    const INTEGER    LDT    = T.leadingDimension();
    const INTEGER    LDQ    = Q.leadingDimension();
    INTEGER          IFST   = iFirst;
    INTEGER          ILST   = iLast;
    INTEGER          INFO;

    if (IsSame<ElementType,DOUBLE>::value) {
        LAPACK_IMPL(dtrexc)(&COMPQ,
                            &N,
                            T.data(),
                            &LDT,
                            Q.data(),
                            &LDQ,
                            &IFST,
                            &ILST,
                            work.data(),
                            &INFO);
    } else {
        ASSERT(0);
    }
    ASSERT(INFO>=0);

    iFirst  = IFST;
    iLast   = ILST;
    return INFO;
}

#endif // CHECK_CXXLAPACK

//== public interface ==========================================================

template <typename MT, typename MQ, typename IndexType, typename VWORK>
IndexType
trexc(bool                          computeQ,
      GeMatrix<MT>                  &T,
      GeMatrix<MQ>                  &Q,
      IndexType                     &iFirst,
      IndexType                     &iLast,
      DenseVector<VWORK>            &work)
{
    LAPACK_DEBUG_OUT("trexc");

//
//  Test the input parameters
//
#   ifndef NDEBUG
    ASSERT(T.firstRow()==1);
    ASSERT(T.firstCol()==1);
    ASSERT(T.numRows()==T.numCols());

    const IndexType n = T.numRows();

    if (computeQ) {
        ASSERT(Q.firstRow()==1);
        ASSERT(Q.firstCol()==1);
        ASSERT(Q.numRows()==Q.numCols());
        ASSERT(Q.numRows()==n);
    }

    ASSERT(work.firstIndex()==1);
    ASSERT(work.length()==n);
#   endif

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of output arguments
//
    typename GeMatrix<MT>::NoView           T_org      = T;
    typename GeMatrix<MQ>::NoView           Q_org      = Q;
    IndexType                               iFirst_org = iFirst;
    IndexType                               iLast_org  = iLast;
    typename DenseVector<VWORK>::NoView     work_org   = work;
#   endif

//
//  Call implementation
//
    IndexType info = trexc_generic(computeQ, T, Q, iFirst, iLast, work);

#   ifdef CHECK_CXXLAPACK
//
//  Make copies of results computed by the generic implementation
//
    typename GeMatrix<MT>::NoView           T_generic      = T;
    typename GeMatrix<MQ>::NoView           Q_generic      = Q;
    IndexType                               iFirst_generic = iFirst;
    IndexType                               iLast_generic  = iLast;
    typename DenseVector<VWORK>::NoView     work_generic   = work;

//
//  restore output arguments
//
    T       = T_org;
    Q       = Q_org;
    iFirst  = iFirst_org;
    iLast   = iLast_org;
    work    = work_org;

//
//  Compare results
//
    IndexType _info = trexc_native(computeQ, T, Q, iFirst, iLast, work);

    bool failed = false;
    if (! isIdentical(T_generic, T, "T_generic", "T")) {
        std::cerr << "CXXLAPACK: T_generic = " << T_generic << std::endl;
        std::cerr << "F77LAPACK: T = " << T << std::endl;
        failed = true;
    }

    if (! isIdentical(Q_generic, Q, "Q_generic", "Q")) {
        std::cerr << "CXXLAPACK: Q_generic = " << Q_generic << std::endl;
        std::cerr << "F77LAPACK: Q = " << Q << std::endl;
        failed = true;
    }

    if (! isIdentical(iFirst_generic, iFirst, "iFirst_generic", "iFirst")) {
        std::cerr << "CXXLAPACK:  iFirst_generic = "
                  << iFirst_generic << std::endl;
        std::cerr << "F77LAPACK: iFirst = "
                  << iFirst << std::endl;
        failed = true;
    }

    if (! isIdentical(iLast_generic, iLast, "iLast_generic", "iLast")) {
        std::cerr << "CXXLAPACK: iLast_generic = "
                  << iLast_generic << std::endl;
        std::cerr << "F77LAPACK: iLast = " << iLast << std::endl;
        failed = true;
    }

    if (! isIdentical(work_generic, work, "work_generic", "work")) {
        std::cerr << "CXXLAPACK: work_generic = " << work_generic << std::endl;
        std::cerr << "F77LAPACK: work = " << work << std::endl;
        failed = true;
    }

    if (! isIdentical(info, _info, " info", "_info")) {
        std::cerr << "CXXLAPACK:  info = " << info << std::endl;
        std::cerr << "F77LAPACK: _info = " << _info << std::endl;
        failed = true;
    }

    if (failed) {
        ASSERT(0);
    }
#   endif

    return info;
}

//-- forwarding ----------------------------------------------------------------
template <typename MT, typename MQ, typename IndexType, typename VWORK>
IndexType
trexc(bool                          computeQ,
      MT                            &&T,
      MQ                            &&Q,
      IndexType                     &iFirst,
      IndexType                     &iLast,
      VWORK                         &&work)
{
    CHECKPOINT_ENTER;
    const IndexType info = trexc(computeQ, T, Q, iFirst, iLast, work);
    CHECKPOINT_LEAVE;

    return info;
}

} } // namespace lapack, flens

#endif // FLENS_LAPACK_EIG_TREVC_TCC
