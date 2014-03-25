/*
 *   Copyright (c) 2012, Michael Lehn, Klaus Pototzky
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

#ifndef FLENS_MATRIXTYPES_GENERAL_IMPL_DIAGMATRIX_H
#define FLENS_MATRIXTYPES_GENERAL_IMPL_DIAGMATRIX_H 1

#include <cxxblas/typedefs.h>
#include <flens/auxiliary/auxiliary.h>
#include <flens/typedefs.h>

namespace flens {
    
    // forward declarations
    template <typename A>
    class DenseVector;
    
    template <typename FS>
    class GbMatrix;
    
    template <typename FS>
    class HbMatrix;
    
    template <typename FS>
    class SbMatrix;
    
    template <typename FS>
    class TbMatrix;
    
    template <typename FS>
    class DiagMatrix;
    
    template <typename FS>
    class DiagMatrix
    : public GeneralMatrix<DiagMatrix<FS> >
    {
    public:
        typedef FS                                  Engine;
        typedef typename Engine::ElementType        ElementType;
        typedef typename Engine::IndexType          IndexType;
        
        // view types from Engine
        typedef typename Engine::ConstView          EngineConstView;
        typedef typename Engine::View               EngineView;
        typedef typename Engine::NoView             EngineNoView;
        
        // view types for
        typedef DenseVector<EngineConstView>        ConstVectorView;
        typedef DenseVector<EngineView>             VectorView;
        typedef DenseVector<EngineNoView>           Vector;
        
        typedef DiagMatrix<EngineConstView>         ConstView;
        typedef DiagMatrix<EngineView>              View;
        typedef DiagMatrix<EngineNoView>            NoView;
        
        DiagMatrix();
        
        DiagMatrix(IndexType dim,
                IndexType firstIndex = Engine::defaultIndexBase);
        
        DiagMatrix(const Engine &engine);
        
        DiagMatrix(const DiagMatrix &rhs);
        
        template <typename RHS>
        DiagMatrix(const DiagMatrix<RHS> &rhs);
        
        template <typename RHS>
        DiagMatrix(const Matrix<RHS> &rhs);
        
        template <typename RHS>
        DiagMatrix(DiagMatrix<RHS> &rhs);
        
        template <typename RHS>
        DiagMatrix(const DenseVector<RHS> &rhs);
        
        // -- operators --------------------------------------------------------
        
        DiagMatrix &
        operator=(const DiagMatrix &rhs);
        
        template <typename RHS>
        DiagMatrix &
        operator=(const Matrix<RHS> &rhs);
        
        template <typename RHS>
        DiagMatrix &
        operator=(const DenseVector<RHS> &rhs);
        
        template <typename RHS>
        DiagMatrix &
        operator+=(const Matrix<RHS> &rhs);
        
        template <typename RHS>
        DiagMatrix &
        operator-=(const Matrix<RHS> &rhs);
        
        DiagMatrix<FS> &
        operator=(const ElementType &alpha);
        
        DiagMatrix<FS> &
        operator+=(const ElementType &alpha);
        
        DiagMatrix<FS> &
        operator-=(const ElementType &alpha);
        
        DiagMatrix &
        operator*=(const ElementType &alpha);
        
        DiagMatrix &
        operator/=(const ElementType &alpha);
        
        const ElementType &
        operator()(IndexType row, IndexType col) const;
        
        ElementType &
        operator()(IndexType row, IndexType col);
        
        // -- views ------------------------------------------------------------

        // diag views
        VectorView
        diag();
        
        const ConstVectorView
        diag() const;
        
        // -- methods ----------------------------------------------------------
        
        IndexType
        dim() const;
        
        IndexType
        numCols() const;
        
        IndexType
        numRows() const;
        
        IndexType
        firstIndex() const;
        
        IndexType
        lastIndex() const;
        
        IndexType
        leadingDimension() const;
        
        const ElementType *
        data() const;
        
        ElementType *
        data();
        
        template <typename RHS>
        bool
        resize(const DiagMatrix<RHS> &rhs,
               const ElementType &value = ElementType());
        
        bool
        resize(IndexType dim,
               IndexType firstIndex = Engine::defaultIndexBase,
               const ElementType &value = ElementType());
        
        bool
        fill(const ElementType &value = ElementType(0));
        
        bool
        fillRandom();
        
        // -- implementation ---------------------------------------------------
        
        const Engine &
        engine() const;
        
        Engine &
        engine();
        
    private:
        Engine       _engine;
        StorageUpLo  _upLo;
    };
    
    //-- Traits --------------------------------------------------------------------
    //
    //  IsDiagMatrix
    //
    struct _DiagMatrixChecker
    {
        
        struct Two {
            char x;
            char y;
        };
        
        static Two
        check(_AnyConversion);
        
        template <typename Any>
        static char
        check(DiagMatrix<Any>);
    };
    
    template <typename T>
    struct IsDiagMatrix
    {
        static T var;
        static const bool value = sizeof(_DiagMatrixChecker::check(var))==1;
    };
    
    //
    //  IsRealDiagMatrix
    //
    template <typename T>
    struct IsRealDiagMatrix
    {
        typedef typename std::remove_reference<T>::type  TT;
        
        static const bool value = IsDiagMatrix<TT>::value
        && IsNotComplex<typename TT::ElementType>::value;
    };
    
    //
    //  IsComplexDiagMatrix
    //
    template <typename T>
    struct IsComplexDiagMatrix
    {
        typedef typename std::remove_reference<T>::type  TT;
        
        static const bool value = IsDiagMatrix<TT>::value
        && IsComplex<typename TT::ElementType>::value;
    };
    
    
} // namespace flens

#endif // FLENS_MATRIXTYPES_GENERAL_IMPL_DIAGMATRIX_H
