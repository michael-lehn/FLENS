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

#ifndef FLENS_STORAGE_TINYFULLSTORAGE_TINYFULLSTORAGEVIEW_H
#define FLENS_STORAGE_TINYFULLSTORAGE_TINYFULLSTORAGEVIEW_H 1

namespace flens {
    
    template <typename T, int m, int n, int ldA, int indexBase=1>
    class TinyFullStorageView
    {
    public:
        typedef T    ElementType;
        typedef int  IndexType;
        
        static const int numRows = m;
        static const int numCols = n;
        static const int leadingDimension = ldA;
        static const int firstRow = indexBase;
        static const int firstCol = indexBase;
        static const int lastRow = firstRow+numRows-1;
        static const int lastCol = firstCol+numCols-1;
        
        TinyFullStorageView(ElementType *data);
        
        ~TinyFullStorageView();
        
        //-- operators ---------------------------------------------------------
        
        const ElementType &
        operator()(IndexType row, IndexType col) const;
        
        ElementType &
        operator()(IndexType row, IndexType col);
        
        //-- methods -----------------------------------------------------------
        
        const ElementType *
        data() const;
        
        ElementType *
        data();
        
        void
        fill(const ElementType &value = ElementType(0));
        
    private:
        ElementType *_data;
    };
    
} // namespace flens

#endif // FLENS_STORAGE_TINYFULLSTORAGE_TINYFULLSTORAGEVIEW_H
