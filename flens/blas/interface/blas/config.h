#include <cmath>
#include <flens/flens.cxx>

#ifndef BLAS
#define BLAS(x)  x##_
#endif

#ifndef INTEGER
#define INTEGER int
#endif

extern "C" {

void
BLAS(xerbla)(const char *SRNAME, const INTEGER *INFO);

} // extern "C"


namespace flens {

template <typename ENUM>
typename RestrictTo<IsSame<ENUM,Transpose>::value,
                    Transpose>::Type
convertTo(const char c)
{
    if (c=='N') {
        return NoTrans;
    } else if (c=='T') {
        return Trans;
    } else if (c=='C') {
        return ConjTrans;
    } else if (c=='R') {
        return Conj;
    }
    std::cerr << "c = " << c << std::endl;
    ASSERT(0);
    return NoTrans;
}

typedef std::complex<float>                   cfloat;
typedef std::complex<double>                  cdouble;

//
//  Set index base and type
//
typedef IndexOptions<INTEGER, 1>              IndexBase;

//
//  Array types
//
typedef Array<float, IndexBase>               SArray;
typedef ArrayView<float, IndexBase>           SArrayView;
typedef ConstArrayView<float, IndexBase>      SConstArrayView;

typedef Array<double, IndexBase>              DArray;
typedef ArrayView<double, IndexBase>          DArrayView;
typedef ConstArrayView<double, IndexBase>     DConstArrayView;

typedef Array<cfloat, IndexBase>              CArray;
typedef ArrayView<cfloat, IndexBase>          CArrayView;
typedef ConstArrayView<cfloat, IndexBase>     CConstArrayView;

typedef Array<cdouble, IndexBase>             ZArray;
typedef ArrayView<cdouble, IndexBase>         ZArrayView;
typedef ConstArrayView<cdouble, IndexBase>    ZConstArrayView;

//
//  Vector types
//
typedef DenseVector<SArray>                   SDenseVector;
typedef DenseVector<SArrayView>               SDenseVectorView;
typedef DenseVector<SConstArrayView>          SDenseVectorConstView;

typedef DenseVector<DArray>                   DDenseVector;
typedef DenseVector<DArrayView>               DDenseVectorView;
typedef DenseVector<DConstArrayView>          DDenseVectorConstView;

typedef DenseVector<CArray>                   CDenseVector;
typedef DenseVector<CArrayView>               CDenseVectorView;
typedef DenseVector<CConstArrayView>          CDenseVectorConstView;

typedef DenseVector<ZArray>                   ZDenseVector;
typedef DenseVector<ZArrayView>               ZDenseVectorView;
typedef DenseVector<ZConstArrayView>          ZDenseVectorConstView;

//
//  FullStorage types
//
typedef FullStorage<float, ColMajor, IndexBase>             SFull;
typedef FullStorageView<float, ColMajor, IndexBase>         SFullView;
typedef ConstFullStorageView<float, ColMajor, IndexBase>    SFullConstView;

typedef FullStorage<double, ColMajor, IndexBase>            DFull;
typedef FullStorageView<double, ColMajor, IndexBase>        DFullView;
typedef ConstFullStorageView<double, ColMajor, IndexBase>   DFullConstView;

typedef FullStorage<cfloat, ColMajor, IndexBase>            CFull;
typedef FullStorageView<cfloat, ColMajor, IndexBase>        CFullView;
typedef ConstFullStorageView<cfloat, ColMajor, IndexBase>   CFullConstView;

typedef FullStorage<cdouble, ColMajor, IndexBase>           ZFull;
typedef FullStorageView<cdouble, ColMajor, IndexBase>       ZFullView;
typedef ConstFullStorageView<cdouble, ColMajor, IndexBase>  ZFullConstView;

//
//  GeMatrix types
//
typedef GeMatrix<SFull>             SGeMatrix;
typedef GeMatrix<SFullView>         SGeMatrixView;
typedef GeMatrix<SFullConstView>    SGeMatrixConstView;

typedef GeMatrix<DFull>             DGeMatrix;
typedef GeMatrix<DFullView>         DGeMatrixView;
typedef GeMatrix<DFullConstView>    DGeMatrixConstView;

typedef GeMatrix<CFull>             CGeMatrix;
typedef GeMatrix<CFullView>         CGeMatrixView;
typedef GeMatrix<CFullConstView>    CGeMatrixConstView;

typedef GeMatrix<ZFull>             ZGeMatrix;
typedef GeMatrix<ZFullView>         ZGeMatrixView;
typedef GeMatrix<ZFullConstView>    ZGeMatrixConstView;

//
//  HeMatrix types
//
typedef HeMatrix<CFull>             CHeMatrix;
typedef HeMatrix<CFullView>         CHeMatrixView;
typedef HeMatrix<CFullConstView>    CHeMatrixConstView;

typedef HeMatrix<ZFull>             ZHeMatrix;
typedef HeMatrix<ZFullView>         ZHeMatrixView;
typedef HeMatrix<ZFullConstView>    ZHeMatrixConstView;

//
//  SyMatrix types
//
typedef SyMatrix<SFull>             SSyMatrix;
typedef SyMatrix<SFullView>         SSyMatrixView;
typedef SyMatrix<SFullConstView>    SSyMatrixConstView;

typedef SyMatrix<DFull>             DSyMatrix;
typedef SyMatrix<DFullView>         DSyMatrixView;
typedef SyMatrix<DFullConstView>    DSyMatrixConstView;

typedef SyMatrix<CFull>             CSyMatrix;
typedef SyMatrix<CFullView>         CSyMatrixView;
typedef SyMatrix<CFullConstView>    CSyMatrixConstView;

typedef SyMatrix<ZFull>             ZSyMatrix;
typedef SyMatrix<ZFullView>         ZSyMatrixView;
typedef SyMatrix<ZFullConstView>    ZSyMatrixConstView;


//
//  TrMatrix types
//
typedef TrMatrix<SFull>             STrMatrix;
typedef TrMatrix<SFullView>         STrMatrixView;
typedef TrMatrix<SFullConstView>    STrMatrixConstView;

typedef TrMatrix<DFull>             DTrMatrix;
typedef TrMatrix<DFullView>         DTrMatrixView;
typedef TrMatrix<DFullConstView>    DTrMatrixConstView;

typedef TrMatrix<CFull>             CTrMatrix;
typedef TrMatrix<CFullView>         CTrMatrixView;
typedef TrMatrix<CFullConstView>    CTrMatrixConstView;

typedef TrMatrix<ZFull>             ZTrMatrix;
typedef TrMatrix<ZFullView>         ZTrMatrixView;
typedef TrMatrix<ZFullConstView>    ZTrMatrixConstView;



} // namespace flens


