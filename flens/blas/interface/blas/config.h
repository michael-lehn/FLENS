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

//
//  BandStorage types
//
typedef BandStorage<float, ColMajor, IndexBase>             SBand;
typedef BandStorageView<float, ColMajor, IndexBase>         SBandView;
typedef ConstBandStorageView<float, ColMajor, IndexBase>    SBandConstView;

typedef BandStorage<double, ColMajor, IndexBase>            DBand;
typedef BandStorageView<double, ColMajor, IndexBase>        DBandView;
typedef ConstBandStorageView<double, ColMajor, IndexBase>   DBandConstView;

typedef BandStorage<cfloat, ColMajor, IndexBase>            CBand;
typedef BandStorageView<cfloat, ColMajor, IndexBase>        CBandView;
typedef ConstBandStorageView<cfloat, ColMajor, IndexBase>   CBandConstView;

typedef BandStorage<cdouble, ColMajor, IndexBase>           ZBand;
typedef BandStorageView<cdouble, ColMajor, IndexBase>       ZBandView;
typedef ConstBandStorageView<cdouble, ColMajor, IndexBase>  ZBandConstView;

//
//  GbMatrix types
//
typedef GbMatrix<SBand>             SGbMatrix;
typedef GbMatrix<SBandView>         SGbMatrixView;
typedef GbMatrix<SBandConstView>    SGbMatrixConstView;

typedef GbMatrix<DBand>             DGbMatrix;
typedef GbMatrix<DBandView>         DGbMatrixView;
typedef GbMatrix<DBandConstView>    DGbMatrixConstView;

typedef GbMatrix<CBand>             CGbMatrix;
typedef GbMatrix<CBandView>         CGbMatrixView;
typedef GbMatrix<CBandConstView>    CGbMatrixConstView;

typedef GbMatrix<ZBand>             ZGbMatrix;
typedef GbMatrix<ZBandView>         ZGbMatrixView;
typedef GbMatrix<ZBandConstView>    ZGbMatrixConstView;

//
//  HbMatrix types
//
typedef HbMatrix<CBand>             CHbMatrix;
typedef HbMatrix<CBandView>         CHbMatrixView;
typedef HbMatrix<CBandConstView>    CHbMatrixConstView;

typedef HbMatrix<ZBand>             ZHbMatrix;
typedef HbMatrix<ZBandView>         ZHbMatrixView;
typedef HbMatrix<ZBandConstView>    ZHbMatrixConstView;

//
//  SbMatrix types
//
typedef SbMatrix<SBand>             SSbMatrix;
typedef SbMatrix<SBandView>         SSbMatrixView;
typedef SbMatrix<SBandConstView>    SSbMatrixConstView;

typedef SbMatrix<DBand>             DSbMatrix;
typedef SbMatrix<DBandView>         DSbMatrixView;
typedef SbMatrix<DBandConstView>    DSbMatrixConstView;

typedef SbMatrix<CBand>             CSbMatrix;
typedef SbMatrix<CBandView>         CSbMatrixView;
typedef SbMatrix<CBandConstView>    CSbMatrixConstView;

typedef SbMatrix<ZBand>             ZSbMatrix;
typedef SbMatrix<ZBandView>         ZSbMatrixView;
typedef SbMatrix<ZBandConstView>    ZSbMatrixConstView;

//
//  TbMatrix types
//
typedef TbMatrix<SBand>             STbMatrix;
typedef TbMatrix<SBandView>         STbMatrixView;
typedef TbMatrix<SBandConstView>    STbMatrixConstView;

typedef TbMatrix<DBand>             DTbMatrix;
typedef TbMatrix<DBandView>         DTbMatrixView;
typedef TbMatrix<DBandConstView>    DTbMatrixConstView;

typedef TbMatrix<CBand>             CTbMatrix;
typedef TbMatrix<CBandView>         CTbMatrixView;
typedef TbMatrix<CBandConstView>    CTbMatrixConstView;

typedef TbMatrix<ZBand>             ZTbMatrix;
typedef TbMatrix<ZBandView>         ZTbMatrixView;
typedef TbMatrix<ZBandConstView>    ZTbMatrixConstView;

//
//  PackedStorage types
//
typedef PackedStorage<float, ColMajor, IndexBase>             SPacked;
typedef PackedStorageView<float, ColMajor, IndexBase>         SPackedView;
typedef ConstPackedStorageView<float, ColMajor, IndexBase>    SPackedConstView;

typedef PackedStorage<double, ColMajor, IndexBase>            DPacked;
typedef PackedStorageView<double, ColMajor, IndexBase>        DPackedView;
typedef ConstPackedStorageView<double, ColMajor, IndexBase>   DPackedConstView;

typedef PackedStorage<cfloat, ColMajor, IndexBase>            CPacked;
typedef PackedStorageView<cfloat, ColMajor, IndexBase>        CPackedView;
typedef ConstPackedStorageView<cfloat, ColMajor, IndexBase>   CPackedConstView;

typedef PackedStorage<cdouble, ColMajor, IndexBase>           ZPacked;
typedef PackedStorageView<cdouble, ColMajor, IndexBase>       ZPackedView;
typedef ConstPackedStorageView<cdouble, ColMajor, IndexBase>  ZPackedConstView;

//
//  HpMatrix types
//
typedef HpMatrix<CPacked>             CHpMatrix;
typedef HpMatrix<CPackedView>         CHpMatrixView;
typedef HpMatrix<CPackedConstView>    CHpMatrixConstView;

typedef HpMatrix<ZPacked>             ZHpMatrix;
typedef HpMatrix<ZPackedView>         ZHpMatrixView;
typedef HpMatrix<ZPackedConstView>    ZHpMatrixConstView;

//
//  SpMatrix typesPacked
//
typedef SpMatrix<SPacked>             SSpMatrix;
typedef SpMatrix<SPackedView>         SSpMatrixView;
typedef SpMatrix<SPackedConstView>    SSpMatrixConstView;

typedef SpMatrix<DPacked>             DSpMatrix;
typedef SpMatrix<DPackedView>         DSpMatrixView;
typedef SpMatrix<DPackedConstView>    DSpMatrixConstView;

typedef SpMatrix<CPacked>             CSpMatrix;
typedef SpMatrix<CPackedView>         CSpMatrixView;
typedef SpMatrix<CPackedConstView>    CSpMatrixConstView;

typedef SpMatrix<ZPacked>             ZSpMatrix;
typedef SpMatrix<ZPackedView>         ZSpMatrixView;
typedef SpMatrix<ZPackedConstView>    ZSpMatrixConstView;

//
//  TpMatrix typesPacked
//
typedef TpMatrix<SPacked>             STpMatrix;
typedef TpMatrix<SPackedView>         STpMatrixView;
typedef TpMatrix<SPackedConstView>    STpMatrixConstView;

typedef TpMatrix<DPacked>             DTpMatrix;
typedef TpMatrix<DPackedView>         DTpMatrixView;
typedef TpMatrix<DPackedConstView>    DTpMatrixConstView;

typedef TpMatrix<CPacked>             CTpMatrix;
typedef TpMatrix<CPackedView>         CTpMatrixView;
typedef TpMatrix<CPackedConstView>    CTpMatrixConstView;

typedef TpMatrix<ZPacked>             ZTpMatrix;
typedef TpMatrix<ZPackedView>         ZTpMatrixView;
typedef TpMatrix<ZPackedConstView>    ZTpMatrixConstView;







} // namespace flens


