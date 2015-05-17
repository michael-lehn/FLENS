#ifndef LU_TILE_H
#define LU_TILE_H 1

#include <flens/flens.h>

namespace flens {

template <typename MA, typename IndexType, typename VB>
void
copyToTiles(IndexType bs, const GeMatrix<MA> &A, DenseVector<VB> &b)
{
    typedef typename GeMatrix<MA>::View   TileView;

    IndexType m = A.numRows();
    IndexType n = A.numCols();

    Underscore<IndexType> _;

    IndexType ib = 0;

    for (IndexType j=1; j<=n; j+=bs) {
        IndexType nb = std::min(bs, n-j+1);

        for (IndexType i=1; i<=m; i+=bs) {
            IndexType mb = std::min(bs, m-i+1);
            TileView  B  = TileView(mb, nb, b(_(ib+1,ib+bs*bs)), bs);

            B = A(_(i,i+mb-1),_(j,j+nb-1));
            ib += bs*bs;
        }
    }
}

template <typename MA, typename IndexType, typename VB>
void
copyFromTiles(IndexType bs, const DenseVector<VB> &b, GeMatrix<MA> &A)
{
    typedef typename GeMatrix<MA>::ConstView   TileView;

    IndexType m = A.numRows();
    IndexType n = A.numCols();

    Underscore<IndexType> _;

    IndexType ib = 0;

    for (IndexType j=1; j<=n; j+=bs) {
        IndexType nb = std::min(bs, n-j+1);

        for (IndexType i=1; i<=m; i+=bs) {
            IndexType mb = std::min(bs, m-i+1);
            TileView  B  = TileView(mb, nb, b(_(ib+1,ib+bs*bs)), bs);

            A(_(i,i+mb-1),_(j,j+nb-1)) = B;
            ib += bs*bs;
        }
    }
}

template <typename GeMatrixType>
class TiledCopy
    : public GeneralMatrix<GeMatrixType>
{
    public:

        typedef typename GeMatrixType::View         TileView;
        typedef typename GeMatrixType::ElementType  ElementType;
        typedef typename GeMatrixType::IndexType    IndexType;

        TiledCopy(const GeMatrixType &A, int bs_)
            : numRows(A.numRows()), numCols(A.numCols()),
              bs(bs_), mr(numRows % bs), nr(numCols % bs),
              tileRows((numRows+bs-1) / bs),
              tileCols((numCols+bs-1) / bs)
        {
            tiles.resize(bs*bs*tileRows*tileCols);
            copyToTiles(bs, A, tiles);
        }

        int
        numTileRows() const
        {
            return tileRows;
        }

        int
        numTileCols() const
        {
            return tileCols;
        }

        int
        blockSize() const
        {
            return bs;
        }

        TileView
        operator()(int i, int j)
        {
            int mb = (i<tileRows || mr==0) ? bs : mr;
            int nb = (j<tileCols || nr==0) ? bs : nr;
            int ib = tileRows*bs*bs*(j-1) + bs*bs*(i-1);

            return TileView(mb, nb, tiles(Range<int>(ib+1,ib+bs*bs)), bs);
        }

        void
        untile(GeMatrixType &A) const
        {
            copyFromTiles(bs, tiles, A);
        }

    private:

        int                             numRows, numCols;
        int                             bs, mr, nr;
        int                             tileRows, tileCols;
        typename GeMatrixType::Vector   tiles;
};

} // namespace flens

#endif // LU_TILE_H
