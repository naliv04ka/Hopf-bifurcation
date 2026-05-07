#ifndef GRID_H
#define GRID_H

#include "common.h"

class Grid {
public:
    Grid(Index nx, Index ny, Real Lx, Real Ly);
    
    // Grid dimensions
    Index getNx() const { return nx_; }
    Index getNy() const { return ny_; }
    Real getLx() const { return Lx_; }
    Real getLy() const { return Ly_; }
    
    // Grid spacing
    Real getDx() const { return dx_; }
    Real getDy() const { return dy_; }
    
    // Grid coordinates
    Real getX(Index i) const { return i * dx_; }
    Real getY(Index j) const { return j * dy_; }
    
    // Total number of cells
    Index getNumCells() const { return nx_ * ny_; }
    
    // Index conversion
    Index getIndex(Index i, Index j) const { return j * nx_ + i; }
    void getIJ(Index idx, Index& i, Index& j) const {
        i = idx % nx_;
        j = idx / nx_;
    }
    
private:
    Index nx_, ny_;
    Real Lx_, Ly_;
    Real dx_, dy_;
};

#endif // GRID_H
