#include "field.h"

// ScalarField implementation
ScalarField::ScalarField(const Grid& grid, Real initValue)
    : grid_(grid), data_(grid.getNumCells(), initValue)
{
}

Real& ScalarField::operator()(Index i, Index j) {
    return data_[grid_.getIndex(i, j)];
}

const Real& ScalarField::operator()(Index i, Index j) const {
    return data_[grid_.getIndex(i, j)];
}

void ScalarField::fill(Real value) {
    std::fill(data_.begin(), data_.end(), value);
}

void ScalarField::copyFrom(const ScalarField& other) {
    data_ = other.data_;
}

void ScalarField::applyDirichletBC(Real value, bool left, bool right, bool bottom, bool top) {
    Index nx = grid_.getNx();
    Index ny = grid_.getNy();
    
    if (left) {
        for (Index j = 0; j < ny; ++j) {
            (*this)(0, j) = value;
        }
    }
    if (right) {
        for (Index j = 0; j < ny; ++j) {
            (*this)(nx-1, j) = value;
        }
    }
    if (bottom) {
        for (Index i = 0; i < nx; ++i) {
            (*this)(i, 0) = value;
        }
    }
    if (top) {
        for (Index i = 0; i < nx; ++i) {
            (*this)(i, ny-1) = value;
        }
    }
}

void ScalarField::applyNeumannBC(const Grid& grid, bool left, bool right, bool bottom, bool top) {
    Index nx = grid.getNx();
    Index ny = grid.getNy();
    
    if (left) {
        for (Index j = 0; j < ny; ++j) {
            (*this)(0, j) = (*this)(1, j);
        }
    }
    if (right) {
        for (Index j = 0; j < ny; ++j) {
            (*this)(nx-1, j) = (*this)(nx-2, j);
        }
    }
    if (bottom) {
        for (Index i = 0; i < nx; ++i) {
            (*this)(i, 0) = (*this)(i, 1);
        }
    }
    if (top) {
        for (Index i = 0; i < nx; ++i) {
            (*this)(i, ny-1) = (*this)(i, ny-2);
        }
    }
}

// VectorField implementation
VectorField::VectorField(const Grid& grid)
    : grid_(grid), u_(grid), v_(grid)
{
}

void VectorField::fill(Real uValue, Real vValue) {
    u_.fill(uValue);
    v_.fill(vValue);
}

void VectorField::copyFrom(const VectorField& other) {
    u_.copyFrom(other.u_);
    v_.copyFrom(other.v_);
}

void VectorField::applyLidDrivenCavityBC(Real lidVelocity) {
    Index nx = grid_.getNx();
    Index ny = grid_.getNy();
    
    // Bottom wall: u = v = 0
    for (Index i = 0; i < nx; ++i) {
        u_(i, 0) = 0.0;
        v_(i, 0) = 0.0;
    }
    
    // Top wall (lid): u = lidVelocity, v = 0
    for (Index i = 0; i < nx; ++i) {
        u_(i, ny-1) = lidVelocity;
        v_(i, ny-1) = 0.0;
    }
    
    // Left wall: u = v = 0
    for (Index j = 0; j < ny; ++j) {
        u_(0, j) = 0.0;
        v_(0, j) = 0.0;
    }
    
    // Right wall: u = v = 0
    for (Index j = 0; j < ny; ++j) {
        u_(nx-1, j) = 0.0;
        v_(nx-1, j) = 0.0;
    }
}
