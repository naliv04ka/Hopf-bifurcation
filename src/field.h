#ifndef FIELD_H
#define FIELD_H

#include "common.h"
#include "grid.h"

// Scalar field class
class ScalarField {
public:
    ScalarField(const Grid& grid, Real initValue = 0.0);
    
    Real& operator()(Index i, Index j);
    const Real& operator()(Index i, Index j) const;
    
    Real& operator[](Index idx) { return data_[idx]; }
    const Real& operator[](Index idx) const { return data_[idx]; }
    
    void fill(Real value);
    void copyFrom(const ScalarField& other);
    
    Index size() const { return data_.size(); }
    const std::vector<Real>& getData() const { return data_; }
    std::vector<Real>& getData() { return data_; }
    
    // Apply boundary conditions
    void applyDirichletBC(Real value, bool left, bool right, bool bottom, bool top);
    void applyNeumannBC(const Grid& grid, bool left, bool right, bool bottom, bool top);
    
private:
    const Grid& grid_;
    std::vector<Real> data_;
};

// Vector field class (stores u and v components)
class VectorField {
public:
    VectorField(const Grid& grid);
    
    ScalarField& u() { return u_; }
    const ScalarField& u() const { return u_; }
    
    ScalarField& v() { return v_; }
    const ScalarField& v() const { return v_; }
    
    void fill(Real uValue, Real vValue);
    void copyFrom(const VectorField& other);
    
    // Apply lid-driven cavity boundary conditions
    void applyLidDrivenCavityBC(Real lidVelocity = 1.0);
    
private:
    const Grid& grid_;
    ScalarField u_;
    ScalarField v_;
};

#endif // FIELD_H
