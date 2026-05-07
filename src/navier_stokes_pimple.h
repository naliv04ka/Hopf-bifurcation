#ifndef NAVIER_STOKES_PIMPLE_H
#define NAVIER_STOKES_PIMPLE_H

#include "common.h"
#include "grid.h"
#include "field.h"

class NavierStokesPIMPLE {
public:
    NavierStokesPIMPLE(const Grid& grid, const SimulationParams& params);
    
    // Time stepping
    void advance();
    
    // Get current fields
    const VectorField& getVelocity() const { return velocity_; }
    const ScalarField& getPressure() const { return pressure_; }
    
    // Get velocity at a point (with interpolation)
    Real getVelocityU(Real x, Real y) const;
    Real getVelocityV(Real x, Real y) const;
    
    // Get current time and time step
    Real getCurrentTime() const { return currentTime_; }
    Index getCurrentStep() const { return currentStep_; }
    
    // Statistics
    Real getMaxVelocity() const;
    Real getMaxDivergence() const;
    
private:
    const Grid& grid_;
    const SimulationParams& params_;
    
    // Fields
    VectorField velocity_;      // Current velocity field
    VectorField velocityOld_;   // Previous time step velocity
    VectorField velocityStar_;  // Predicted velocity
    ScalarField pressure_;      // Pressure field
    ScalarField pressureCorr_;  // Pressure correction
    
    // Time tracking
    Real currentTime_;
    Index currentStep_;
    
    // PIMPLE algorithm methods
    void momentumPredictor();
    void pressureCorrection();
    void velocityCorrection();
    
    // Discretization methods
    void computeConvectiveTerm(const VectorField& u, VectorField& conv);
    void computeDiffusiveTerm(const VectorField& u, VectorField& diff);
    void solvePressurePoisson(ScalarField& p);
    
    // Linear algebra solvers
    void solveGaussSeidel(ScalarField& field, const std::vector<Real>& rhs,
                         const std::vector<Real>& coeffDiag,
                         const std::vector<Real>& coeffOffDiag,
                         Real tolerance, Index maxIter);
    
    // Gradient and divergence
    void computeGradient(const ScalarField& phi, VectorField& grad);
    Real computeDivergence(const VectorField& u, Index i, Index j) const;
    
    // Interpolation
    Real interpolate(const ScalarField& field, Real x, Real y) const;
};

#endif // NAVIER_STOKES_PIMPLE_H
