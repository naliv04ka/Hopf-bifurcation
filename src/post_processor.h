#ifndef POST_PROCESSOR_H
#define POST_PROCESSOR_H

#include "common.h"
#include "grid.h"
#include "field.h"
#include "navier_stokes_pimple.h"

class PostProcessor {
public:
    PostProcessor(const Grid& grid);
    
    // Write VTK file for visualization
    void writeVTK(const std::string& filename,
                  const VectorField& velocity,
                  const ScalarField& pressure,
                  Real time);
    
    // Compute and write streamlines
    void writeStreamlines(const std::string& filename,
                         const VectorField& velocity);
    
    // Write velocity field data
    void writeVelocityField(const std::string& filename,
                           const VectorField& velocity);
    
    // Write pressure field data
    void writePressureField(const std::string& filename,
                           const ScalarField& pressure);
    
    // Compute streamfunction (for 2D flows)
    void computeStreamfunction(const VectorField& velocity, ScalarField& psi);
    
private:
    const Grid& grid_;
    
    // Helper function for streamline integration
    void integrateStreamline(const VectorField& velocity,
                            Real x0, Real y0,
                            std::vector<std::array<Real, 2>>& points,
                            Real stepSize, Index maxSteps);
    
    Real interpolateU(const VectorField& velocity, Real x, Real y) const;
    Real interpolateV(const VectorField& velocity, Real x, Real y) const;
};

#endif // POST_PROCESSOR_H
