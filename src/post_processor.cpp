#include "post_processor.h"
#include <fstream>
#include <iomanip>

PostProcessor::PostProcessor(const Grid& grid)
    : grid_(grid)
{
}

void PostProcessor::writeVTK(const std::string& filename,
                            const VectorField& velocity,
                            const ScalarField& pressure,
                            Real time)
{
    std::ofstream file(filename);
    
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return;
    }
    
    Index nx = grid_.getNx();
    Index ny = grid_.getNy();
    
    // VTK header
    file << "# vtk DataFile Version 3.0\n";
    file << "Cavity flow at t=" << time << "\n";
    file << "ASCII\n";
    file << "DATASET STRUCTURED_POINTS\n";
    file << "DIMENSIONS " << nx << " " << ny << " 1\n";
    file << "ORIGIN 0 0 0\n";
    file << "SPACING " << grid_.getDx() << " " << grid_.getDy() << " 1\n";
    file << "POINT_DATA " << nx * ny << "\n";
    
    // Velocity field
    file << "VECTORS velocity double\n";
    for (Index j = 0; j < ny; ++j) {
        for (Index i = 0; i < nx; ++i) {
            file << velocity.u()(i, j) << " " << velocity.v()(i, j) << " 0.0\n";
        }
    }
    
    // Pressure field
    file << "SCALARS pressure double 1\n";
    file << "LOOKUP_TABLE default\n";
    for (Index j = 0; j < ny; ++j) {
        for (Index i = 0; i < nx; ++i) {
            file << pressure(i, j) << "\n";
        }
    }
    
    // Velocity magnitude
    file << "SCALARS velocity_magnitude double 1\n";
    file << "LOOKUP_TABLE default\n";
    for (Index j = 0; j < ny; ++j) {
        for (Index i = 0; i < nx; ++i) {
            Real u = velocity.u()(i, j);
            Real v = velocity.v()(i, j);
            Real mag = std::sqrt(u*u + v*v);
            file << mag << "\n";
        }
    }
    
    file.close();
}

Real PostProcessor::interpolateU(const VectorField& velocity, Real x, Real y) const {
    Real dx = grid_.getDx();
    Real dy = grid_.getDy();
    
    Index i = static_cast<Index>(x / dx);
    Index j = static_cast<Index>(y / dy);
    
    i = std::max(0, std::min(i, grid_.getNx() - 2));
    j = std::max(0, std::min(j, grid_.getNy() - 2));
    
    Real xi = (x - i * dx) / dx;
    Real eta = (y - j * dy) / dy;
    
    return (1.0 - xi) * (1.0 - eta) * velocity.u()(i, j) +
           xi * (1.0 - eta) * velocity.u()(i+1, j) +
           (1.0 - xi) * eta * velocity.u()(i, j+1) +
           xi * eta * velocity.u()(i+1, j+1);
}

Real PostProcessor::interpolateV(const VectorField& velocity, Real x, Real y) const {
    Real dx = grid_.getDx();
    Real dy = grid_.getDy();
    
    Index i = static_cast<Index>(x / dx);
    Index j = static_cast<Index>(y / dy);
    
    i = std::max(0, std::min(i, grid_.getNx() - 2));
    j = std::max(0, std::min(j, grid_.getNy() - 2));
    
    Real xi = (x - i * dx) / dx;
    Real eta = (y - j * dy) / dy;
    
    return (1.0 - xi) * (1.0 - eta) * velocity.v()(i, j) +
           xi * (1.0 - eta) * velocity.v()(i+1, j) +
           (1.0 - xi) * eta * velocity.v()(i, j+1) +
           xi * eta * velocity.v()(i+1, j+1);
}

void PostProcessor::integrateStreamline(const VectorField& velocity,
                                       Real x0, Real y0,
                                       std::vector<std::array<Real, 2>>& points,
                                       Real stepSize, Index maxSteps)
{
    points.clear();
    
    Real x = x0;
    Real y = y0;
    
    Real Lx = grid_.getLx();
    Real Ly = grid_.getLy();
    
    for (Index step = 0; step < maxSteps; ++step) {
        // Check if still in domain
        if (x < 0.0 || x > Lx || y < 0.0 || y > Ly) break;
        
        points.push_back({x, y});
        
        // RK4 integration
        Real u1 = interpolateU(velocity, x, y);
        Real v1 = interpolateV(velocity, x, y);
        
        Real u2 = interpolateU(velocity, x + 0.5*stepSize*u1, y + 0.5*stepSize*v1);
        Real v2 = interpolateV(velocity, x + 0.5*stepSize*u1, y + 0.5*stepSize*v1);
        
        Real u3 = interpolateU(velocity, x + 0.5*stepSize*u2, y + 0.5*stepSize*v2);
        Real v3 = interpolateV(velocity, x + 0.5*stepSize*u2, y + 0.5*stepSize*v2);
        
        Real u4 = interpolateU(velocity, x + stepSize*u3, y + stepSize*v3);
        Real v4 = interpolateV(velocity, x + stepSize*u3, y + stepSize*v3);
        
        Real u = (u1 + 2.0*u2 + 2.0*u3 + u4) / 6.0;
        Real v = (v1 + 2.0*v2 + 2.0*v3 + v4) / 6.0;
        
        Real mag = std::sqrt(u*u + v*v);
        if (mag < 1e-10) break; // Stagnation point
        
        x += stepSize * u / mag;
        y += stepSize * v / mag;
    }
}

void PostProcessor::writeStreamlines(const std::string& filename,
                                    const VectorField& velocity)
{
    std::ofstream file(filename);
    
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return;
    }
    
    Real Lx = grid_.getLx();
    Real Ly = grid_.getLy();
    
    // Generate streamlines from different starting points
    Index numStreamlines = 20;
    Real stepSize = 0.01;
    Index maxSteps = 1000;
    
    file << "# Streamlines\n";
    file << "# Number of streamlines: " << numStreamlines << "\n";
    
    for (Index n = 0; n < numStreamlines; ++n) {
        Real y0 = Ly * (n + 1.0) / (numStreamlines + 1.0);
        Real x0 = 0.05 * Lx; // Start near left wall
        
        std::vector<std::array<Real, 2>> points;
        integrateStreamline(velocity, x0, y0, points, stepSize, maxSteps);
        
        file << "# Streamline " << n << "\n";
        file << "# Points: " << points.size() << "\n";
        
        for (const auto& pt : points) {
            file << pt[0] << "\t" << pt[1] << "\n";
        }
        
        file << "\n\n";
    }
    
    file.close();
}

void PostProcessor::writeVelocityField(const std::string& filename,
                                      const VectorField& velocity)
{
    std::ofstream file(filename);
    
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return;
    }
    
    Index nx = grid_.getNx();
    Index ny = grid_.getNy();
    
    file << "# X\tY\tU\tV\n";
    
    for (Index j = 0; j < ny; ++j) {
        for (Index i = 0; i < nx; ++i) {
            Real x = grid_.getX(i);
            Real y = grid_.getY(j);
            file << x << "\t" << y << "\t" 
                 << velocity.u()(i, j) << "\t" 
                 << velocity.v()(i, j) << "\n";
        }
    }
    
    file.close();
}

void PostProcessor::writePressureField(const std::string& filename,
                                      const ScalarField& pressure)
{
    std::ofstream file(filename);
    
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return;
    }
    
    Index nx = grid_.getNx();
    Index ny = grid_.getNy();
    
    file << "# X\tY\tP\n";
    
    for (Index j = 0; j < ny; ++j) {
        for (Index i = 0; i < nx; ++i) {
            Real x = grid_.getX(i);
            Real y = grid_.getY(j);
            file << x << "\t" << y << "\t" 
                 << pressure(i, j) << "\n";
        }
    }
    
    file.close();
}

void PostProcessor::computeStreamfunction(const VectorField& velocity, ScalarField& psi) {
    Index nx = grid_.getNx();
    Index ny = grid_.getNy();
    Real dx = grid_.getDx();
    Real dy = grid_.getDy();
    
    // Initialize
    psi.fill(0.0);
    
    // Integrate from bottom-left corner
    // ψ is defined such that u = ∂ψ/∂y and v = -∂ψ/∂x
    
    // Integrate along bottom edge (j=0)
    for (Index i = 1; i < nx; ++i) {
        psi(i, 0) = psi(i-1, 0) - velocity.v()(i, 0) * dx;
    }
    
    // Integrate upward for each column
    for (Index i = 0; i < nx; ++i) {
        for (Index j = 1; j < ny; ++j) {
            psi(i, j) = psi(i, j-1) + velocity.u()(i, j) * dy;
        }
    }
}
