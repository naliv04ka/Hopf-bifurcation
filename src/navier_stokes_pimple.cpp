#include "navier_stokes_pimple.h"
#include <cmath>
#include <algorithm>

NavierStokesPIMPLE::NavierStokesPIMPLE(const Grid& grid, const SimulationParams& params)
    : grid_(grid), params_(params),
      velocity_(grid), velocityOld_(grid), velocityStar_(grid),
      pressure_(grid), pressureCorr_(grid),
      currentTime_(0.0), currentStep_(0)
{
    // Initialize fields
    velocity_.fill(0.0, 0.0);
    velocityOld_.fill(0.0, 0.0);
    pressure_.fill(0.0);
    
    // Apply initial boundary conditions
    velocity_.applyLidDrivenCavityBC();
}

void NavierStokesPIMPLE::advance() {
    // Store old velocity
    velocityOld_.copyFrom(velocity_);
    
    // PIMPLE loop
    for (Index corr = 0; corr < params_.nCorrectors; ++corr) {
        // Momentum predictor
        momentumPredictor();
        
        // Pressure-velocity coupling
        for (Index nonOrth = 0; nonOrth <= params_.nNonOrthCorr; ++nonOrth) {
            pressureCorrection();
        }
        
        // Velocity correction
        velocityCorrection();
    }
    
    // Update time
    currentTime_ += params_.dt;
    currentStep_++;
}

void NavierStokesPIMPLE::momentumPredictor() {
    Index nx = grid_.getNx();
    Index ny = grid_.getNy();
    Real dx = grid_.getDx();
    Real dy = grid_.getDy();
    Real dt = params_.dt;
    Real Re = params_.Reynolds;
    Real nu = 1.0 / Re;
    
    VectorField convTerm(grid_);
    VectorField diffTerm(grid_);
    
    // Compute convective and diffusive terms
    computeConvectiveTerm(velocityOld_, convTerm);
    computeDiffusiveTerm(velocityOld_, diffTerm);
    
    // Predict velocity (without pressure gradient)
    for (Index j = 1; j < ny-1; ++j) {
        for (Index i = 1; i < nx-1; ++i) {
            // u-component: ∂u/∂t = -∇·(uu) + ν∇²u - ∂p/∂x
            Real dpdx = (pressure_(i+1, j) - pressure_(i-1, j)) / (2.0 * dx);
            
            velocityStar_.u()(i, j) = velocityOld_.u()(i, j) + dt * (
                -convTerm.u()(i, j) +
                nu * diffTerm.u()(i, j) -
                dpdx
            );
            
            // v-component: ∂v/∂t = -∇·(uv) + ν∇²v - ∂p/∂y
            Real dpdy = (pressure_(i, j+1) - pressure_(i, j-1)) / (2.0 * dy);
            
            velocityStar_.v()(i, j) = velocityOld_.v()(i, j) + dt * (
                -convTerm.v()(i, j) +
                nu * diffTerm.v()(i, j) -
                dpdy
            );
        }
    }
    
    // Apply boundary conditions to predicted velocity
    velocityStar_.applyLidDrivenCavityBC();
}

void NavierStokesPIMPLE::computeConvectiveTerm(const VectorField& u, VectorField& conv) {
    Index nx = grid_.getNx();
    Index ny = grid_.getNy();
    Real dx = grid_.getDx();
    Real dy = grid_.getDy();
    
    // Use central differencing for convective term
    for (Index j = 1; j < ny-1; ++j) {
        for (Index i = 1; i < nx-1; ++i) {
            // ∂(u²)/∂x
            Real dudx = (u.u()(i+1, j) * u.u()(i+1, j) - u.u()(i-1, j) * u.u()(i-1, j)) / (2.0 * dx);
            
            // ∂(uv)/∂y
            Real uvAvg_jp = 0.25 * (u.u()(i, j) + u.u()(i, j+1)) * (u.v()(i, j) + u.v()(i+1, j));
            Real uvAvg_jm = 0.25 * (u.u()(i, j-1) + u.u()(i, j)) * (u.v()(i, j-1) + u.v()(i+1, j-1));
            Real duvdy = (uvAvg_jp - uvAvg_jm) / dy;
            
            conv.u()(i, j) = dudx + duvdy;
            
            // ∂(uv)/∂x
            Real uvAvg_ip = 0.25 * (u.u()(i, j) + u.u()(i, j+1)) * (u.v()(i, j) + u.v()(i+1, j));
            Real uvAvg_im = 0.25 * (u.u()(i-1, j) + u.u()(i-1, j+1)) * (u.v()(i-1, j) + u.v()(i, j));
            Real duvdx = (uvAvg_ip - uvAvg_im) / dx;
            
            // ∂(v²)/∂y
            Real dvdy = (u.v()(i, j+1) * u.v()(i, j+1) - u.v()(i, j-1) * u.v()(i, j-1)) / (2.0 * dy);
            
            conv.v()(i, j) = duvdx + dvdy;
        }
    }
}

void NavierStokesPIMPLE::computeDiffusiveTerm(const VectorField& u, VectorField& diff) {
    Index nx = grid_.getNx();
    Index ny = grid_.getNy();
    Real dx = grid_.getDx();
    Real dy = grid_.getDy();
    Real dx2 = dx * dx;
    Real dy2 = dy * dy;
    
    // Laplacian using central differences
    for (Index j = 1; j < ny-1; ++j) {
        for (Index i = 1; i < nx-1; ++i) {
            // ∇²u
            diff.u()(i, j) = (u.u()(i+1, j) - 2.0*u.u()(i, j) + u.u()(i-1, j)) / dx2 +
                             (u.u()(i, j+1) - 2.0*u.u()(i, j) + u.u()(i, j-1)) / dy2;
            
            // ∇²v
            diff.v()(i, j) = (u.v()(i+1, j) - 2.0*u.v()(i, j) + u.v()(i-1, j)) / dx2 +
                             (u.v()(i, j+1) - 2.0*u.v()(i, j) + u.v()(i, j-1)) / dy2;
        }
    }
}

void NavierStokesPIMPLE::pressureCorrection() {
    Index nx = grid_.getNx();
    Index ny = grid_.getNy();
    Real dx = grid_.getDx();
    Real dy = grid_.getDy();
    Real dt = params_.dt;
    
    // Build right-hand side: RHS = -∇·u* / dt
    std::vector<Real> rhs(grid_.getNumCells(), 0.0);
    
    for (Index j = 1; j < ny-1; ++j) {
        for (Index i = 1; i < nx-1; ++i) {
            Real divU = computeDivergence(velocityStar_, i, j);
            Index idx = grid_.getIndex(i, j);
            rhs[idx] = -divU / dt;
        }
    }
    
    // Solve pressure Poisson equation: ∇²p' = RHS
    solvePressurePoisson(pressureCorr_);
    
    // Update pressure: p = p + αₚ·p'
    for (Index idx = 0; idx < grid_.getNumCells(); ++idx) {
        pressure_[idx] += params_.pRelaxation * pressureCorr_[idx];
    }
    
    // Apply pressure boundary conditions (Neumann)
    pressure_.applyNeumannBC(grid_, true, true, true, true);
}

void NavierStokesPIMPLE::solvePressurePoisson(ScalarField& p) {
    Index nx = grid_.getNx();
    Index ny = grid_.getNy();
    Real dx = grid_.getDx();
    Real dy = grid_.getDy();
    Real dx2 = dx * dx;
    Real dy2 = dy * dy;
    
    // Gauss-Seidel iterations
    Real omega = 1.5; // SOR parameter
    Index maxIter = 1000;
    Real tolerance = params_.pTolerance;
    
    // Build divergence of predicted velocity as RHS
    std::vector<Real> rhs(grid_.getNumCells(), 0.0);
    for (Index j = 1; j < ny-1; ++j) {
        for (Index i = 1; i < nx-1; ++i) {
            Index idx = grid_.getIndex(i, j);
            rhs[idx] = -computeDivergence(velocityStar_, i, j) / params_.dt;
        }
    }
    
    p.fill(0.0);
    
    for (Index iter = 0; iter < maxIter; ++iter) {
        Real maxResidual = 0.0;
        
        for (Index j = 1; j < ny-1; ++j) {
            for (Index i = 1; i < nx-1; ++i) {
                Index idx = grid_.getIndex(i, j);
                
                Real pOld = p[idx];
                
                // Solve: (2/dx² + 2/dy²)p(i,j) = p(i+1,j)/dx² + p(i-1,j)/dx² + 
                //                                  p(i,j+1)/dy² + p(i,j-1)/dy² - rhs(i,j)
                Real coeff = 2.0 / dx2 + 2.0 / dy2;
                Real sum = (p(i+1, j) + p(i-1, j)) / dx2 +
                          (p(i, j+1) + p(i, j-1)) / dy2 -
                          rhs[idx];
                
                Real pNew = sum / coeff;
                
                // SOR
                p[idx] = pOld + omega * (pNew - pOld);
                
                maxResidual = std::max(maxResidual, std::abs(pNew - pOld));
            }
        }
        
        // Apply boundary conditions
        p.applyNeumannBC(grid_, true, true, true, true);
        
        if (maxResidual < tolerance) {
            break;
        }
    }
}

void NavierStokesPIMPLE::velocityCorrection() {
    Index nx = grid_.getNx();
    Index ny = grid_.getNy();
    Real dx = grid_.getDx();
    Real dy = grid_.getDy();
    Real dt = params_.dt;
    
    // Correct velocity: u = u* - dt·∇p'
    for (Index j = 1; j < ny-1; ++j) {
        for (Index i = 1; i < nx-1; ++i) {
            Real dpdx = (pressureCorr_(i+1, j) - pressureCorr_(i-1, j)) / (2.0 * dx);
            Real dpdy = (pressureCorr_(i, j+1) - pressureCorr_(i, j-1)) / (2.0 * dy);
            
            velocity_.u()(i, j) = velocityStar_.u()(i, j) - dt * dpdx;
            velocity_.v()(i, j) = velocityStar_.v()(i, j) - dt * dpdy;
        }
    }
    
    // Apply boundary conditions
    velocity_.applyLidDrivenCavityBC();
}

Real NavierStokesPIMPLE::computeDivergence(const VectorField& u, Index i, Index j) const {
    Real dx = grid_.getDx();
    Real dy = grid_.getDy();
    
    Real dudx = (u.u()(i+1, j) - u.u()(i-1, j)) / (2.0 * dx);
    Real dvdy = (u.v()(i, j+1) - u.v()(i, j-1)) / (2.0 * dy);
    
    return dudx + dvdy;
}

Real NavierStokesPIMPLE::interpolate(const ScalarField& field, Real x, Real y) const {
    Real dx = grid_.getDx();
    Real dy = grid_.getDy();
    
    Index i = static_cast<Index>(x / dx);
    Index j = static_cast<Index>(y / dy);
    
    // Clamp to valid range
    i = std::max(0, std::min(i, grid_.getNx() - 2));
    j = std::max(0, std::min(j, grid_.getNy() - 2));
    
    Real xi = (x - i * dx) / dx;
    Real eta = (y - j * dy) / dy;
    
    // Bilinear interpolation
    return (1.0 - xi) * (1.0 - eta) * field(i, j) +
           xi * (1.0 - eta) * field(i+1, j) +
           (1.0 - xi) * eta * field(i, j+1) +
           xi * eta * field(i+1, j+1);
}

Real NavierStokesPIMPLE::getVelocityU(Real x, Real y) const {
    return interpolate(velocity_.u(), x, y);
}

Real NavierStokesPIMPLE::getVelocityV(Real x, Real y) const {
    return interpolate(velocity_.v(), x, y);
}

Real NavierStokesPIMPLE::getMaxVelocity() const {
    Real maxU = 0.0;
    const auto& uData = velocity_.u().getData();
    const auto& vData = velocity_.v().getData();
    
    for (size_t i = 0; i < uData.size(); ++i) {
        Real magU = std::sqrt(uData[i]*uData[i] + vData[i]*vData[i]);
        maxU = std::max(maxU, magU);
    }
    
    return maxU;
}

Real NavierStokesPIMPLE::getMaxDivergence() const {
    Real maxDiv = 0.0;
    Index nx = grid_.getNx();
    Index ny = grid_.getNy();
    
    for (Index j = 1; j < ny-1; ++j) {
        for (Index i = 1; i < nx-1; ++i) {
            Real div = std::abs(computeDivergence(velocity_, i, j));
            maxDiv = std::max(maxDiv, div);
        }
    }
    
    return maxDiv;
}
