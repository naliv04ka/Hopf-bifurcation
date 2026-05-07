#include "pimple_solver.h"
#include "fft_analyzer.h"
#include "post_processor.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <sstream>
#include <filesystem>
#include <map>

namespace fs = std::filesystem;

struct SimulationResult {
    Real Reynolds;
    Index gridSize;
    bool isOscillatory;
    Real frequency;
    Real amplitude;
    std::string outputDir;
};

// Адаптер для PostProcessor - конвертирует поля PIMPLE в формат для визуализации
class VTKWriter {
public:
    static void writeVTK(const std::string& filename, const PIMPLESolver& solver,
                        Index nx, Index ny, Real time) {
        std::vector<Real> u_cell, v_cell, p_cell;
        solver.getVelocity(u_cell, v_cell);
        solver.getPressure(p_cell);
        
        std::ofstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Failed to open file: " << filename << std::endl;
            return;
        }
        
        Real dx = 1.0 / nx;
        Real dy = 1.0 / ny;
        
        // VTK header
        file << "# vtk DataFile Version 3.0\n";
        file << "Cavity flow at t=" << time << "\n";
        file << "ASCII\n";
        file << "DATASET STRUCTURED_POINTS\n";
        file << "DIMENSIONS " << nx << " " << ny << " 1\n";
        file << "ORIGIN 0 0 0\n";
        file << "SPACING " << dx << " " << dy << " 1\n";
        file << "POINT_DATA " << nx * ny << "\n";
        
        // Velocity field
        file << "VECTORS velocity double\n";
        for (Index j = 0; j < ny; ++j) {
            for (Index i = 0; i < nx; ++i) {
                file << u_cell[i * ny + j] << " " << v_cell[i * ny + j] << " 0.0\n";
            }
        }
        
        // Pressure field
        file << "SCALARS pressure double 1\n";
        file << "LOOKUP_TABLE default\n";
        for (Index j = 0; j < ny; ++j) {
            for (Index i = 0; i < nx; ++i) {
                file << p_cell[i * ny + j] << "\n";
            }
        }
        
        // Velocity magnitude
        file << "SCALARS velocity_magnitude double 1\n";
        file << "LOOKUP_TABLE default\n";
        for (Index j = 0; j < ny; ++j) {
            for (Index i = 0; i < nx; ++i) {
                Real u = u_cell[i * ny + j];
                Real v = v_cell[i * ny + j];
                Real mag = std::sqrt(u*u + v*v);
                file << mag << "\n";
            }
        }
        
        file.close();
    }
};

class CavityFlowSimulator {
public:
    SimulationResult runSimulation(Real Re, Index gridSize, const std::string& baseDir) {
        std::cout << "\n========================================\n";
        std::cout << "Running simulation:\n";
        std::cout << "  Re = " << Re << "\n";
        std::cout << "  Grid = " << gridSize << "x" << gridSize << "\n";
        std::cout << "========================================\n";
        
        SimulationResult result;
        result.Reynolds = Re;
        result.gridSize = gridSize;
        result.isOscillatory = false;
        result.frequency = 0.0;
        result.amplitude = 0.0;
        
        // Create output directory
        std::ostringstream oss;
        oss << baseDir << "/Re" << static_cast<int>(Re) << "_grid" << gridSize;
        result.outputDir = oss.str();
        fs::create_directories(result.outputDir);
        
        // Simulation parameters
        Real dt = 0.005;
        Real totalTime = 200.0;
        Index numSteps = static_cast<Index>(totalTime / dt);
        
        // Create solver
        PIMPLESolver solver(gridSize, gridSize, Re);
        FFTAnalyzer fft;
        
        // Monitoring point
        Real monitorX = 0.5;
        Real monitorY = 0.5;
        
        Index outputInterval = 500;  // Output every 1000 steps
        Index monitorInterval = 10;   // Monitor every 10 steps
        
        // Monitoring file
        std::ofstream monitorFile(result.outputDir + "/velocity_monitor.dat");
        monitorFile << "# Time\tU\tV\tU_magnitude\n";
        
        auto startTime = std::chrono::high_resolution_clock::now();
        
        // Time integration
        for (Index step = 0; step < numSteps; ++step) {
            solver.advance(dt);
            
            Real currentTime = step * dt;
            
            // Monitor velocity at (0.5, 0.5)
            if (step % monitorInterval == 0) {
                Real u = solver.getVelocityU(monitorX, monitorY);
                Real v = solver.getVelocityV(monitorX, monitorY);
                Real uMag = std::sqrt(u*u + v*v);
                
                monitorFile << currentTime << "\t" << u << "\t" << v << "\t" << uMag << "\n";
                
                // Collect data for FFT (after initial transient)
                if (currentTime > 50.0) {
                    fft.addDataPoint(currentTime, u);
                }
            }
            
            // Progress output
            if (step % outputInterval == 0) {
                Real progress = 100.0 * step / numSteps;
                Real maxU = solver.getMaxVelocity();
                Real maxDiv = solver.getMaxDivergence();
                
                std::cout << "  Step " << step << "/" << numSteps 
                         << " (" << std::fixed << std::setprecision(1) << progress << "%)"
                         << " | t=" << std::setprecision(2) << currentTime
                         << " | max|U|=" << std::scientific << std::setprecision(3) << maxU
                         << " | max|div|=" << maxDiv << "\n";
            }
            
            // Save VTK files periodically
            if (step % (outputInterval) == 0) {
                std::ostringstream vtkName;
                vtkName << result.outputDir << "/flow_" 
                       << std::setfill('0') << std::setw(6) << step << ".vtk";
                VTKWriter::writeVTK(vtkName.str(), solver, gridSize, gridSize, currentTime);
            }
        }
        
        monitorFile.close();
        
        auto endTime = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(endTime - startTime);
        
        std::cout << "\nSimulation completed in " << duration.count() << " seconds\n";
        
        // Save final state
        VTKWriter::writeVTK(result.outputDir + "/final_state.vtk",
                           solver, gridSize, gridSize, totalTime);
        
        // Perform FFT analysis
        std::cout << "\nPerforming FFT analysis...\n";
        fft.computeFFT();
        
        fft.saveTimeSeries(result.outputDir + "/time_series.dat");
        fft.saveSpectrum(result.outputDir + "/spectrum.dat");
        
        result.frequency = fft.getDominantFrequency();
        result.amplitude = fft.getDominantAmplitude();
        
        // Determine if oscillatory
        result.isOscillatory = (result.amplitude > 0.01);
        
        std::cout << "Dominant frequency: " << result.frequency << " Hz\n";
        std::cout << "Dominant amplitude: " << result.amplitude << "\n";
        std::cout << "Flow type: " << (result.isOscillatory ? "OSCILLATORY" : "STATIONARY") << "\n";
        
        return result;
    }
    
    void runParametricStudy() {
        std::string baseDir = "results";
        fs::create_directories(baseDir);
        
        // Reynolds numbers to test
        std::vector<Real> ReynoldsNumbers = {1000, 2000, 3000, 5000, 7500, 10000};
        
        // Grid sizes to test
        std::vector<Index> gridSizes = {81, 129, 257};
        
        std::vector<SimulationResult> allResults;
        
        // Run all simulations
        for (Index gridSize : gridSizes) {
            for (Real Re : ReynoldsNumbers) {
                SimulationResult result = runSimulation(Re, gridSize, baseDir);
                allResults.push_back(result);
            }
        }
        
        // Analyze results
        std::cout << "\n\n========================================\n";
        std::cout << "SUMMARY OF RESULTS\n";
        std::cout << "========================================\n\n";
        
        std::ofstream summaryFile(baseDir + "/summary.dat");
        summaryFile << "# Re\tGrid\tOscillatory\tFrequency\tAmplitude\n";
        
        for (const auto& result : allResults) {
            summaryFile << result.Reynolds << "\t" 
                       << result.gridSize << "\t"
                       << (result.isOscillatory ? "YES" : "NO") << "\t"
                       << result.frequency << "\t"
                       << result.amplitude << "\n";
            
            std::cout << "Re=" << result.Reynolds 
                     << " | Grid=" << result.gridSize << "x" << result.gridSize
                     << " | " << (result.isOscillatory ? "OSCILLATORY" : "STATIONARY");
            
            if (result.isOscillatory) {
                std::cout << " | f=" << result.frequency << " Hz";
            }
            std::cout << "\n";
        }
        
        summaryFile.close();
        
        // Find critical Reynolds number for each grid
        std::cout << "\n========================================\n";
        std::cout << "CRITICAL REYNOLDS NUMBERS\n";
        std::cout << "========================================\n";
        
        std::ofstream criticalFile(baseDir + "/critical_re.dat");
        criticalFile << "# GridSize\tRe_cr_lower\tRe_cr_upper\tRe_cr_estimate\tdx\n";
        
        for (Index gridSize : gridSizes) {
            Real ReLower = 0.0;
            Real ReUpper = 20000.0;
            
            for (const auto& result : allResults) {
                if (result.gridSize == gridSize) {
                    if (!result.isOscillatory && result.Reynolds > ReLower) {
                        ReLower = result.Reynolds;
                    }
                    if (result.isOscillatory && result.Reynolds < ReUpper) {
                        ReUpper = result.Reynolds;
                    }
                }
            }
            
            Real ReCritical = (ReLower + ReUpper) / 2.0;
            Real dx = 1.0 / (gridSize - 1);
            
            std::cout << "Grid " << gridSize << "x" << gridSize 
                     << ": Re_cr ∈ [" << ReLower << ", " << ReUpper << "]"
                     << " → Re_cr ≈ " << ReCritical << "\n";
            
            criticalFile << gridSize << "\t" << ReLower << "\t" << ReUpper 
                        << "\t" << ReCritical << "\t" << dx << "\n";
        }
        
        criticalFile.close();
        
        std::cout << "\n========================================\n";
        std::cout << "Results saved to: " << baseDir << "\n";
        std::cout << "========================================\n";
    }
};

int main(int argc, char** argv) {
    std::cout << "===================================\n";
    std::cout << "  Lid-Driven Cavity Flow Solver   \n";
    std::cout << "  PIMPLE Algorithm                \n";
    std::cout << "===================================\n\n";
    
    CavityFlowSimulator simulator;
    
    if (argc > 1) {
        // Single simulation mode
        Real Re = std::atof(argv[1]);
        Index gridSize = (argc > 2) ? std::atoi(argv[2]) : 129;
        
        simulator.runSimulation(Re, gridSize, "results");
    } else {
        // Full parametric study
        std::cout << "Running full parametric study...\n";
        std::cout << "This will take several hours.\n\n";
        
        simulator.runParametricStudy();
    }
    
    return 0;
}
