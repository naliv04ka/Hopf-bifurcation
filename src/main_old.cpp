#include "common.h"
#include "grid.h"
#include "field.h"
#include "navier_stokes_pimple.h"
#include "fft_analyzer.h"
#include "post_processor.h"
#include <chrono>
#include <sstream>
#include <map>
#include <filesystem>

namespace fs = std::filesystem;

struct SimulationResult {
    Real Reynolds;
    Index gridSize;
    bool isOscillatory;
    Real frequency;
    Real amplitude;
    std::string outputDir;
};

class CavityFlowSimulator {
public:
    CavityFlowSimulator() {}
    
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
        
        // Create directories using C++17 filesystem (cross-platform)
        fs::create_directories(result.outputDir);
        
        // Setup simulation parameters
        SimulationParams params;
        params.Reynolds = Re;
        params.nx = gridSize;
        params.ny = gridSize;
        params.dt = 0.005;
        params.totalTime = 200.0;
        
        // Create grid and solver
        Grid grid(params.nx, params.ny, params.Lx, params.Ly);
        NavierStokesPIMPLE solver(grid, params);
        PostProcessor postProc(grid);
        FFTAnalyzer fft;
        
        // Monitoring point (0.5, 0.5)
        Real monitorX = 0.5;
        Real monitorY = 0.5;
        
        Index numSteps = static_cast<Index>(params.totalTime / params.dt);
        Index outputInterval = 1000; // Output every 1000 steps
        Index monitorInterval = 10;  // Monitor every 10 steps
        
        // File for monitoring velocity
        std::ofstream monitorFile(result.outputDir + "/velocity_monitor.dat");
        monitorFile << "# Time\tU\tV\tU_magnitude\n";
        
        auto startTime = std::chrono::high_resolution_clock::now();
        
        // Time integration
        for (Index step = 0; step < numSteps; ++step) {
            solver.advance();
            
            // Monitor velocity at (0.5, 0.5)
            if (step % monitorInterval == 0) {
                Real u = solver.getVelocityU(monitorX, monitorY);
                Real v = solver.getVelocityV(monitorX, monitorY);
                Real uMag = std::sqrt(u*u + v*v);
                
                Real time = solver.getCurrentTime();
                monitorFile << time << "\t" << u << "\t" << v << "\t" << uMag << "\n";
                
                // Collect data for FFT (after initial transient)
                if (time > 50.0) {
                    fft.addDataPoint(time, u);
                }
            }
            
            // Progress output
            if (step % outputInterval == 0) {
                Real progress = 100.0 * step / numSteps;
                Real maxU = solver.getMaxVelocity();
                Real maxDiv = solver.getMaxDivergence();
                
                std::cout << "  Step " << step << "/" << numSteps 
                         << " (" << std::fixed << std::setprecision(1) << progress << "%)"
                         << " | t=" << solver.getCurrentTime()
                         << " | max|U|=" << std::scientific << std::setprecision(3) << maxU
                         << " | max|div|=" << maxDiv << "\n";
            }
            
            // Save VTK files periodically
            if (step % (outputInterval * 5) == 0) {
                std::ostringstream vtkName;
                vtkName << result.outputDir << "/flow_" 
                       << std::setfill('0') << std::setw(6) << step << ".vtk";
                postProc.writeVTK(vtkName.str(), solver.getVelocity(), 
                                 solver.getPressure(), solver.getCurrentTime());
            }
        }
        
        monitorFile.close();
        
        auto endTime = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(endTime - startTime);
        
        std::cout << "\nSimulation completed in " << duration.count() << " seconds\n";
        
        // Save final state
        postProc.writeVTK(result.outputDir + "/final_state.vtk",
                         solver.getVelocity(), solver.getPressure(), 
                         solver.getCurrentTime());
        
        postProc.writeStreamlines(result.outputDir + "/streamlines.dat",
                                 solver.getVelocity());
        
        // Perform FFT analysis
        std::cout << "\nPerforming FFT analysis...\n";
        fft.computeFFT();
        
        fft.saveTimeSeries(result.outputDir + "/time_series.dat");
        fft.saveSpectrum(result.outputDir + "/spectrum.dat");
        
        result.frequency = fft.getDominantFrequency();
        result.amplitude = fft.getDominantAmplitude();
        
        // Determine if oscillatory (check if dominant frequency is significant)
        // A flow is considered oscillatory if the amplitude is > 1% of mean flow
        result.isOscillatory = (result.amplitude > 0.01);
        
        std::cout << "Dominant frequency: " << result.frequency << " Hz\n";
        std::cout << "Dominant amplitude: " << result.amplitude << "\n";
        std::cout << "Flow type: " << (result.isOscillatory ? "OSCILLATORY" : "STATIONARY") << "\n";
        
        return result;
    }
    
    void runParametricStudy() {
        std::string baseDir = "results";
        
        // Create base directory using C++17 filesystem (cross-platform)
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
        
        // Analyze results and find critical Reynolds number for each grid
        std::cout << "\n\n========================================\n";
        std::cout << "SUMMARY OF RESULTS\n";
        std::cout << "========================================\n\n";
        
        std::ofstream summaryFile(baseDir + "/summary.dat");
        summaryFile << "# Re\tGrid\tOscillatory\tFrequency\tAmplitude\n";
        
        std::map<Index, std::pair<Real, Real>> criticalRe; // gridSize -> (Re_lower, Re_upper)
        
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
