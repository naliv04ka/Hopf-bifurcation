#include "pimple_solver.h"
#include "fft_analyzer.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <chrono>

void saveResults(const PIMPLESolver& solver, Index nx, Index ny, 
                 const std::string& filename) {
    std::vector<Real> u_cell, v_cell, p_cell;
    solver.getVelocity(u_cell, v_cell);
    solver.getPressure(p_cell);
    
    std::ofstream out(filename);
    out << "x,y,u,v,p\n";
    
    Real dx = 1.0 / nx;
    Real dy = 1.0 / ny;
    
    for (Index i = 0; i < nx; ++i) {
        for (Index j = 0; j < ny; ++j) {
            Real x = (i + 0.5) * dx;
            Real y = (j + 0.5) * dy;
            
            out << x << "," << y << "," 
                << u_cell[i * ny + j] << "," 
                << v_cell[i * ny + j] << ","
                << p_cell[i * ny + j] << "\n";
        }
    }
    
    out.close();
    std::cout << "Results saved to " << filename << std::endl;
}

int main(int argc, char** argv) {
    std::cout << "===================================\n";
    std::cout << "  Lid-Driven Cavity Flow Solver   \n";
    std::cout << "  Improved PIMPLE Algorithm        \n";
    std::cout << "===================================\n\n";
    
    // Параметры симуляции
    Real Re = 1000.0;
    Index gridSize = 81;
    
    if (argc > 1) {
        Re = std::atof(argv[1]);
    }
    if (argc > 2) {
        gridSize = std::atoi(argv[2]);
    }
    
    std::cout << "Running simulation:\n";
    std::cout << "  Re = " << Re << "\n";
    std::cout << "  Grid = " << gridSize << "x" << gridSize << "\n\n";
    
    // Параметры времени
    Real dt = 0.02;
    Real totalTime = 10.0;
    Index numSteps = static_cast<Index>(totalTime / dt);
    Index outputInterval = 20;
    
    // Создать решатель
    PIMPLESolver solver(gridSize, gridSize, Re);
    
    // Мониторинг в центре
    Real monitorX = 0.5;
    Real monitorY = 0.5;
    
    std::ofstream monitorFile("velocity_monitor.dat");
    monitorFile << "# Time\tU\tV\tU_magnitude\n";
    
    auto startTime = std::chrono::high_resolution_clock::now();
    
    // Временной цикл
    for (Index step = 0; step < numSteps; ++step) {
        solver.advance(dt);
        
        Real currentTime = step * dt;
        
        // Мониторинг
        if (step % 10 == 0) {
            Real u = solver.getVelocityU(monitorX, monitorY);
            Real v = solver.getVelocityV(monitorX, monitorY);
            Real uMag = std::sqrt(u*u + v*v);
            
            monitorFile << currentTime << "\t" << u << "\t" << v << "\t" << uMag << "\n";
        }
        
        // Прогресс
        if (step % outputInterval == 0) {
            Real maxU = solver.getMaxVelocity();
            Real maxDiv = solver.getMaxDivergence();
            
            std::cout << "Step " << step << "/" << numSteps 
                     << " (" << std::fixed << std::setprecision(1) 
                     << (100.0 * step / numSteps) << "%)"
                     << " | t=" << std::setprecision(2) << currentTime
                     << " | max|U|=" << std::scientific << std::setprecision(3) << maxU
                     << " | max|div|=" << maxDiv << "\n";
        }
    }
    
    monitorFile.close();
    
    auto endTime = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(endTime - startTime);
    
    std::cout << "\nSimulation completed in " << duration.count() << " seconds\n";
    
    // Сохранить результаты
    saveResults(solver, gridSize, gridSize, "cavity_results.csv");
    
    // Вывести профиль U
    std::cout << "\nU-velocity profile at x=0.5:\n";
    std::cout << "y\tu\n";
    
    Index midX = gridSize / 2;
    Real dx = 1.0 / gridSize;
    Real dy = 1.0 / gridSize;
    
    std::vector<Real> u_cell, v_cell, p_cell;
    solver.getVelocity(u_cell, v_cell);
    
    for (Index j = gridSize - 1; j >= 0; --j) {
        Real y = (j + 0.5) * dy;
        std::cout << std::fixed << std::setprecision(4) << y << "\t" 
                 << u_cell[midX * gridSize + j] << "\n";
    }
    
    return 0;
}
