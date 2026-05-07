#ifndef COMMON_H
#define COMMON_H

#include <vector>
#include <array>
#include <cmath>
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <complex>
#include <memory>

// Common type definitions
using Real = double;
using Index = int;

// Mathematical constants
constexpr Real PI = 3.14159265358979323846;

// Simulation parameters structure
struct SimulationParams {
    Real Reynolds;          // Reynolds number
    Real dt;                // Time step
    Real totalTime;         // Total simulation time
    Index nx, ny;           // Grid dimensions
    Real Lx, Ly;            // Domain size
    
    // PIMPLE parameters
    Index nCorrectors;      // Number of PIMPLE correctors
    Index nNonOrthCorr;     // Number of non-orthogonal correctors
    Real pRelaxation;       // Pressure under-relaxation
    Real uRelaxation;       // Velocity under-relaxation
    
    // Convergence criteria
    Real pTolerance;
    Real uTolerance;
    
    SimulationParams() :
        Reynolds(1000.0),
        dt(0.002),              // Уменьшен шаг для стабильности
        totalTime(200.0),
        nx(129), ny(129),
        Lx(1.0), Ly(1.0),
        nCorrectors(3),         // Увеличено для лучшей сходимости
        nNonOrthCorr(1),        // Добавлена ортогональная коррекция
        pRelaxation(0.3),       // Хорошо для давления
        uRelaxation(0.5),       // Уменьшена релаксация скорости!
        pTolerance(1e-6),
        uTolerance(1e-6)
    {}
};

#endif // COMMON_H
