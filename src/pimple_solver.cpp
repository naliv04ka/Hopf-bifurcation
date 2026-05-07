#include "pimple_solver.h"
#include <iostream>
#include <iomanip>

const Real RHO = 1.0;
const Real U_LID = 1.0;

PIMPLESolver::PIMPLESolver(Index nx, Index ny, Real Re)
    : nx_(nx), ny_(ny), Re_(Re)
{
    dx_ = 1.0 / nx_;
    dy_ = 1.0 / ny_;
    nu_ = U_LID * 1.0 / Re_;
    
    // Параметры PIMPLE (консервативные для стабильности)
    nOuter_ = 3;     // 3 внешних итерации SIMPLE
    nCorr_ = 1;      // 1 коррекция PISO
    alphaU_ = 0.7;   // Релаксация скорости
    alphaP_ = 0.3;   // Релаксация давления
    
    // Размещенная сетка:
    // u: (nx+1) × ny точек на вертикальных гранях
    // v: nx × (ny+1) точек на горизонтальных гранях
    // p: nx × ny точек в центрах ячеек
    
    u_.resize((nx_ + 1) * ny_, 0.0);
    v_.resize(nx_ * (ny_ + 1), 0.0);
    p_.resize(nx_ * ny_, 0.0);
    
    u_n_.resize((nx_ + 1) * ny_, 0.0);
    v_n_.resize(nx_ * (ny_ + 1), 0.0);
    
    u_work_.resize((nx_ + 1) * ny_, 0.0);
    v_work_.resize(nx_ * (ny_ + 1), 0.0);
    p_work_.resize(nx_ * ny_, 0.0);
    
    auP_.resize((nx_ + 1) * ny_, 0.0);
    avP_.resize(nx_ * (ny_ + 1), 0.0);
}

void PIMPLESolver::momentumPredictorU(Real alphaU) {
    for (Index i = 1; i < nx_; ++i) {
        for (Index j = 0; j < ny_; ++j) {
            Index idx = uIndex(i, j);
            
            // Интерполяция скоростей
            Real uc = u_work_[idx];
            Real uE = u_work_[uIndex(i + 1, j)];
            Real uW = u_work_[uIndex(i - 1, j)];
            
            Real ue = 0.5 * (uc + uE);
            Real uw = 0.5 * (uc + uW);
            
            Real vn = 0.5 * (v_work_[vIndex(i - 1, j + 1)] + v_work_[vIndex(i, j + 1)]);
            Real vs = 0.5 * (v_work_[vIndex(i - 1, j)] + v_work_[vIndex(i, j)]);
            
            // Массовые потоки
            Real Fe = RHO * ue * dy_;
            Real Fw = RHO * uw * dy_;
            Real Fn = RHO * vn * dx_;
            Real Fs = RHO * vs * dx_;
            
            // Коэффициенты (upwind + диффузия)
            Real aE = std::max(-Fe, 0.0) + nu_ * dy_ / dx_;
            Real aW = std::max(Fw, 0.0) + nu_ * dy_ / dx_;
            
            Real aN, aS, uN_ghost, uS_ghost;
            
            // Граничные условия
            if (j == ny_ - 1) {  // Верхняя крышка
                aN = std::max(-Fn, 0.0) + nu_ * dx_ / (0.5 * dy_);
                uN_ghost = 2.0 * U_LID - uc;
            } else {
                aN = std::max(-Fn, 0.0) + nu_ * dx_ / dy_;
                uN_ghost = u_work_[uIndex(i, j + 1)];
            }
            
            if (j == 0) {  // Нижняя стенка
                aS = std::max(Fs, 0.0) + nu_ * dx_ / (0.5 * dy_);
                uS_ghost = -uc;
            } else {
                aS = std::max(Fs, 0.0) + nu_ * dx_ / dy_;
                uS_ghost = u_work_[uIndex(i, j - 1)];
            }
            
            // Временной член (будет использовать переданный dt через advance())
            // Пока используем фиксированный dt=0.02
            Real dt_fixed = 0.02;
            Real aP_base = aE + aW + aN + aS + RHO * dx_ * dy_ / dt_fixed;
            
            // Оператор H (включая временной член от u_n)
            Real H = aE * uE + aW * uW + aN * uN_ghost + aS * uS_ghost
                    + (RHO * dx_ * dy_ / dt_fixed) * u_n_[idx];
            
            // Градиент давления
            Real dP = (p_work_[pIndex(i, j)] - p_work_[pIndex(i - 1, j)]) * dy_;
            
            // Предикторная скорость с под-релаксацией
            u_work_[idx] = (1.0 - alphaU) * u_work_[idx]
                          + (alphaU / aP_base) * (H - dP);
            
            // Эффективная диагональ
            auP_[idx] = aP_base / alphaU;
        }
    }
    
    // Граничные условия
    for (Index j = 0; j < ny_; ++j) {
        u_work_[uIndex(0, j)] = 0.0;      // Левая стенка
        u_work_[uIndex(nx_, j)] = 0.0;    // Правая стенка
    }
}

void PIMPLESolver::momentumPredictorV(Real alphaU) {
    for (Index i = 0; i < nx_; ++i) {
        for (Index j = 1; j < ny_; ++j) {
            Index idx = vIndex(i, j);
            
            // Интерполяция скоростей
            Real vc = v_work_[idx];
            Real vN = v_work_[vIndex(i, j + 1)];
            Real vS = v_work_[vIndex(i, j - 1)];
            
            Real vn = 0.5 * (vc + vN);
            Real vsf = 0.5 * (vc + vS);
            
            Real ue = 0.5 * (u_work_[uIndex(i + 1, j - 1)] + u_work_[uIndex(i + 1, j)]);
            Real uw = 0.5 * (u_work_[uIndex(i, j - 1)] + u_work_[uIndex(i, j)]);
            
            // Массовые потоки
            Real Fe = RHO * ue * dy_;
            Real Fw = RHO * uw * dy_;
            Real Fn = RHO * vn * dx_;
            Real Fs = RHO * vsf * dx_;
            
            // Коэффициенты
            Real aN = std::max(-Fn, 0.0) + nu_ * dx_ / dy_;
            Real aS = std::max(Fs, 0.0) + nu_ * dx_ / dy_;
            
            Real aE, aW, vE_ghost, vW_ghost;
            
            // Граничные условия
            if (i == nx_ - 1) {  // Правая стенка
                aE = std::max(-Fe, 0.0) + nu_ * dy_ / (0.5 * dx_);
                vE_ghost = -vc;
            } else {
                aE = std::max(-Fe, 0.0) + nu_ * dy_ / dx_;
                vE_ghost = v_work_[vIndex(i + 1, j)];
            }
            
            if (i == 0) {  // Левая стенка
                aW = std::max(Fw, 0.0) + nu_ * dy_ / (0.5 * dx_);
                vW_ghost = -vc;
            } else {
                aW = std::max(Fw, 0.0) + nu_ * dy_ / dx_;
                vW_ghost = v_work_[vIndex(i - 1, j)];
            }
            
            // Временной член
            Real dt_fixed = 0.02;
            Real aP_base = aE + aW + aN + aS + RHO * dx_ * dy_ / dt_fixed;
            
            // Оператор H (включая временной член от v_n)
            Real H = aE * vE_ghost + aW * vW_ghost + aN * vN + aS * vS
                    + (RHO * dx_ * dy_ / dt_fixed) * v_n_[idx];
            
            // Градиент давления
            Real dP = (p_work_[pIndex(i, j)] - p_work_[pIndex(i, j - 1)]) * dx_;
            
            // Предикторная скорость
            v_work_[idx] = (1.0 - alphaU) * v_work_[idx]
                          + (alphaU / aP_base) * (H - dP);
            
            // Эффективная диагональ
            avP_[idx] = aP_base / alphaU;
        }
    }
    
    // Граничные условия
    for (Index i = 0; i < nx_; ++i) {
        v_work_[vIndex(i, 0)] = 0.0;      // Нижняя стенка
        v_work_[vIndex(i, ny_)] = 0.0;    // Верхняя стенка
    }
}

void PIMPLESolver::solvePressure(Real alphaP) {
    std::vector<Real> p_corr(nx_ * ny_, 0.0);
    
    // Гаусс-Зейдель с 100 итерациями
    for (Index iter = 0; iter < 100; ++iter) {
        for (Index i = 0; i < nx_; ++i) {
            for (Index j = 0; j < ny_; ++j) {
                Index idx = pIndex(i, j);
                
                // Фиксация давления в точке (0,0)
                if (i == 0 && j == 0) {
                    p_corr[idx] = 0.0;
                    continue;
                }
                
                // Коэффициенты уравнения Пуассона
                Real aE = (i < nx_ - 1) ? (dy_ * dy_ / auP_[uIndex(i + 1, j)]) : 0.0;
                Real aW = (i > 0) ? (dy_ * dy_ / auP_[uIndex(i, j)]) : 0.0;
                Real aN = (j < ny_ - 1) ? (dx_ * dx_ / avP_[vIndex(i, j + 1)]) : 0.0;
                Real aS = (j > 0) ? (dx_ * dx_ / avP_[vIndex(i, j)]) : 0.0;
                Real aP = aE + aW + aN + aS;
                
                // Невязка массы
                Real b = -((u_work_[uIndex(i + 1, j)] - u_work_[uIndex(i, j)]) * dy_
                         + (v_work_[vIndex(i, j + 1)] - v_work_[vIndex(i, j)]) * dx_);
                
                Real pE = (i < nx_ - 1) ? p_corr[pIndex(i + 1, j)] : 0.0;
                Real pW = (i > 0) ? p_corr[pIndex(i - 1, j)] : 0.0;
                Real pN = (j < ny_ - 1) ? p_corr[pIndex(i, j + 1)] : 0.0;
                Real pS = (j > 0) ? p_corr[pIndex(i, j - 1)] : 0.0;
                
                if (aP > 1e-12) {
                    p_corr[idx] = (aE * pE + aW * pW + aN * pN + aS * pS + b) / aP;
                }
            }
        }
    }
    
    // Коррекция давления с релаксацией
    for (Index i = 0; i < nx_; ++i) {
        for (Index j = 0; j < ny_; ++j) {
            p_work_[pIndex(i, j)] += alphaP * p_corr[pIndex(i, j)];
        }
    }
    
    // Коррекция скоростей (без релаксации!)
    for (Index i = 1; i < nx_; ++i) {
        for (Index j = 0; j < ny_; ++j) {
            Index idx = uIndex(i, j);
            u_work_[idx] -= (dy_ / auP_[idx]) * (p_corr[pIndex(i, j)] - p_corr[pIndex(i - 1, j)]);
        }
    }
    
    for (Index i = 0; i < nx_; ++i) {
        for (Index j = 1; j < ny_; ++j) {
            Index idx = vIndex(i, j);
            v_work_[idx] -= (dx_ / avP_[idx]) * (p_corr[pIndex(i, j)] - p_corr[pIndex(i, j - 1)]);
        }
    }
}

void PIMPLESolver::advance(Real dt) {
    // Сохранить поля на начало шага
    u_n_ = u_;
    v_n_ = v_;
    
    // Инициализировать рабочие массивы
    u_work_ = u_;
    v_work_ = v_;
    p_work_ = p_;
    
    // ВНЕШНИЙ ЦИКЛ (SIMPLE итерации)
    for (Index outer = 0; outer < nOuter_; ++outer) {
        Real alphaU_current, alphaP_current;
        
        if (outer < nOuter_ - 1) {
            // Промежуточные итерации - SIMPLE (с релаксацией)
            alphaU_current = alphaU_;
            alphaP_current = alphaP_;
        } else {
            // Последняя итерация - PISO (без релаксации)
            alphaU_current = 1.0;
            alphaP_current = 1.0;
        }
        
        // Предиктор импульса
        momentumPredictorU(alphaU_current);
        momentumPredictorV(alphaU_current);
        
        // ВНУТРЕННИЙ ЦИКЛ (PISO коррекции)
        for (Index corr = 0; corr < nCorr_; ++corr) {
            solvePressure(alphaP_current);
        }
    }
    
    // Обновить результаты
    u_ = u_work_;
    v_ = v_work_;
    p_ = p_work_;
}

void PIMPLESolver::getVelocity(std::vector<Real>& u_cell, std::vector<Real>& v_cell) const {
    u_cell.resize(nx_ * ny_);
    v_cell.resize(nx_ * ny_);
    
    for (Index i = 0; i < nx_; ++i) {
        for (Index j = 0; j < ny_; ++j) {
            // Интерполяция в центр ячейки
            u_cell[i * ny_ + j] = 0.5 * (u_[uIndex(i, j)] + u_[uIndex(i + 1, j)]);
            v_cell[i * ny_ + j] = 0.5 * (v_[vIndex(i, j)] + v_[vIndex(i, j + 1)]);
        }
    }
}

void PIMPLESolver::getPressure(std::vector<Real>& p_cell) const {
    p_cell = p_;
}

Real PIMPLESolver::getMaxDivergence() const {
    Real maxDiv = 0.0;
    
    for (Index i = 0; i < nx_; ++i) {
        for (Index j = 0; j < ny_; ++j) {
            Real div = (u_[uIndex(i + 1, j)] - u_[uIndex(i, j)]) * dy_
                     + (v_[vIndex(i, j + 1)] - v_[vIndex(i, j)]) * dx_;
            maxDiv = std::max(maxDiv, std::abs(div));
        }
    }
    
    return maxDiv;
}

Real PIMPLESolver::getMaxVelocity() const {
    Real maxU = 0.0;
    
    for (Index i = 0; i <= nx_; ++i) {
        for (Index j = 0; j < ny_; ++j) {
            maxU = std::max(maxU, std::abs(u_[uIndex(i, j)]));
        }
    }
    
    for (Index i = 0; i < nx_; ++i) {
        for (Index j = 0; j <= ny_; ++j) {
            maxU = std::max(maxU, std::abs(v_[vIndex(i, j)]));
        }
    }
    
    return maxU;
}

Real PIMPLESolver::getVelocityU(Real x, Real y) const {
    // Простая билинейная интерполяция
    Index i = static_cast<Index>(x / dx_);
    Index j = static_cast<Index>(y / dy_);
    
    i = std::max(0, std::min(i, nx_ - 1));
    j = std::max(0, std::min(j, ny_ - 1));
    
    return 0.5 * (u_[uIndex(i, j)] + u_[uIndex(i + 1, j)]);
}

Real PIMPLESolver::getVelocityV(Real x, Real y) const {
    Index i = static_cast<Index>(x / dx_);
    Index j = static_cast<Index>(y / dy_);
    
    i = std::max(0, std::min(i, nx_ - 1));
    j = std::max(0, std::min(j, ny_ - 1));
    
    return 0.5 * (v_[vIndex(i, j)] + v_[vIndex(i, j + 1)]);
}
