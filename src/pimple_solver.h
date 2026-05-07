#ifndef PIMPLE_SOLVER_H
#define PIMPLE_SOLVER_H

#include "common.h"
#include <vector>
#include <cmath>
#include <algorithm>

class PIMPLESolver {
public:
    PIMPLESolver(Index nx, Index ny, Real Re);
    
    // Продвижение по времени
    void advance(Real dt);
    
    // Получение полей
    void getVelocity(std::vector<Real>& u_cell, std::vector<Real>& v_cell) const;
    void getPressure(std::vector<Real>& p_cell) const;
    
    // Статистика
    Real getMaxDivergence() const;
    Real getMaxVelocity() const;
    
    // Интерполяция для мониторинга
    Real getVelocityU(Real x, Real y) const;
    Real getVelocityV(Real x, Real y) const;
    
private:
    // Размеры сетки
    Index nx_, ny_;
    Real dx_, dy_;
    Real Re_, nu_;
    
    // Параметры PIMPLE
    Index nOuter_;      // Внешние итерации SIMPLE
    Index nCorr_;       // Внутренние коррекции PISO
    Real alphaU_;       // Релаксация скорости
    Real alphaP_;       // Релаксация давления
    
    // Поля на размещенной сетке
    std::vector<Real> u_;      // u на вертикальных гранях (nx+1) × ny
    std::vector<Real> v_;      // v на горизонтальных гранях nx × (ny+1)
    std::vector<Real> p_;      // p в центрах ячеек nx × ny
    
    std::vector<Real> u_n_;    // u на начало шага
    std::vector<Real> v_n_;    // v на начало шага
    
    std::vector<Real> u_work_; // Рабочие массивы
    std::vector<Real> v_work_;
    std::vector<Real> p_work_;
    
    std::vector<Real> auP_;    // Эффективные диагонали
    std::vector<Real> avP_;
    
    // Вспомогательные функции
    void momentumPredictorU(Real alphaU);
    void momentumPredictorV(Real alphaU);
    void solvePressure(Real alphaP);
    
    // Индексация
    Index uIndex(Index i, Index j) const { return i * ny_ + j; }
    Index vIndex(Index i, Index j) const { return i * (ny_ + 1) + j; }
    Index pIndex(Index i, Index j) const { return i * ny_ + j; }
};

#endif // PIMPLE_SOLVER_H
