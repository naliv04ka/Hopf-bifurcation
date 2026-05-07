#include "grid.h"

Grid::Grid(Index nx, Index ny, Real Lx, Real Ly)
    : nx_(nx), ny_(ny), Lx_(Lx), Ly_(Ly)
{
    dx_ = Lx_ / (nx_ - 1);
    dy_ = Ly_ / (ny_ - 1);
}
