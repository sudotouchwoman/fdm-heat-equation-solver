#pragma once

#include "model.hpp"

namespace solver {
    using diagonal = std::vector<double>;
    using tridiagonal_mx_extended = std::array<diagonal, 4>;

    // TDMA stands for tridiagonal matrix algorithm
    // aka Thomas algorithm in en literature
    // this one can solve sets of linear equations (SLEs)
    // with tridiagonal matrix in linear time
    class TDMA {
    private:
        diagonal c_star;
        diagonal d_star;
    private:
    public:
        TDMA() = default;
        TDMA(const size_t diagonal_length);
        ~TDMA() = default;
        void solve(const tridiagonal_mx_extended & newSLE, diagonal & storage);
    };

    // Problem entity wraps everything, i.e. the model and solvers
    // (also allocates some auxiliary storage once for the run)
    // problem is solved in an iterative manner, with grid's current
    // values updated at each step
    class Problem {
    private:
        size_t current_step = 0;
        const size_t n_iters = 0;
        model::IModel & m;
        TDMA solver_x;
        TDMA solver_y;
        tridiagonal_mx_extended mx_x;
        tridiagonal_mx_extended mx_y;
        diagonal f_x;
        diagonal f_y;
    private:
        void update_grid_row(const size_t y);
        void update_grid_col(const size_t x);
    public:
        Problem() = delete;
        Problem(model::IModel & model, const size_t n_iters);
        void step();
    };
}
