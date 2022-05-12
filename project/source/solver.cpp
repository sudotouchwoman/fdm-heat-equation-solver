#include "solver.hpp"

#include <algorithm>
#include <iostream>

constexpr bool VERBOSE = false;

namespace solver {

    TDMA::TDMA(const size_t diagonal_length):
    c_star(diagonal_length, 0), d_star(diagonal_length, 0) {}

    void TDMA::solve(
        const tridiagonal_mx_extended & SLE,
        diagonal & storage
    ) {
        const diagonal a = SLE[0], b = SLE[1], c = SLE[2], d = SLE[3];
        const size_t N = b.size();

        if (N != a.size())
            throw std::runtime_error("dimension mismatch for a");
        if (N != c.size())
            throw std::runtime_error("dimension mismatch for c");
        if (N != d.size())
            throw std::runtime_error("dimension mismatch for d");        
        if (N != storage.size())
            throw std::runtime_error("dimension mismatch for storage");

        // reallocate memory for auxiliary
        // collections if the dimension has changed
        if (c_star.size() != N) {
            std::cerr
            << "changes c^* size: was "
            << c_star.size()
            << " now: "
            << N;
            c_star = diagonal(N, 0);
            d_star = diagonal(N, 0);
        }

        // update the coefficients in the first row
        c_star[0] = c[0] / b[0];
        d_star[0] = d[0] / b[0];

        // update other coefficients iteratively
        for (size_t i = 1; i < N; ++i) {
            const double w = 1.0 / (b[i] - a[i] * c_star[i-1]);
            c_star[i] = c[i] * w;
            d_star[i] = (d[i] - a[i] * d_star[i-1]) * w;
        }

        // store the solution in the storage
        for (size_t i = N - 1; i-- > 0; ) {
            storage[i] = d_star[i] - c_star[i] * d[i+1];
        }
    }

    static void pprint_tridiag_matrix(const tridiagonal_mx_extended & mx, std::ostream & out) {
        const size_t diag_length = mx[0].size();
        if (diag_length == 0) throw std::runtime_error("zero-length matrix");

        auto add_spaces = [&out](const size_t n) {
            for (size_t tabs = 0; tabs < n; ++tabs) {
                out << " ";
            }
        };

        out << "0-" << mx[3][0] << ":\t\t" << mx[1][0] << ' ' << mx[2][0] << '\n';
        for (size_t i = 1; i < diag_length - 1; ++i) {
            out << i << '-' << mx[3][i] << ":\t\t";
            add_spaces(i);
            out << mx[0][i] << ' ' << mx[1][i] << ' ' << mx[2][i] << '\n';
        }
        out << diag_length - 1 << '-' << mx[3][diag_length - 1] << ":\t\t";
        add_spaces(diag_length - 1);
        out << mx[0][diag_length - 1] << ' ' << mx[1][diag_length - 1] << '\n';

    }

    static void pprint_solution_row(const diagonal & d, std::ostream & os) {
        os << "solution: ";
        for (const auto & e: d)
            os << e << ' ';
        os << '\n';
    }

    Problem::Problem(model::IModel & model, const size_t n_iters):
        n_iters(n_iters),
        m(model),
        solver_x(model.x_dim()),
        solver_y(model.y_dim()),
        mx_x({
            // SLE contains 3 diagonals (a, b & c) and the RHS (d)
            // see https://quantstart.com/articles/Tridiagonal-Matrix-Solver-via-Thomas-Algorithm/
            diagonal(model.x_dim(), 0),
            diagonal(model.x_dim(), 0),
            diagonal(model.x_dim(), 0),
            diagonal(model.x_dim(), 0)
        }),
        mx_y({
            diagonal(model.y_dim(), 0),
            diagonal(model.y_dim(), 0),
            diagonal(model.y_dim(), 0),
            diagonal(model.y_dim(), 0)
        }),
        f_x(model.x_dim()),
        f_y(model.y_dim()) {}

    void Problem::step() {
        if (current_step++ == n_iters) throw std::runtime_error("out of iterations");
        // performs simulation step and stores the
        // result in the grid of the model
        const size_t x_dim = m.x_dim(), y_dim = m.y_dim();

        // first, solve the 1D subproblems in the horizontal direction
        // --> y_dim systems for each grid row
        // and update the current values at each node
        for (size_t y = 0; y < y_dim; ++y) {
            // solve SLE for row y
            // boundary conditions on edge:
            model::boundary_coefs bc = m.get_x_first_coefs(y);
            // unpack the duple into diagonals
            // once again, see https://quantstart.com/articles/Tridiagonal-Matrix-Solver-via-Thomas-Algorithm/
            mx_x[1][0] = bc[0];
            mx_x[2][0] = bc[1];
            mx_x[3][0] = m.get_RHS_coefs_x(0, y);

            for (size_t x = 1; x < x_dim - 1; ++x) {
                // unpack triples into diagonals + right-hand side into d
                const model::tridiag_coefs tc = m.get_x_coefs(x, y);
                mx_x[0][x] = tc[0];
                mx_x[1][x] = tc[1];
                mx_x[2][x] = tc[2];
                mx_x[3][x] = m.get_RHS_coefs_x(x, y);
            }

            bc = m.get_x_last_coefs(y);
            mx_x[0][x_dim - 1] = bc[0];
            mx_x[1][x_dim - 1] = bc[1];
            mx_x[3][x_dim - 1] = m.get_RHS_coefs_x(x_dim - 1, y);

            // call solver and update current values in the row
            // std::cout << "matrix " << y << "\n";
            solver_x.solve(mx_x, f_x);
            if (VERBOSE) {
                pprint_tridiag_matrix(mx_x, std::cout);
                pprint_solution_row(f_x, std::cout);
            }
            // std::getchar();
            update_grid_row(y);
        }

        // then solve in the vertical direction -->
        // x_dim systems for each grid column
        // and update the current values at each node
        for (size_t x = 0; x < x_dim; ++x) {
            // solve SLE for column x
            // boundary conditions on edge:
            model::boundary_coefs bc = m.get_y_first_coefs(x);
            // unpack the duple into diagonals
            // once again, see https://quantstart.com/articles/Tridiagonal-Matrix-Solver-via-Thomas-Algorithm/
            mx_y[1][0] = bc[0];
            mx_y[2][0] = bc[1];
            mx_y[3][0] = m.get_RHS_coefs_y(x, 0);

            for (size_t y = 1; y < y_dim - 1; ++y) {
                const model::tridiag_coefs tc = m.get_y_coefs(x, y);
                mx_y[0][y] = tc[0];
                mx_y[1][y] = tc[1];
                mx_y[2][y] = tc[2];
                mx_y[3][y] = m.get_RHS_coefs_y(x, y);
            }

            bc = m.get_y_last_coefs(x);
            mx_y[0][y_dim - 1] = bc[0];
            mx_y[1][y_dim - 1] = bc[1];
            mx_y[3][y_dim - 1] = m.get_RHS_coefs_y(x, y_dim - 1);

            // call solver and update current values in the column
            // std::cout << "matrix " << x << "\n";
            // std::getchar();
            solver_y.solve(mx_y, f_y);
            if (VERBOSE) {
                pprint_tridiag_matrix(mx_y, std::cout);
                pprint_solution_row(f_y, std::cout);
            }
            update_grid_col(x);
        }
    }

    void Problem::update_grid_row(const size_t y) {
        size_t i = 0;
        std::for_each(
            f_x.cbegin(), f_x.cend(),
            [&](const double & f) { m.set_current_value(i++, y, f); }
        );
    }

    void Problem::update_grid_col(const size_t x) {
        size_t i = 0;
        std::for_each(
            f_y.cbegin(), f_y.cend(),
            [&](const double & f) { m.set_current_value(x, i++, f); }
        );
    }
}
