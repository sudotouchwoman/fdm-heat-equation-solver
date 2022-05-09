#pragma once

#include <array>
#include <vector>

#define T_FLOOR 50.0
#define T_CEIL 80.0

namespace model {
    using tridiag_coefs = std::array<double, 3>;
    using boundary_coefs = std::array<double, 2>;

    enum condition {
        NO_CONDITION,
        BOUNDARY_1TYPE,
        BOUNDARY_2TYPE,
        BOUNDARY_3TYPE,
        OUTER_NODE
    };

    enum direction {
        VERTICAL,
        HORIZONTAL
    };

    // model instance should provide sets of coefficients
    // for each grid point in horizontal and vertical directions
    // to apply the Thompson's tridiagonal matrix algorithm
    // additionally, model can store required parameters like dt, dx, dy
    // and a_1, a_2
    // the heat PDE is defined as follows:
    // \frac{\partial{T}^2}{\partial{t}^2} = a_1 * \frac{\partial^2{T}}{\partial^2{x}} + a_2 * \frac{\partial^2{T}}{\partial^2{y}}
    class IModel {
    public:
        virtual tridiag_coefs get_x_coefs(const size_t x, const size_t y) const = 0;
        virtual tridiag_coefs get_y_coefs(const size_t x, const size_t y) const = 0;
        virtual boundary_coefs get_x_upper_coefs(const size_t y) const = 0;
        virtual boundary_coefs get_y_upper_coefs(const size_t y) const = 0;
        virtual boundary_coefs get_x_lower_coefs(const size_t x) const = 0;
        virtual boundary_coefs get_y_lower_coefs(const size_t x) const = 0;
        virtual size_t x_dim() const = 0;
        virtual size_t y_dim() const = 0;
    };

    // Model79 implements IModel interface and stands for my particular problem setup
    // thus such methods as is_inner and is_border are present to deduce
    // the geometry. This is not the most elegant approach, however...
    class Model79: public IModel {
    private:
        // holds boundary condition type if present
        // and the default temperature value in case of
        // Dirichlet's BC (T = const)
        struct Node {
            condition condition_type = NO_CONDITION;
            double current_value = 0;
            double initial_value = 0;
        };
        // Grid entity contains matrix of condition flags
        // for ease boundary condition detection
        struct Grid {
            const size_t width;
            const size_t height;
            std::vector<std::vector<Node>> nodes;
            Grid(const size_t width, const size_t height):
                width(width), height(height),
                nodes(height, std::vector<Node>(width)) {};
            ~Grid() = default;
        };
    private:
        const double dt;
        const double dx;
        const double dy;
        const double a;
        const std::pair<size_t, size_t> dims;
        Grid grid;
    private:
        void throw_on_bounds(const size_t x, const size_t y) const;
        void grid_set_up();  // init grid with required flags + default values
        bool is_inner(const size_t x, const size_t y) const;
    public:
        tridiag_coefs get_x_coefs(const size_t x, const size_t y) const override;
        tridiag_coefs get_y_coefs(const size_t x, const size_t y) const override;
        boundary_coefs get_x_upper_coefs(const size_t y) const override;
        boundary_coefs get_y_upper_coefs(const size_t y) const override;
        boundary_coefs get_x_lower_coefs(const size_t x) const override;
        boundary_coefs get_y_lower_coefs(const size_t x) const override;
        size_t x_dim() const override { return dims.first; }
        size_t y_dim() const override { return dims.second; }

        ~Model79() = default;
        Model79() = delete;
        Model79(
            const double dt,
            const double dx,
            const double dy,
            const double a,
            const size_t x_nodes,
            const size_t y_nodes) :
            dt(dt), dx(dx), dy(dy), a(a),
            dims(std::make_pair(x_nodes, y_nodes)),
            grid(x_nodes, y_nodes) { grid_set_up(); };

        friend void pprint_grid(const Model79 & m, std::ostream & out);
    };
}
