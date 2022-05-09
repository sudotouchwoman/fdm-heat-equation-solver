#include "model.hpp"

#include <iostream>

/*
*    in schematic, the geometry of the plate
*    can be drawn like this (god forgive me)
*
*    +---------+
*    |           \
*    |   +-----+  \
*    |   |     |   \
*    |   |     |    \
*    |   +-----+     \
*    |                \
*    +-----------------+
*/

namespace model {
    void Model79::throw_on_bounds(const size_t x, const size_t y) const {
        if (x > dims.first) throw std::runtime_error("X index exceding grid bounds");
        if (y > dims.second) throw std::runtime_error("Y index exceding grid bounds");
    }

    void Model79::grid_set_up() {
        const size_t x_dim = dims.first;
        const size_t y_dim = dims.second;

        // trust me these define the relative
        // positions of key points
        const size_t x_half = x_dim / 2;
        const size_t x_hole_left = x_half * 2 / 5;
        const size_t x_hole_right = x_half;
        const size_t y_hole_lower = y_dim * 1.5 / 5;
        const size_t y_hole_upper = y_dim * 3.5 / 5;

        // left side
        for (size_t i = 0; i < y_dim; ++i) {
            grid.nodes[i][0].condition_type = BOUNDARY_2TYPE;
        }

        // floor
        for (size_t i = 0; i < x_dim; ++i) {
            grid.nodes[y_dim - 1][i].condition_type = BOUNDARY_1TYPE;
            grid.nodes[y_dim - 1][i].initial_value = T_FLOOR;
        }

        // ceiling
        for (size_t i = 0; i < x_half; ++i) {
            grid.nodes[0][i].condition_type = BOUNDARY_1TYPE;
            grid.nodes[0][i].initial_value = T_CEIL;
        }

        // right side
        for (size_t i = 0; i < x_half; ++i) {
            // nodes to the right to the inclined side
            for (size_t j = i; j < y_dim; ++j) {
                grid.nodes[y_dim - j - 1][x_dim - i].condition_type = OUTER_NODE;
            }
            grid.nodes[x_half - i - 1][x_dim - i - 1].condition_type = BOUNDARY_1TYPE;
            grid.nodes[x_half - i - 1][x_dim - i - 1].initial_value = T_CEIL;
        }

        // hole
        for (size_t i = x_hole_left + 1; i < x_hole_right; ++i) {
            for (size_t j = y_hole_lower + 1; j < y_hole_upper; ++j) {
                // nodes inside the hole are marked as outer
                grid.nodes[j][i].condition_type = OUTER_NODE;
            }
        }
        for (size_t i = x_hole_left + 1; i < x_hole_right; ++i) {
            grid.nodes[y_hole_lower][i].condition_type = BOUNDARY_3TYPE;
            grid.nodes[y_hole_upper][i].condition_type = BOUNDARY_3TYPE;
        }
        for (size_t i = y_hole_lower + 1; i < y_hole_upper; ++i) {
            grid.nodes[i][x_hole_left].condition_type = BOUNDARY_3TYPE;
            grid.nodes[i][x_hole_right].condition_type = BOUNDARY_3TYPE;
        }
    }

    bool Model79::is_inner(const size_t x, const size_t y) const {
        throw_on_bounds(x, y);
        return not (grid.nodes[y][x].condition_type == OUTER_NODE);
    }

    tridiag_coefs Model79::get_x_coefs(const size_t x, const size_t y) const {
        tridiag_coefs coefs = {0, 0, 0};
        const double c = (-1.0) / (dx + 1);
        const double R = a * dt / (dx * dx);

        const condition cond = grid.nodes[y][x].condition_type;
        const condition cond_up = grid.nodes[y + 1][x].condition_type;

        switch (cond) {
            case OUTER_NODE:
            case BOUNDARY_1TYPE:
                // must be constant value
                coefs[1] = 1;
                return coefs;
            // case BOUNDARY_2TYPE: {
            //     // if node above is inner (has no condition applied)
            //     if (cond_up == NO_CONDITION) {
            //         coefs = {0, 1, -1};
            //     } else {
            //         // otherwise, node below must be inner
            //         coefs = {-1, 1, 0};
            //     }
            //     return coefs;
            // }
            case BOUNDARY_3TYPE: {
                // same story as with the 2 type condition
                if (cond_up == NO_CONDITION) {
                    coefs = {0, 1, c};
                } else {
                    coefs = {c, 0, 1};
                }
                return coefs;
            }
            case NO_CONDITION: {
                coefs = {-R, R + 1, -R};
                return coefs;
            }
            default:
                throw std::runtime_error("unknown condition type");
        }
    }

    tridiag_coefs Model79::get_y_coefs(const size_t x, const size_t y) const {
        tridiag_coefs coefs = {0, 0, 0};
        const double c = (-1.0) / (dx + 1);
        const double R = a * dt / (dy * dy);

        const condition cond = grid.nodes[y][x].condition_type;
        const condition cond_right = grid.nodes[y][x + 1].condition_type;

        switch (cond) {
            case OUTER_NODE:
            case BOUNDARY_1TYPE:
                // must be constant value
                coefs[1] = 1;
                return coefs;
            // case BOUNDARY_2TYPE: {
            //     // if node above is inner (has no condition applied)
            //     if (cond_right == NO_CONDITION) {
            //         coefs = {0, 1, -1};
            //     } else {
            //         // otherwise, node below must be inner
            //         coefs = {-1, 1, 0};
            //     }
            //     return coefs;
            // }
            case BOUNDARY_3TYPE: {
                // same story as with the 2 type condition
                if (cond_right == NO_CONDITION) {
                    coefs = {0, 1, c};
                } else {
                    coefs = {c, 0, 1};
                }
                return coefs;
            }
            case NO_CONDITION: {
                coefs = {-R, R + 1, -R};
                return coefs;
            }
            default:
                throw std::runtime_error("unknown condition type");
        }
    }

    // numeration order is as follows:
    // upmost nodes = 0, lowest nodes = x_max, leftmost nodes = 0, rightmost_nodes = 0
    // these functions would only apply to my particuler problem
    boundary_coefs Model79::get_x_upper_coefs(const size_t y) const {
        boundary_coefs coefs = {0, 1};
        return coefs;
    }

    boundary_coefs Model79::get_x_lower_coefs(const size_t x) const {
        // in this case, the 2 TYPE boundary is defined on the
        // leftmost side of the plate  
        boundary_coefs coefs = {1, -1};
        return coefs;
    }

    boundary_coefs Model79::get_y_upper_coefs(const size_t y) const {
        boundary_coefs coefs = {0, 1};
        return coefs;
    }

    boundary_coefs Model79::get_y_lower_coefs(const size_t x) const {
        boundary_coefs coefs = {1, 0};
        return coefs;
    }

    static char cond_to_symbol(const condition c) {
        switch (c) {
            case NO_CONDITION:
            return '#';
            case BOUNDARY_1TYPE:
            return '*';
            case BOUNDARY_2TYPE:
            return '+';
            case BOUNDARY_3TYPE:
            return '@';
            case OUTER_NODE:
            return ' ';
            default:
            throw std::runtime_error("unknown condition type");
        }
    }

    void pprint_grid(const Model79 & m, std::ostream & out) {
        const size_t x_dim = m.x_dim();
        const size_t y_dim = m.y_dim();

        for (size_t i = 0; i < y_dim; ++i) {
            for (size_t j = 0; j < x_dim; ++j) {
                out << cond_to_symbol(m.grid.nodes[i][j].condition_type) << ' ';
            }
            out << '\n';
        }
    }
}
