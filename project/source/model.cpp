#include "model.hpp"

struct Point {
    const double x;
    const double y;
    Point(const double x, const double y) :
        x(x), y(y) {};
    ~Point() = default;
};

// main geometry points
const Point LOWER_LEFT(0.0, 0.0);
const Point LOWER_RIGHT(10.0, 0.0);
const Point UPPER_LEFT(0.0, 5.0);
const Point UPPER_RIGHT(5.0, 5.0);

const Point HOLE_LOWER_LEFT(2.0, 1.5);
const Point HOLE_LOWER_RIGHT(5.0, 1.5);
const Point HOLE_UPPER_LEFT(2.0, 3.5);
const Point HOLE_UPPER_RIGHT(5.0, 3.5);

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

    static bool between(const double i, const double a, const double b) {
        return a < i and i < b;
    }

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
            grid.nodes[0][i].condition_type = BOUNDARY_2TYPE;
        }

        // floor
        for (size_t i = 0; i < x_dim; ++i) {
            grid.nodes[i][0].condition_type = BOUNDARY_1TYPE;
            grid.nodes[i][0].initial_value = T_FLOOR;
        }

        // ceiling
        for (size_t i = 0; i < x_half; ++i) {
            grid.nodes[i][y_dim - 1].condition_type = BOUNDARY_1TYPE;
            grid.nodes[i][y_dim - 1].initial_value = T_CEIL;
        }

        // right side
        for (size_t i = 0; i < x_half; ++i) {
            // nodes to the right to the inclined side
            for (size_t j = i; j < y_dim; ++j) {
                grid.nodes[x_half + i][j].condition_type = OUTER_NODE;
            }
            grid.nodes[x_half + i][x_half - i - 1].condition_type = BOUNDARY_1TYPE;
            grid.nodes[x_half + i][x_half - i - 1].initial_value = T_CEIL;
        }

        // hole
        for (size_t i = x_hole_left; i < x_hole_right; ++i) {
            for (size_t j = y_hole_lower; j < y_hole_upper - 1; ++j) {
                // nodes inside the hole are marked as outer
                grid.nodes[i][j].condition_type = OUTER_NODE;
            }
            grid.nodes[i][y_hole_lower].condition_type = BOUNDARY_3TYPE;
            grid.nodes[i][y_hole_upper].condition_type = BOUNDARY_3TYPE;
        }
        for (size_t i = y_hole_lower; i < y_hole_upper; ++i) {
            grid.nodes[x_hole_left][i].condition_type = BOUNDARY_3TYPE;
            grid.nodes[x_hole_right][i].condition_type = BOUNDARY_3TYPE;
        }
    }

    bool Model79::is_inner(const size_t x_node, const size_t y_node) const {
        throw_on_bounds(x_node, y_node);
        const double x = LOWER_LEFT.x + dx * x_node;
        const double y = LOWER_LEFT.y + dy * y_node;

        // must be out there
        if (not between(x, LOWER_LEFT.x, LOWER_RIGHT.x)) return false;
        if (not between(y, LOWER_LEFT.y, UPPER_LEFT.y)) return false;

        // above the right inclined side
        if (x > UPPER_RIGHT.x and x - UPPER_RIGHT.x > y) return false;

        // inside the square hole
        if (
            between(x, HOLE_LOWER_LEFT.x, HOLE_UPPER_RIGHT.x) and
            between(y, HOLE_LOWER_LEFT.y, HOLE_UPPER_RIGHT.y)
        ) return false;

        return true;
    }

    tridiag_coefs Model79::get_x_coefs(const size_t x, const size_t y) const {
        tridiag_coefs coefs = {0, 0, 0};
        const double c = (-1.0) / (dx + 1);

        const condition cond = grid.nodes[x][y].condition_type;
        const condition cond_up = grid.nodes[x][y + 1].condition_type;

        switch (cond) {
            case NO_CONDITION:
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
            default:
                throw std::runtime_error("unknown condition type");
        }
    }

    tridiag_coefs Model79::get_y_coefs(const size_t x, const size_t y) const {
        tridiag_coefs coefs = {0, 0, 0};
        const double c = (-1.0) / (dx + 1);

        const condition cond = grid.nodes[x][y].condition_type;
        const condition cond_right = grid.nodes[x + 1][y].condition_type;

        switch (cond) {
            case NO_CONDITION:
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
}
