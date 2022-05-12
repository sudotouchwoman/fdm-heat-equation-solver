#include "model.hpp"

#include <iostream>
#include <iomanip>

/*
*    in schematic, the geometry of the plate
*    can be drawn like this (god forgive me)
*
*   # - no condition
*   * - 1st type condition
*   + - 2nd type condition
*   @ - 3rd type condition
*
*   * * * * * * * * * * * * * * * * * * * * * * * * * *                                                 
*   + # # # # # # # # # # # # # # # # # # # # # # # # # *                                               
*   + # # # # # # # # # # # # # # # # # # # # # # # # # # *                                             
*   + # # # # # # # # # # # # # # # # # # # # # # # # # # # *                                           
*   + # # # # # # # # # # # # # # # # # # # # # # # # # # # # *                                         
*   + # # # # # # # # # # # # # # # # # # # # # # # # # # # # # *                                       
*   + # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # *                                     
*   + # # # # # # # # # # @ @ @ @ @ @ @ @ @ @ @ @ @ @ # # # # # # # *                                   
*   + # # # # # # # # # @                             @ # # # # # # # *                                 
*   + # # # # # # # # # @                             @ # # # # # # # # *                               
*   + # # # # # # # # # @                             @ # # # # # # # # # *                             
*   + # # # # # # # # # @                             @ # # # # # # # # # # *                           
*   + # # # # # # # # # @                             @ # # # # # # # # # # # *                         
*   + # # # # # # # # # @                             @ # # # # # # # # # # # # *                       
*   + # # # # # # # # # @                             @ # # # # # # # # # # # # # *                     
*   + # # # # # # # # # @                             @ # # # # # # # # # # # # # # *                   
*   + # # # # # # # # # @                             @ # # # # # # # # # # # # # # # *                 
*   + # # # # # # # # # # @ @ @ @ @ @ @ @ @ @ @ @ @ @ # # # # # # # # # # # # # # # # # *               
*   + # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # *             
*   + # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # *           
*   + # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # *         
*   + # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # *       
*   + # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # *     
*   + # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # *   
*   * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
*/

namespace model {
    void Model79::dump(std::ostream & os) const {
        for (const auto & row: grid.nodes) {
            for (const auto & e: row) {
                os << e.current_value << ' ';
            }
            os << '\n';
        }
    }

    void Model79::throw_on_bounds(const size_t x, const size_t y) const {
        if (x >= dims.first) throw std::runtime_error("X index exceeding grid bounds");
        if (y >= dims.second) throw std::runtime_error("Y index exceeding grid bounds");
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
            grid.nodes[i][0].condition_type = BOUNDARY_2TYPE_X;
        }

        // floor
        for (size_t i = 0; i < x_dim; ++i) {
            grid.nodes[y_dim - 1][i].condition_type = BOUNDARY_1TYPE;
            grid.nodes[y_dim - 1][i].current_value = T_FLOOR;
            grid.nodes[y_dim - 1][i].initial_value = T_FLOOR;
        }

        // ceiling
        for (size_t i = 0; i < x_half; ++i) {
            grid.nodes[0][i].condition_type = BOUNDARY_1TYPE;
            grid.nodes[0][i].current_value = T_CEIL;
            grid.nodes[0][i].initial_value = T_CEIL;
        }

        // right side
        for (size_t i = 0; i < x_half; ++i) {
            // nodes to the right to the inclined side
            for (size_t j = i; j < y_dim; ++j) {
                grid.nodes[y_dim - j - 1][x_dim - i - 1].condition_type = OUTER_NODE;
            }
            grid.nodes[x_half - i - 1][x_dim - i - 1].condition_type = BOUNDARY_1TYPE;
            grid.nodes[x_half - i - 1][x_dim - i - 1].current_value = T_CEIL;
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
            grid.nodes[y_hole_lower][i].condition_type = BOUNDARY_3TYPE_Y;
            grid.nodes[y_hole_upper][i].condition_type = BOUNDARY_3TYPE_Y;
        }
        for (size_t i = y_hole_lower + 1; i < y_hole_upper; ++i) {
            grid.nodes[i][x_hole_left].condition_type = BOUNDARY_3TYPE_X;
            grid.nodes[i][x_hole_right].condition_type = BOUNDARY_3TYPE_X;
        }

        // corner nodes
        grid.nodes[y_hole_upper][x_hole_left].condition_type = BOUNDARY_3TYPE_XY;
        grid.nodes[y_hole_upper][x_hole_right].condition_type = BOUNDARY_3TYPE_XY;
        grid.nodes[y_hole_lower][x_hole_right].condition_type = BOUNDARY_3TYPE_XY;
        grid.nodes[y_hole_lower][x_hole_left].condition_type = BOUNDARY_3TYPE_XY;
    }

    bool Model79::is_inner(const size_t x, const size_t y) const {
        throw_on_bounds(x, y);
        return not (grid.nodes[y][x].condition_type == OUTER_NODE);
    }

    void Model79::set_current_value(const size_t x, const size_t y, const double value) {
        throw_on_bounds(x, y);
        const condition cond = grid.nodes[y][x].condition_type;
        if (cond == BOUNDARY_1TYPE or cond == OUTER_NODE) return;
        grid.nodes[y][x].current_value = value;
    }

    double Model79::get_RHS_coefs_x(const size_t x, const size_t y) const {
        throw_on_bounds(x, y);
        const Node grid_node = grid.nodes[y][x];

        switch (grid_node.condition_type) {
            case OUTER_NODE:
                return grid_node.initial_value;
            case BOUNDARY_3TYPE_XY:
            case BOUNDARY_2TYPE_X:
            case BOUNDARY_3TYPE_X:
                return 0;
            case BOUNDARY_2TYPE_Y:
            case BOUNDARY_3TYPE_Y:
            case BOUNDARY_1TYPE:
            case NO_CONDITION:
                return grid_node.current_value;
            default:
                throw std::runtime_error("unknown condition type");
        }
    }

    double Model79::get_RHS_coefs_y(const size_t x, const size_t y) const {
        throw_on_bounds(x, y);
        const Node grid_node = grid.nodes[y][x];

        switch (grid_node.condition_type) {
            case OUTER_NODE:
                return grid_node.initial_value;
            case BOUNDARY_3TYPE_XY:
            case BOUNDARY_2TYPE_Y:
            case BOUNDARY_3TYPE_Y:
                return 0;
            case BOUNDARY_2TYPE_X:
            case BOUNDARY_3TYPE_X:
            case BOUNDARY_1TYPE:
            case NO_CONDITION:
                return grid_node.current_value;
            default:
                throw std::runtime_error("unknown condition type");
        }
    }

    tridiag_coefs Model79::get_x_coefs(const size_t x, const size_t y) const {
        throw_on_bounds(x, y);
        const double R = a * dt / (dx * dx);    // needed for nodes with no boundary
        const condition cond = grid.nodes[y][x].condition_type;

        switch (cond) {
            // must be constant value
            case OUTER_NODE:
            case BOUNDARY_1TYPE:
                return {0, 1.0, 0};
            case BOUNDARY_3TYPE_XY:
            case BOUNDARY_3TYPE_X: {
                if (grid.nodes[y][x + 1].condition_type == NO_CONDITION) {
                    const double c = (-1.0) / (1.0 + dx);
                    return {0, 1.0, c};
                }
                if (grid.nodes[y][x - 1].condition_type == NO_CONDITION) {
                    const double c = (-1.0) / (1.0 + dx);
                    return {c, 1.0, 0};
                }
            }
            case BOUNDARY_3TYPE_Y:
            case BOUNDARY_2TYPE_Y:
            case NO_CONDITION: {
                return {-R, 2.0*R + 1.0, -R};
            }
            default:
                throw std::runtime_error("unknown condition type");
        }
    }

    tridiag_coefs Model79::get_y_coefs(const size_t x, const size_t y) const {
        throw_on_bounds(x, y);
        const double R = a * dt / (dy * dy);
        const condition cond = grid.nodes[y][x].condition_type;

        switch (cond) {
            // must be constant value
            case OUTER_NODE:
            case BOUNDARY_1TYPE:
                return {0, 1.0, 0};
            case BOUNDARY_3TYPE_XY:
            case BOUNDARY_3TYPE_Y: {
                // in this particular case, 3rd type boundaries
                // are only defined at inner nodes -> if one has these at the edge
                // of the mesh, it is a good idea to bound-check the axis first
                // in order to avoid segfaults
                if (grid.nodes[y + 1][x].condition_type == NO_CONDITION) {
                    const double c = (-1.0) / (1.0 + dx);
                    return {0, 1.0, c};
                }
                if (grid.nodes[y - 1][x].condition_type == NO_CONDITION) {
                    const double c = (-1.0) / (1.0 + dx);
                    return {c, 1.0, 0};
                }
            }
            // in this particular problem, 2nd type boundaries
            // only appears on the edge
            // in X-direction thus are treated as ordinary inner nodes
            case BOUNDARY_3TYPE_X:
            case BOUNDARY_2TYPE_X:
            case NO_CONDITION: {
                return {-R, 2.0*R + 1.0, -R};
            }
            default:
                throw std::runtime_error("unknown condition type");
        }
    }

    // numeration order is as follows:
    // upmost nodes = 0, lowest nodes = x_max, leftmost nodes = 0, rightmost_nodes = y_max
    // these functions would only apply to my particular problem
    boundary_coefs Model79::get_x_last_coefs(const size_t x) const {
        throw_on_bounds(x, 0);
        return {0, 1.0};
    }

    boundary_coefs Model79::get_x_first_coefs(const size_t y) const {
        throw_on_bounds(0, y);
        // in this case, the 2 TYPE boundary is defined on the
        // leftmost side of the plate  
        return {-1.0, 1.0};
    }

    boundary_coefs Model79::get_y_last_coefs(const size_t x) const {
        throw_on_bounds(x, 0);
        return {0.0, 1.0};
    }

    boundary_coefs Model79::get_y_first_coefs(const size_t x) const {
        throw_on_bounds(x, 0);
        return {1.0, 0.0};
    }

    static char cond_to_symbol(const condition c) {
        switch (c) {
            case NO_CONDITION:
            return '#';
            case BOUNDARY_1TYPE:
            return '*';
            case BOUNDARY_2TYPE_X:
            case BOUNDARY_2TYPE_Y:
            return '@';
            case BOUNDARY_3TYPE_X:
            return '-';
            case BOUNDARY_3TYPE_Y:
            return '|';
            case BOUNDARY_3TYPE_XY:
            return 'x';
            case OUTER_NODE:
            return '.';
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
