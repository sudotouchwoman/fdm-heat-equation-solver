#include <iostream>
#include <memory>

#include "solver.hpp"
#include "plotter.hpp"

constexpr size_t TIMESTEPS = 1000;
constexpr double TIME = 15.0;

int main(int argc, char* argv[]) {
    const double dt = TIME / static_cast<double>(TIMESTEPS);
    const double a = 10.0;
    const double x_len = 10.0;
    const double y_len = 5.0;
    const size_t x_nodes = 200;
    const size_t y_nodes = 100;

    const double dx = x_len / static_cast<double>(x_nodes);
    const double dy = y_len / static_cast<double>(y_nodes);

    model::Model79 m(dt, dx, dy, a, x_nodes, y_nodes);
    solver::Problem problem(m, TIMESTEPS);
    plt::GNUPlotWriter plotter(plt::GNUPlotWriter::basic_config.data());

    pprint_grid(m, std::cout);
    // plotter.reciever() << m;
    // plotter.flush_buffer();
    for (size_t i = 0; i < TIMESTEPS; ++i) {
        problem.step();
        // plotter.reciever() << m;
        // plotter.flush_buffer();
    }

    plotter.reciever() << m;
    plotter.flush_buffer();
    return 0;
}