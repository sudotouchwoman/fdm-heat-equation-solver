#include <iostream>
#include <memory>

#include "solver.hpp"
#include "plotter.hpp"

constexpr size_t DEF_TIMESTEPS = 1000;
constexpr double DEF_TIME = 15.0;
constexpr size_t X_NODES = 200;
constexpr size_t Y_NODES = 100;
constexpr double X_LEN = 10.0;
constexpr double Y_LEN = 5.0;

constexpr std::string_view running = "Performing computations: ";
constexpr std::string_view usage = "Usage: <simulation time:double> [<timesteps:uint> [<x_nodes:uint> <y_nodes:uint>]]\n";

int main(int argc, char* argv[]) {
    double time = DEF_TIME;
    double timesteps = DEF_TIMESTEPS;

    const double a = 1.0;

    size_t x_nodes = X_NODES;
    size_t y_nodes = Y_NODES;

    // not the most versatile solution, however
    // it is OK for this case
    if (argc == 1) {
        std::cout << usage;
        return EXIT_FAILURE;
    }
    if (argc >= 2) {
        time = std::stod(argv[1]);
    }
    if (argc >= 3) {
        timesteps = std::stoul(argv[2]);
    }
    if (argc == 5) {
        x_nodes = std::stoul(argv[3]);
        y_nodes = std::stoul(argv[4]);

        if (x_nodes != y_nodes * 2) {
            std::cerr << "X to Y ratio must be 2:1\n";
            return EXIT_FAILURE;
        }
    }

    std::cout << "Simulation time set to " << time << '\n';
    std::cout << "Timesteps set to " << timesteps << '\n';
    std::cout << "Mesh size: [" << x_nodes << ':' << y_nodes << "]\n";

    const double dt = time / static_cast<double>(timesteps);
    const double dx = X_LEN / static_cast<double>(x_nodes);
    const double dy = Y_LEN / static_cast<double>(y_nodes);

    // instantiate model for my case, set up problem
    // environment (e.g., allocate memory for solvers)
    // and gnuplot wrapper to create heatmap gif
    model::Model79 m(dt, dx, dy, a, x_nodes, y_nodes);
    solver::Problem problem(m, timesteps);
    plt::GNUPlotWriter plotter(plt::GNUPlotWriter::basic_gif_config.data());

    std::cout << "The problem schematic (may not fit into the terminal entirely)\n";
    pprint_grid(m, std::cout);
    plotter.reciever() << m;
    plotter.flush_buffer();

    std::cout << running;
    for (size_t i = 0; i < timesteps; ++i) {
        std::cout << '\r' << running << "iter [" << i << '/' << timesteps << ']';
        std::cout.flush();

        problem.step();
        plotter.reciever() << m;
        plotter.flush_buffer();
    }

    std::cout << " Done, OK\n";
    return EXIT_SUCCESS;
}