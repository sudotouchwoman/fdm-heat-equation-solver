#include <iostream>

#include "model.hpp"

const double TIMESTEPS = 1000;
const double TIME = 15;

int main(int argc, char* argv[]) {
    using namespace model;

    const double dt = (double)TIME / TIMESTEPS;
    const double dx = 0.1;
    const double dy = 0.1;
    const double a = 1.0;
    const size_t x_nodes = 50;
    const size_t y_nodes = 25;

    Model79 model(dt, dx, dy, a, x_nodes, y_nodes);
    pprint_grid(model, std::cout);
    return 0;
}