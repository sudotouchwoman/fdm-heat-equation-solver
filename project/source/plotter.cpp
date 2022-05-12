#include "plotter.hpp"

namespace plt {
    GNUPlotWriter::GNUPlotWriter(const std::string & config) {
        pipe = popen(command.data(), "w");
        if (pipe == nullptr) throw std::runtime_error("failed to run GNUplot");
        fputs(config.c_str(), pipe);
    }

    GNUPlotWriter::~GNUPlotWriter() {
        pclose(pipe);
    }

    void GNUPlotWriter::flush_buffer() {
        throw_on_bad_pipe();
        fputs("splot '-' matrix with image\n", pipe);
        fputs(payload_buffer.str().c_str(), pipe);
        fputs("e\n", pipe);
        std::stringstream().swap(payload_buffer);
    }

    void GNUPlotWriter::throw_on_bad_pipe() {
        if (pipe == nullptr)
            throw std::runtime_error("failed to open GNUPlot subprocess");
    }
}
