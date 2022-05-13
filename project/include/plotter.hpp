#pragma once

#include <cctype>
#include <sstream>

#include "shared.hpp"

namespace plt {
    class GNUPlotWriter {
    public:
        GNUPlotWriter(const std::string & config);
        ~GNUPlotWriter();
        std::ostream & reciever() { return payload_buffer; }
        void flush_buffer();
        constexpr static std::string_view basic_gif_config = 
            "set terminal gif size 800 800 animate delay 2 enhanced font"
            "'Verdana, 14'\n"
            "set output 'map.gif'\n"
            "set title 'Heat equation solution using TDMA solver'\n"
            "set xlabel 'X'\n"
            "set ylabel 'Y'\n"
            // "set xrange [0:199]\n"
            // "set yrange [0:99]\n"
            "set view map scale 1\n"
            "set palette color\n"
            "set pm3d map\n";
    private:
        FILE * pipe = nullptr;
        std::stringstream payload_buffer;
        // explicitly prohibit writing anything to the terminal
        constexpr static std::string_view command = "gnuplot 2> /dev/null";
    private:
        void throw_on_bad_pipe();
    };
}
