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
        constexpr static std::string_view basic_config = 
        "set terminal gif size 800 800 animate delay 2 enhanced font"
        "'Verdana, 14'\n"
        "set output 'map.gif'\n"
        "set title 'Heat equation solution using TDMA solver'\n"
        "set xlabel 'X'\n"
        "set ylabel 'Y'\n"
        "set xrange [0:199]\n"
        "set yrange [0:99]\n"
        "set view map scale 1\n"
        // "set size ratio 1\n"
        "set palette color\n"
        // "set pm3d interpolate 1,1\n"
        "set pm3d map\n";
        // "set dgrid3d\n"
        // "set colorbox vertical origin screen 0.0, 0.0"
        // "size screen 0.0, 0.0 front noinvert noborder\n";
    private:
        FILE * pipe = nullptr;
        std::stringstream payload_buffer;
        constexpr static std::string_view command = "gnuplot";
    private:
        void throw_on_bad_pipe();
    };
}
