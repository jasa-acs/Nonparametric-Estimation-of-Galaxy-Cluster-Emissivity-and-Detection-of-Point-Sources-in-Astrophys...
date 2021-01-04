///
/// \file src/main.cpp
/// \brief Launcher for the ASTROQUT solver.
/// \details Handle the user input, calls the preparation tools and the solver.
/// \author Philippe Ganz <philippe.ganz@gmail.com> 2017-2019
/// \version 1.0.1
/// \date August 2019
/// \copyright GPL-3.0
///

#include "const.hpp"
#include "WS/astroQUT.hpp"

#include <cstdlib>
#include <getopt.h>
#include <iostream>
#include <string>

void usage()
{
    std::cerr << std::endl << "usage : ASTROQUT -f|--source SOURCE -e|--sensitivity SENSITIVITY -o|--background BACKGROUND -b|--blurring BLURRING -r|--result RESULT -s|--size SIZE -x|--bootstrap BOOTSTRAP" << std::endl << std::endl;
    std::cerr << "  SOURCE - Path to the source image;" << std::endl;
    std::cerr << "  SENSITIVITY - Path to the sensitivity image;" << std::endl;
    std::cerr << "  BACKGROUND - Path to the background image;" << std::endl;
    std::cerr << "  BLURRING - Path to the blurring filter, defaults to data/blurring.data;" << std::endl;
    std::cerr << "  RESULT - Path to the solution file;" << std::endl;
    std::cerr << "  SIZE - Width of the picture;" << std::endl;
    std::cerr << "  BOOTSTRAP - Amount of bootstraps to perform." << std::endl << std::endl;
}

int main( int argc, char **argv )
{
    alias::WS::Parameters<double> options;
    std::string source;
    std::string sensitivity;
    std::string background;
    std::string blurring;
    std::string result;
    size_t pic_size = 0;
    size_t bootstrap_max = 0;

    int c;

    while (1)
    {
        static struct option long_options[] =
        {
            {"source",      required_argument, nullptr, 'f'},
            {"sensitivity", required_argument, nullptr, 'e'},
            {"background",  required_argument, nullptr, 'o'},
            {"blurring",    required_argument, nullptr, 'b'},
            {"result",      required_argument, nullptr, 'r'},
            {"size",        required_argument, nullptr, 's'},
            {"bootstrap",   required_argument, nullptr, 'x'},
            {nullptr,       0,                 nullptr, 0}
        };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "f:e:o:b:r:s:x:", long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c)
        {
        case 'f':
        {
            source = std::string(optarg);
            break;
        }

        case 'e':
        {
            sensitivity = std::string(optarg);
            break;
        }

        case 'o':
        {
            background = std::string(optarg);
            break;
        }

        case 'b':
        {
            blurring = std::string(optarg);
            break;
        }

        case 'r':
        {
            result = std::string(optarg);
            break;
        }

        case 's':
        {
            pic_size = strtoul(optarg, nullptr, 0);
            break;
        }

        case 'x':
        {
            bootstrap_max = strtoul(optarg, nullptr, 0);
            break;
        }

        default:
        {
            usage();
            return EXIT_FAILURE;
        }
        }
    }

    if( source.compare("") == 0 ||
        sensitivity.compare("") == 0 ||
        background.compare("") == 0 ||
        blurring.compare("") == 0 ||
        result.compare("") == 0 ||
        pic_size == 0 ||
        bootstrap_max == 0)
    {
        usage();
        return EXIT_FAILURE;
    }

    options.blurring_filter = blurring;
    options.pic_size = pic_size;
    options.bootstrap_max = bootstrap_max;

    alias::WS::Solve(source, sensitivity, background, result, options);

    return EXIT_SUCCESS;
}
