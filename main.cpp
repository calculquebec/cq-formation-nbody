#include <fstream>
#include <iostream>

#include "nbody.h"


Params read_parameters(const char * filename)
{
    std::ifstream s(filename);

    if (!s.is_open()) {
        // If the file doesn't exist, we need to exit...
        std::cout << "The file " << filename << " cannot be found!" << std::endl;
        std::exit(1);
    }

    Params params;
    std::string line, key, value;

    // Loop through all lines in the parameter file
    while(std::getline(s, line)) {
        // If it's an empty line, continue
        if (line.empty()) continue;

        // If the line begins with a #, ignore it
        if (line[0] == '#') continue;

        // Assumes that the equals sign can only occur once in the line
        size_t equalPos = line.find('=');

        // If there's no equals sign in this line, continue
        if (equalPos == std::string::npos) continue;

        key = "";
        for (size_t i = 0; i < equalPos; ++i) {
            if (line[i] == ' ') continue;
            key += line[i];
        }

        value = "";
        for (size_t i = equalPos + 1; i < line.size(); ++i) {
            if (line[i] == ' ') continue;
            value += line[i];
        }

        // Save the key,value pair
        params[key] = value;
    }

    s.close();

    return params;
}


void debug_parameters(const Params &params)
{
    for (Params::const_iterator it = params.begin(); it != params.end(); it++) {
        std::cout << it->first << " = " << it->second << std::endl;
    }
}


int main(int argc, char** argv)
{
    if (argc != 2) {
        std::cerr << "Usage: ./nbody parameters.txt" << std::endl;
        return 1;
    }

    Params params = read_parameters(argv[1]);
    //debug_parameters(params);

#ifdef VERLET
    NBody model(NBody::Verlet);
#else
    NBody model(NBody::Runge_Kutta);
#endif

    //model.printConfig();
    model.configure(params);
    //model.printConfig();

    model.integrate();

    return 0;
}

