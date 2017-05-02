#include <iostream>

#include "nbody.h"


Params read_parameters(const char * filename)
{
    Params params;

    return params;
}


int main(int argc, char** argv)
{
    if (argc != 2) {
        std::cerr << "Usage: ./nbody parameters.txt" << std::endl;
        return 1;
    }

    Params params = read_parameters(argv[1]);

#ifdef VERLET
    NBody model(NBody::Verlet);
#else
    NBody model(NBody::Runge_Kutta);
#endif

    model.configure(params);
    model.integrate();

    return 0;
}

