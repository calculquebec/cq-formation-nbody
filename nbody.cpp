#include <cassert>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>

#include "nbody.h"


NBody::NBody(const Algo algo_): rd(), gen(rd()), VRG(0.0, 1.0), algo(algo_)
{
    // All particules
    particules.resize(DEFAULT_NP);

    // Time attributes: Delta T, N steps, write result freq, epsilon limit
    dt = DEFAULT_DT;
    NT = DEFAULT_NT;
    write_freq = DEFAULT_WRITE_FREQ;
    epsilon = DEFAULT_EPSILON;

    // Limits for the masses
    low_mass = DEFAULT_LOW_MASS;
    high_mass = DEFAULT_HIGH_MASS;

    // Options
    center_masses = DEFAULT_CENTER_MASSES;
    bounded_state = DEFAULT_BOUNDED_STATE;
    finite_domain = DEFAULT_FINITE_DOMAIN;

    // Limits for the 3D space
    L.min = DEFAULT_L_MIN;
    L.max = DEFAULT_L_MAX;
}


void NBody::configure(const Params &params)
{
    std::string key;

    // Now that we have all input parameters, see if they match known parameters
    // If so, read in the value and assign it

    key = "seed";
    if (params.find(key) != params.end()) {
        int seed = stoi(params.at(key));
        assert(seed >= 0);
        if (seed == 0) seed = std::time(NULL);
        gen.seed(seed);
    }

    key = "nparticle";
    if (params.find(key) != params.end()) {
        particules.resize(stoi(params.at(key)));
        assert(particules.size() > 1);
    }

    key = "timestep";
    if (params.find(key) != params.end()) {
        dt = stod(params.at(key));
        assert(dt > std::numeric_limits<double>::epsilon());
    }

    key = "max_time";
    if (params.find(key) != params.end()) {
        double totalTime = stod(params.at(key));
        assert(totalTime > std::numeric_limits<double>::epsilon());
        NT = (int)(totalTime / dt);
    }

    key = "write_frequency";
    if (params.find(key) != params.end()) {
        write_freq = stoi(params.at(key));
        assert(write_freq > 0);
    }

    key = "epsilon";
    if (params.find(key) != params.end()) {
        epsilon = stod(params.at(key));
        assert(epsilon > std::numeric_limits<double>::epsilon() && epsilon < 0.1);
    }

    key = "min_mass";
    if (params.find(key) != params.end()) {
        low_mass = stod(params.at(key));
        assert(low_mass > std::numeric_limits<double>::epsilon());
    }

    key = "max_mass";
    if (params.find(key) != params.end()) {
        high_mass = stod(params.at(key));
        assert(high_mass >= low_mass);
    }

    key = "center_of_mass";
    if (params.find(key) != params.end()) {
        center_masses = (params.at(key) == "yes") ? true : false;
    }

    key = "bound_state";
    if (params.find(key) != params.end()) {
        bounded_state = (params.at(key) == "yes") ? true : false;
    }

    key = "finite_domain";
    if (params.find(key) != params.end()) {
        finite_domain = (params.at(key) == "yes") ? true : false;
        assert((finite_domain && center_masses) == false);
    }

    key = "xmin";
    if (params.find(key) != params.end()) {
        L.min.x = stod(params.at(key));
    }

    key = "xmax";
    if (params.find(key) != params.end()) {
        L.max.x = stod(params.at(key));
        assert(L.min.x < L.max.x);
    }

    key = "ymin";
    if (params.find(key) != params.end()) {
        L.min.y = stod(params.at(key));
    }

    key = "ymax";
    if (params.find(key) != params.end()) {
        L.max.y = stod(params.at(key));
        assert(L.min.y < L.max.y);
    }

    key = "zmin";
    if (params.find(key) != params.end()) {
        L.min.z = stod(params.at(key));
    }

    key = "zmax";
    if (params.find(key) != params.end()) {
        L.max.z = stod(params.at(key));
        assert(L.min.z < L.max.z);
    }
}


void NBody::printConfig() const
{
    std::cout << "nparticle = " << particules.size() << std::endl;

    std::cout << "timestep = " << dt << std::endl;
    std::cout << "max_time = " << NT * dt << std::endl;
    std::cout << "write_frequency = " << write_freq << std::endl;
    std::cout << "epsilon = " << epsilon << std::endl;

    std::cout << "min_mass = " << low_mass << std::endl;
    std::cout << "max_mass = " << high_mass << std::endl;

    std::cout << "centre_of_mass = " << (center_masses ? "yes" : "no") << std::endl;
    std::cout << "bound_state = " << (bounded_state ? "yes" : "no") << std::endl;
    std::cout << "finite_domain = " << (finite_domain ? "yes" : "no") << std::endl;

    std::cout << "xmin = " << L.min.x << std::endl;
    std::cout << "xmax = " << L.max.x << std::endl;
    std::cout << "ymin = " << L.min.y << std::endl;
    std::cout << "ymax = " << L.max.y << std::endl;
    std::cout << "zmin = " << L.min.z << std::endl;
    std::cout << "zmax = " << L.max.z << std::endl;
}


void NBody::integrate()
{
    const size_t NP = particules.size();
    const Vect3d V_min(-0.2, -0.2, -0.2);
    const Vect3d V_max( 0.2,  0.2,  0.2);

    // Assign initial values...
    for (size_t i = 0; i < NP; ++i) {
        Vect3d newPos, newVel;

        // For x, y and z
        for (size_t c = 0; c < 3; c++) {
            newPos.comp[c] = random1d(L.min.comp[c], L.max.comp[c]);
            newVel.comp[c] = random1d(V_min.comp[c], V_max.comp[c]);
        }

        // Initial position
        particules[i].setPos(newPos);

        // Initial speed - add a rotation around the z axis
        particules[i].setVel(newVel + Vect3d(-newPos.y, newPos.x, 0) / 10);
    }

    // For each particule
    for (size_t i = 0; i < NP; ++i) {
        // Assign random mass
        particules[i].setMass(random1d(low_mass, high_mass));
    }

    if (center_masses) {
        // Set the center of mass and it's speed to 0
        centerAll();
    }

    if (bounded_state) {
        // Make sure that the total energy of the system is negative
        // so particle don't fly in the distance
        // Set the kinetic energy to half the potential energy
        double U = totalPotentialEnergy();
        double K = totalKineticEnergy();
        double alpha = std::sqrt(0.5 * U / K);

        for (size_t i = 0; i < NP; ++i) {
             particules[i].v() *= alpha;
        }
    }

    if (algo == Verlet) {
        integrateVerlet();
    }
    else if (algo == Runge_Kutta) {
        integrateRungeKutta();
    }
    else {
        std::cerr << "Unsupported algorithm : " << algo << std::endl;
    }
}


void NBody::integrateVerlet()
{
    const size_t NP = particules.size();

    computeAccelerations();

    for (int iter = 0; iter < NT; ++iter) {
        writeState(iter);

        for (size_t i = 0; i < NP; ++i) {
            particules[i].p() +=       particules[i].v() * dt +
                                 0.5 * particules[i].a() * dt * dt;
        }

        boundaryConditions();

        for (size_t i = 0; i < NP; ++i) {
            particules[i].v() += 0.5 * particules[i].a() * dt;
        }

        computeAccelerations();

        for (size_t i = 0; i < NP; ++i) {
            particules[i].v() += 0.5 * particules[i].a() * dt;
        }
    }

    writeState(NT);
}


void NBody::integrateRungeKutta()
{
    const size_t NP = particules.size();
    std::vector<FourthOrderRungeKutta> r(NP);

    for (int iter = 0; iter < NT; ++iter) {
        writeState(iter);

        // Save all positions
        for (size_t i = 0; i < NP; ++i) {
            r[i].pos = particules[i].p();
        }
        computeAccelerations();

        // Compute k1's
        for (size_t i = 0; i < NP; ++i) {
            r[i].k1.vel = particules[i].v();
            r[i].k1.acc = particules[i].a();
            particules[i].p() = r[i].pos + 0.5 * dt * r[i].k1.vel;
        }
        computeAccelerations();

        // Compute k2's
        for (size_t i = 0; i < NP; ++i) {
            r[i].k2.vel = (1.0 + 0.5 * dt) * particules[i].v();
            r[i].k2.acc = particules[i].a();
            particules[i].p() = r[i].pos + 0.5 * dt * r[i].k2.vel;
        }
        computeAccelerations();

        // Compute k3's
        for (size_t i = 0; i < NP; ++i) {
            r[i].k3.vel = (1.0 + 0.5 * dt * (1.0 + 0.5 * dt)) * particules[i].v();
            r[i].k3.acc = particules[i].a();
            particules[i].p() = r[i].pos + dt * r[i].k3.vel;
        }
        computeAccelerations();

        // Compute k4's and final position and velocity
        for (size_t i = 0; i < NP; ++i) {
            r[i].k4.vel = (1.0 + dt * (1.0 + 0.5 * dt * (1.0 + 0.5 * dt))) * particules[i].v();
            r[i].k4.acc = particules[i].a();

            particules[i].p() = r[i].pos +
                dt * (r[i].k1.vel + 2.0 * r[i].k2.vel + 2.0 * r[i].k3.vel + r[i].k4.vel) / 6.0;
            particules[i].v() +=
                dt * (r[i].k1.acc + 2.0 * r[i].k2.acc + 2.0 * r[i].k3.acc + r[i].k4.acc) / 6.0;
        }

        boundaryConditions();
    }

    writeState(NT);
}


double NBody::random1d(double a, double b)
{
    return (a + (b - a) * VRG(gen));
}


Vect3d NBody::random3d(const Vect3d &a, const Vect3d &b)
{
    return Vect3d(random1d(a.x, b.x), random1d(a.y, b.y), random1d(a.z, b.z));
}


void NBody::boundaryConditions()
{
    if (!finite_domain) return;

    const size_t NP = particules.size();

    // Block_Size = Limit.max - Limit.min
    const Vect3d bSize = L.max - L.min;

    for (size_t i = 0; i < NP; ++i) {
        // Relative_Position = Position - Limit.min
        Vect3d relPos = particules[i].p() - L.min;

        // Nb_Blocks = Relative_Position / Block_Size
        // Nb_Blocks_I = floor(Nb_Blocks)
        // Delta = Nb_Blocks_I * Block_Size
        // New_Position = Position - Delta
        particules[i].p() -= bSize * floor(relPos / bSize);
    }
}


void NBody::centerAll()
{
    const size_t NP = particules.size();
    double total_mass = 0.0;
    Vect3d centralPos(0.0, 0.0, 0.0);
    Vect3d centralVel(0.0, 0.0, 0.0);

    for (size_t i = 0; i < NP; ++i) {
        total_mass += particules[i].m();
        centralPos += particules[i].p() * particules[i].m();
        centralVel += particules[i].v() * particules[i].m();
    }

    centralPos /= total_mass;
    centralVel /= total_mass;

    for (size_t i = 0; i < NP; ++i) {
        particules[i].p() -= centralPos;
        particules[i].v() -= centralVel;
    }
}


double NBody::totalPotentialEnergy() const
{
    const size_t NP = particules.size();
    double U = 0.0;

    for (size_t i = 0; i < NP; ++i) {
        for (size_t j = i + 1; j < NP; ++j) {
            U += particules[i].potentialEnergy(particules[j], epsilon);
        }
    }

    return U;
}


double NBody::totalKineticEnergy() const
{
    const size_t NP = particules.size();
    double K = 0.0;

    for (size_t i = 0; i < NP; ++i) {
        K += particules[i].kineticEnergy();
    }

    return K;
}


void NBody::computeAccelerations()
{
    const size_t NP = particules.size();

    for (size_t i = 0; i < NP; ++i) {
        Vect3d sum(0.0, 0.0, 0.0);

        for (size_t j = 0; j < NP; ++j) {
            if (i == j) continue;

            Vect3d delta = particules[j].p() - particules[i].p();
            double rij = std::sqrt(epsilon + delta.dotProd(delta));
            sum += delta * particules[j].m() / (rij * rij * rij);
        }

        particules[i].setAcc(sum);
    }
}


void NBody::writeState(const int iter)
{
    const size_t NP = particules.size();

    if (iter % 100 == 0) {
        std::cout << double(iter) * dt << "  "
            << (totalKineticEnergy() - totalPotentialEnergy()) / NP
            << std::endl;
    }

    if (iter % write_freq == 0) {
        int zero = 0;
        std::stringstream sstream;
        sstream << "nbody_" << iter << ".mol";

        std::string filename = sstream.str();
        std::ofstream s(filename.c_str());

        s << "nbody_" << iter << std::endl;
        s << "  MOE2000" << std::endl;
        s << std::endl;
        s << std::setw(3) << NP;
        s << std::setw(3) << zero;
        s << " 0  0  0  0  0  0  0  0   1 V2000" << std::endl;

        for (size_t i = 0; i < NP; ++i) {
            s << particules[i].p() << " ";
            s << std::setw(3) << std::setiosflags(std::ios::left) << "C";
            s << std::resetiosflags(std::ios::left) << " 0  0  0  0  0  0  0  0  0  0  0  0" << std::endl;
        }

        s << "M  END" << std::endl;
        s << "$$$$" << std::endl;
        s.close();
    }
}


