#include "global.h"
#include <omp.h>


double drandom(double x,double y)
{
  double out = VRG(gen);
  out = x + (y - x)*out;
  return out;
}

void write_state(int timestep,const double* x)
{
  int i,zero = 0;
  std::stringstream sstream;
  sstream << "nbody_" << timestep << ".mol";
  std::string filename = sstream.str();
  std::ofstream s(filename.c_str());
  s << "nbody_" << timestep << std::endl;
  s << "  MOE2000" << std::endl;
  s << std::endl;
  s << std::setw(3) << NP << std::setw(3) << zero << " 0  0  0  0  0  0  0  0   1 V2000" << std::endl;
  for(i=0; i<NP; ++i) {
    s << std::setw(10) << std::setprecision(4) << std::setiosflags(std::ios::fixed) << x[3*i] << std::setw(10) << std::setprecision(4) << std::setiosflags(std::ios::fixed) << x[3*i+1] << std::setw(10) << std::setprecision(4) << std::setiosflags(std::ios::fixed) << x[3*i+2] << " " << std::setw(3) << std::setiosflags(std::ios::left) << "C" << std::resetiosflags(std::ios::left) << " 0  0  0  0  0  0  0  0  0  0  0  0" << std::endl;
  }
  s << "M  END" << std::endl;
  s << "$$$$" << std::endl;
  s.close();
}

void boundary_conditions(double* x)
{
  int i,j;
  double xv,dsize;

  if (!finite_domain) return;

  for(i=0; i<NP; ++i) {
    for(j=0; j<3; ++j) {
      xv = x[3*i+j];
      if (xv < L[2*j]) {
        dsize = L[2*j+1] - L[2*j];
        do {
          xv += dsize;
          if (xv > L[2*j]) break;
        } while(true);
        x[3*i+j] = xv;
      }
      else if (xv > L[2*j+1]) {
        dsize = L[2*j+1] - L[2*j];
        do {
          xv -= dsize;
          if (xv < L[2*j+1]) break;
        } while(true);
        x[3*i+j] = xv;
      }
    }
  }
}

void compute_acceleration(const double* x,const double* mass,double* acc,const int* range)
{
  int i,j,k;
  double delta,rij,pfactor,sum[3];
  for(i=range[0]; i<range[1]; ++i) {
    for(j=0; j<3; ++j) {
      sum[j] = 0.0;
    }
    for(j=0; j<NP; ++j) {
      if (i == j) continue;
      delta = 0.0;
      for(k=0; k<3; ++k) {
        delta += (x[3*i+k] - x[3*j+k])*(x[3*i+k] - x[3*j+k]);
      }
      rij = std::sqrt(epsilon + delta);
      pfactor = mass[j]/(rij*rij*rij);
      for(k=0; k<3; ++k) {
        sum[k] += pfactor*(x[3*i+k] - x[3*j+k]);
      }
    }
    for(j=0; j<3; ++j) {
      acc[3*i+j] = -sum[j];
    }
  }
}


double compute_energy(const double* x,const double* v,const double* mass,const int* range)
{
  int i,j,k;
  double delta,T = 0.0,U = 0.0;

  // First the kinetic energy...
  for(i=range[0]; i<range[1]; ++i) {
    delta = 0.0;
    for(j=0; j<3; ++j) {
      delta += v[3*i+j]*v[3*i+j];
    }
    T += mass[i]*delta;
  }

  // Now the potential energy
  for(i=range[0]; i<range[1]; ++i) {
    for(j=1+i; j<NP; ++j) {
      delta = 0.0;
      for(k=0; k<3; ++k) {
        delta += (x[3*i+k] - x[3*j+k])*(x[3*i+k] - x[3*j+k]);
      }
      U += mass[i]*mass[j]/std::sqrt(epsilon + delta); 
    }
  }
  return (0.5*T - U);  
}

double compute_kinetic_energy(const double* v,const double* mass)
{
  int i,j;
  double delta,T = 0.0;

  for(i=0; i<NP; ++i) {
    delta = 0.0;
    for(j=0; j<3; ++j) {
      delta += v[3*i+j]*v[3*i+j];
    }
    T += mass[i]*delta;
  }
  return 0.5*T;  
}

double compute_potential_energy(const double* x, const double* mass)
{
  int i,j,k;
  double delta,U = 0.0;

  for(i=0; i<NP; ++i) {
    for(j=1+i; j<NP; ++j) {
      delta = 0.0;
      for(k=0; k<3; ++k) {
        delta += (x[3*i+k] - x[3*j+k])*(x[3*i+k] - x[3*j+k]);
      }
      U += mass[i]*mass[j]/std::sqrt(epsilon + delta); 
    }
  }
  return U;  
}

void compute_center_of_mass(const double* x, const double* mass, double* center)
{
  int i,j;
  double total_mass = 0.0;

  for(j=0; j<3; ++j) {
    center[j] = 0.0;
  }
  for(i=0; i<NP; ++i) {
    for(j=0; j<3; ++j) {
      center[j] += x[3*i+j]*mass[i];
    }
    total_mass += mass[i];
  }
  for(j=0; j<3; ++j) {
    center[j] /= total_mass;
  }
}

void center_particles(double* x, const double* mass)
{
  int i,j;
  double center[3];

  compute_center_of_mass(x,mass,center);
  for(i=0; i<NP; ++i) {
    for(j=0; j<3; ++j) {
      x[3*i+j] -= center[j];
    }
  }
}

void integrate()
{
  int i,j,l;
  double x[3*NP],xnew[3*NP],v[3*NP],vnew[3*NP];
  double mass[NP],acc[3*NP],temp[3*NP];
  double K,U,alpha,E;
  int *range;
  int n_threads, tread_n;

  // Assign initial values...
  for(i=0; i<NP; ++i) {
    for(j=0; j<3; ++j) {
      // Initial position
      x[3*i+j] = drandom(L[2*j],L[2*j+1]);
      v[3*i+j] = drandom(-0.2,0.2);
    }
  }
  for(i=0; i<NP; ++i) {
    mass[i] = drandom(low_mass,high_mass);
  }

  if (center_masses) {
    // Set the center of mass and it's speed to 0
    center_particles(x,mass);
    center_particles(v,mass);
  }
  
  for(i=0; i<NP; ++i) {
    v[3*i+1] += x[3*i+0]/10.0;
    v[3*i+0] -= x[3*i+1]/10.0;
  }
  
  if (bounded_state) {
    // Make sure that the total energy of the system is negative so particle don't fly in the distance
    // set the kinetic energy to half the potential energy
    U = compute_potential_energy(x,mass);
    K = compute_kinetic_energy(v,mass);
    alpha = std::sqrt(U/(2.0*K));

    for(i=0; i<NP; ++i) {
      for(j=0; j<3; ++j) {
        v[3*i+j] *= alpha;
      }
    }
  }
  
  write_state(0,x);
  std::cout << "0.0  " << (compute_energy_k(v,mass)-compute_energy_p(x,mass))/double(NP) << std::endl;

#ifdef VERLET
    #pragma omp parallel private(i,j,tread_n)
    {
      tread_n = omp_get_thread_num();
      n_threads = omp_get_num_threads();
      #pragma omp single
      {
        range = new int[n_threads+1];
        for(i=0;i<n_threads+1;++i){
          range[i] = NP*i/n_threads;
        }
      }

      compute_acceleration(x,mass,acc, range+tread_n);
      #pragma omp barrier
      l=1;
      while(l<=NT){

        // Print out the system's total energy per particle (should be fairly constant)
        #pragma omp master
        if (l%100 == 0) {
            E = 0;
            //#pragma omp for reduction(+:E)
            for(i=0;i<n_threads;++i){
                E += compute_energy(x,v,mass, &range[i]);
            } 
            std::cout << dt*double(l) << "  " << E/double(NP) << std::endl;
        }

        #pragma omp for collapse(2) schedule(static)
        for(i=0; i<NP; ++i) {
          for(j=0; j<3; ++j) {
            xnew[3*i+j] = x[3*i+j] + dt*v[3*i+j] + 0.5*dt*dt*acc[3*i+j];
          }
        }
        if (finite_domain) boundary_conditions(xnew);

        compute_acceleration(xnew,mass,temp, range+tread_n);

        #pragma omp barrier
        #pragma omp for collapse(2) schedule(static) 
        for(i=0; i<NP; ++i) {
          for(j=0; j<3; ++j) {
            vnew[3*i+j] = v[3*i+j] + 0.5*dt*(acc[3*i+j] + temp[3*i+j]);
          }
        }
        
        #pragma omp master 
        if ((l % write_freq) == 0) write_state(l,xnew);
        
        // Now update the arrays
        #pragma omp for collapse(2) schedule(static)
        for(i=0; i<NP; ++i) {
          for(j=0; j<3; ++j) {
            x[3*i+j] = xnew[3*i+j];
            v[3*i+j] = vnew[3*i+j];
            acc[3*i+j] = temp[3*i+j];
          }
        }
               
        #pragma omp master 
        l++;
        #pragma omp barrier
      }
  }
#else
  // Fourth-order Runge-Kutta
  double k1[6*NP],k2[6*NP],k3[6*NP],k4[6*NP];
  
  for(l=1; l<=NT; ++l) {
    compute_acceleration(x,mass,acc);
    for(i=0; i<NP; ++i) {
      for(j=0; j<3; ++j) {
        k1[3*i+j] = v[3*i+j];
        k1[3*NP+3*i+j] = acc[3*i+j];
        temp[3*i+j] = x[3*i+j] + 0.5*dt*k1[3*i+j];
      }
    }
    compute_acceleration(temp,mass,acc);
    for(i=0; i<NP; ++i) {
      for(j=0; j<3; ++j) {
        k2[3*i+j] = (1.0 + 0.5*dt)*v[3*i+j];
        k2[3*NP+3*i+j] = acc[3*i+j];
        temp[3*i+j] = x[3*i+j] + 0.5*dt*k2[3*i+j];
      }
    }
    compute_acceleration(temp,mass,acc);
    for(i=0; i<NP; ++i) {
      for(j=0; j<3; ++j) {
        k3[3*i+j] = (1.0 + 0.5*dt + 0.25*dt*dt)*v[3*i+j];
        k3[3*NP+3*i+j] = acc[3*i+j];
        temp[3*i+j] = x[3*i+j] + dt*k3[3*i+j];
      }
    }
    compute_acceleration(temp,mass,acc);
    for(i=0; i<NP; ++i) {
      for(j=0; j<3; ++j) {
        k4[3*i+j] = (1.0 + dt + 0.5*dt*dt + 0.25*dt*dt*dt)*v[3*i+j];
        k4[3*NP+3*i+j] = acc[3*i+j];
      }
    }
    for(i=0; i<NP; ++i) {
      for(j=0; j<3; ++j) {
        xnew[3*i+j] = x[3*i+j] + dt*(k1[3*i+j] + 2.0*k2[3*i+j] + 2.0*k3[3*i+j]+k4[3*i+j])/6.0;
        vnew[3*i+j] = v[3*i+j] + dt*(k1[3*NP+3*i+j] + 2.0*k2[3*NP+3*i+j] + 2.0*k3[3*NP+3*i+j] + k4[3*NP+3*i+j])/6.0;
      }
    }
    boundary_conditions(xnew);
    // Print out the system's total energy per particle (should be fairly constant)
    if (l%100 == 0) {
      std::cout << dt*double(l) << "  " << compute_energy(xnew,vnew,mass)/double(NP) << std::endl;
    }
    if ((l % write_freq) == 0) write_state(l,xnew);
    // Now update the arrays
    for(i=0; i<NP; ++i) {
      for(j=0; j<3; ++j) {
        x[3*i+j] = xnew[3*i+j];
        v[3*i+j] = vnew[3*i+j];
      }
    }
  }
#endif
  write_state(NT,x);
}

void read_parameters(const char* filename)
{
  unsigned int i,bpoint;
  double tvalue;
  std::string line,name,value;

  std::ifstream s(filename);
  if (!s.is_open()) {
    // If the file doesn't exist, we need to exit...
    std::cout << "The file " << filename << " cannot be found!" << std::endl;
    std::exit(1);
  }

  // Loop through all lines in the parameter file
  while(std::getline(s,line)) {
    // If it's an empty line, continue
    if (line.empty()) continue;
    // If the line begins with a #, ignore it
    if (line[0] == '#') continue;
    // If there's no equals sign in this line, continue
    if (line.find('=') == std::string::npos) continue;
    // Assumes that the equals sign can only occur once in 
    // the line
    bpoint = 0;
    name = "";
    for(i=0; i<line.size(); ++i) {
      if (line[i] == ' ') continue;
      if (line[i] == '=') {
        bpoint = i;
        break;
      }
      name += line[i];
    }
    value = "";
    for(i=1+bpoint; i<line.size(); ++i) {
      if (line[i] == ' ') continue;
      value += line[i];
    }
    // Now that we have the parameter name, see if it matches
    // any of the known parameters. If so, read in the value and
    // assign it
    if (name == "nparticle") {
      NP = stoi(value);
    }
    else if (name == "max_time") {
      tvalue = stod(value);
    }
    else if (name == "seed") {
      seed = stoi(value);
    }
    else if (name == "timestep") {
      dt = stod(value);
    }
    else if (name == "epsilon") {
      epsilon = stod(value);
    }
    else if (name == "min_mass") {
      low_mass = stod(value);
    }
    else if (name == "max_mass") {
      high_mass = stod(value);
    }
    else if (name == "write_frequency") {
      write_freq = stoi(value);
    }
    else if (name == "finite_domain") {
      finite_domain = (value == "yes") ? true : false;
    }
    else if (name == "center_of_mass") {
      center_masses = (value == "yes") ? true : false;
    }
    else if (name == "bound_state") {
      bounded_state = (value == "yes") ? true : false;
    }
    else if (name == "xmin") {
      L[0] = stod(value);
    }
    else if (name == "xmax") {
      L[1] = stod(value);
    }
    else if (name == "ymin") {
      L[2] = stod(value);
    }
    else if (name == "ymax") {
      L[3] = stod(value);
    }
    else if (name == "zmin") {
      L[4] = stod(value);
    }
    else if (name == "zmax") {
      L[5] = stod(value);
    }
  }
  s.close();
  // Sanity checks
  assert(tvalue > std::numeric_limits<double>::epsilon());
  assert(NP > 1);
  assert(dt > std::numeric_limits<double>::epsilon());
  assert(epsilon > std::numeric_limits<double>::epsilon() && epsilon < 0.1);
  assert(write_freq > 0);
  assert(low_mass > std::numeric_limits<double>::epsilon());
  assert(high_mass >= low_mass);
  assert(seed >= 0);
  for(int i=0; i<3; ++i) {
    assert(L[2*i+1] > L[2*i]);
  }
  if (seed == 0) seed = std::time(NULL);
  gen.seed(seed);
  NT = int(tvalue/dt);
}

int main(int argc,char** argv)
{
  if (argc > 2) {
    std::cerr << "Usage: ./nbody parameters.txt" << std::endl;
    return 0;
  }
  if (argc == 2) read_parameters(argv[1]);

  integrate();

  return 0;
}
