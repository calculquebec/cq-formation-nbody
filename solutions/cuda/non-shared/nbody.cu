#include "global.h"

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

__global__ void boundary_conditions_gpu(int NP, bool finite_domain, const double *L, double *x)
{
   int i,j;
   double xv,dsize;
   i=blockIdx.x*blockDim.x + threadIdx.x;

   if(finite_domain){
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

   __syncthreads();
  }
}

void compute_acceleration(const double* x,const double* mass,double* acc)
{
  int i,j,k;
  double delta,rij,pfactor,sum[3];

  for(i=0; i<NP; ++i) {
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

__global__ void compute_acceleration_gpu(int NP,double epsilon,const double *x,const double *mass,double *acc)
{
   int i,j,k;
   double delta,rij,pfactor,sum[3];
   i=blockIdx.x*blockDim.x + threadIdx.x;
   for(k=0; k<3; ++k) sum[k] = 0.0;
   for(j=0;j<NP;++j){
	if(i!=j){
		delta = 0.0;
		for(k=0; k<3; ++k) delta += (x[3*i+k] - x[3*j+k])*(x[3*i+k] - x[3*j+k]);
	//	rij = std::sqrt(epsilon + delta);
		rij = __dsqrt_rd (epsilon + delta);
		pfactor = mass[j]/(rij*rij*rij);
		for(k=0; k<3; ++k) sum[k] += pfactor*(x[3*i+k] - x[3*j+k]);	
	}
   }
   __syncthreads();
   for(k=0; k<3; ++k) acc[3*i+k] = -sum[k];
   __syncthreads();
}

double compute_energy(const double* x,const double* v,const double* mass)
{
  int i,j,k;
  double delta,T = 0.0,U = 0.0;

  // First the kinetic energy...
  for(i=0; i<NP; ++i) {
    delta = 0.0;
    for(j=0; j<3; ++j) {
      delta += v[3*i+j]*v[3*i+j];
    }
    T += mass[i]*delta;
  }
 
  // Now the potential energy
  for(i=0; i<NP; ++i) {
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


__global__ void compute_energy_gpu(int NP,const double epsilon,const double *x,const double *v,const double *mass,double *E)
{
  int i,j,k;
  double delta;
  i=blockIdx.x*blockDim.x + threadIdx.x;

  // First the kinetic energy...
  delta = 0.0;
  for(j=0; j<3; ++j) delta += v[3*i+j]*v[3*i+j];
  E[i] = 0.5*mass[i]*delta;

  // Now the potential energy
  for(j=1+i; j<NP; ++j) {
  	delta = 0.0;
	for(k=0; k<3; ++k) delta += (x[3*i+k] - x[3*j+k])*(x[3*i+k] - x[3*j+k]);
	E[i] -= mass[i]*mass[j]/std::sqrt(epsilon + delta);
  }
  __syncthreads(); 
}

double compute_kinetic_energy(const double* x,const double* v,const double* mass)
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

double compute_potential_energy(const double* x,const double* v,const double* mass)
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

__global__ void get_new_coordinates_gpu(const double dt,const double *x,const double *v, const double *acc, double *xnew)
{
  int i=blockIdx.x*blockDim.x + threadIdx.x;
  for(int j=0; j<3; ++j) {
  	xnew[3*i+j] = x[3*i+j] + dt*v[3*i+j] + 0.5*dt*dt*acc[3*i+j];
//	xnew[3*i+j] = x[3*i+j] dt*v[3*i+j];
  }
  __syncthreads();
  
}

__global__ void get_new_velocities_gpu(const double dt,const double *v, const double *temp,const double *acc, double *vnew)
{
  int i=blockIdx.x*blockDim.x + threadIdx.x; 
      for(int j=0; j<3; ++j) {
        vnew[3*i+j] = v[3*i+j] + 0.5*dt*(acc[3*i+j] + temp[3*i+j]);
      }
  __syncthreads();
}

__global__ void update_variables_gpu(const double *xnew, const double *vnew, const double *temp, double *x, double *v, double *acc)
{
  int i=blockIdx.x*blockDim.x + threadIdx.x;
      for(int j=0; j<3; ++j) {
        x[3*i+j] = xnew[3*i+j];
        v[3*i+j] = vnew[3*i+j];
        acc[3*i+j] = temp[3*i+j];
      }
  __syncthreads();  
}

__global__ void kernel_RK1_gpu(const int NP,const double dt,const double *x,const double *v,const double *acc,double *k1,double *temp)
{
  int i=blockIdx.x*blockDim.x + threadIdx.x;
      for(int j=0; j<3; ++j) {
        k1[3*i+j] = v[3*i+j];
        k1[3*NP+3*i+j] = acc[3*i+j];
        temp[3*i+j] = x[3*i+j] + 0.5*dt*k1[3*i+j];
      }
  __syncthreads();
}

__global__ void kernel_RK2_gpu(const int NP,const double dt,const double *x,const double *v,const double *acc,double *k2,double *temp)
{ 
  int i=blockIdx.x*blockDim.x + threadIdx.x;
     for(int j=0; j<3; ++j) {
        k2[3*i+j] = (1.0 + 0.5*dt)*v[3*i+j];
        k2[3*NP+3*i+j] = acc[3*i+j];
        temp[3*i+j] = x[3*i+j] + 0.5*dt*k2[3*i+j];
      }
  __syncthreads();
}

__global__ void kernel_RK3_gpu(const int NP,const double dt,const double *x,const double *v,const double *acc,double *k3,double *temp)
{
  int i=blockIdx.x*blockDim.x + threadIdx.x;
  for(int j=0; j<3; ++j) {
        k3[3*i+j] = (1.0 + 0.5*dt + 0.25*dt*dt)*v[3*i+j];
        k3[3*NP+3*i+j] = acc[3*i+j];
        temp[3*i+j] = x[3*i+j] + dt*k3[3*i+j];
      }
  __syncthreads();
}

__global__ void kernel_RK4_gpu(const int NP,const double dt,const double *x,const double *v,const double *acc,double *k4)
{
  int i=blockIdx.x*blockDim.x + threadIdx.x;
  for(int j=0; j<3; ++j) {
        k4[3*i+j] = (1.0 + dt + 0.5*dt*dt + 0.25*dt*dt*dt)*v[3*i+j];
        k4[3*NP+3*i+j] = acc[3*i+j];
      }
  __syncthreads();
}

__global__ void get_newXV_RK_gpu(const int NP,const double dt,const double *x,const double *v,const double *k1,const double *k2,const double *k3,const double *k4, double *xnew, double *vnew)
{
  int i=blockIdx.x*blockDim.x + threadIdx.x;
      for(int j=0; j<3; ++j) {
        xnew[3*i+j] = x[3*i+j] + dt*(k1[3*i+j] + 2.0*k2[3*i+j] + 2.0*k3[3*i+j]+k4[3*i+j])/6.0;
        vnew[3*i+j] = v[3*i+j] + dt*(k1[3*NP+3*i+j] + 2.0*k2[3*NP+3*i+j] + 2.0*k3[3*NP+3*i+j] + k4[3*NP+3*i+j])/6.0;
      }
  __syncthreads();
}

__global__ void update_variables_RK_gpu(const double *xnew, const double *vnew, double *x,double *v)
{
  int i=blockIdx.x*blockDim.x + threadIdx.x;
  for(int j=0; j<3; ++j) {
        x[3*i+j] = xnew[3*i+j];
        v[3*i+j] = vnew[3*i+j];
      }
  __syncthreads();
}


void integrate()
{
  int i,j,l;
  double x[3*NP],xnew[3*NP],v[3*NP];
//  double vnew[3*NP],acc[3*NP],temp[3*NP];
  double mass[NP];
  double K,U,alpha;

  // Define GPU variables
  double *d_x, *d_xnew, *d_v, *d_vnew;
  double *d_mass, *d_acc, *d_temp;
  double *E, *d_E;
  double *d_L;
  double Ef;

  // Allocate GPU variables
  size_t memSize=sizeof(double)*NP*3;
  cudaMalloc((void**)&d_x,sizeof(double)*NP*3);
  cudaMalloc((void**)&d_xnew,sizeof(double)*NP*3);
  cudaMalloc((void**)&d_v,sizeof(double)*NP*3);
  cudaMalloc((void**)&d_vnew,sizeof(double)*NP*3);
  cudaMalloc((void**)&d_acc,sizeof(double)*NP*3);
  cudaMalloc((void**)&d_temp,sizeof(double)*NP*3);
  cudaMalloc((void**)&d_mass,sizeof(double)*NP);
  cudaMalloc((void**)&d_E,sizeof(double)*NP);
  cudaMalloc((void**)&d_L,6*sizeof(double));

  E = (double*)malloc(sizeof(double)*NP);

  // Assign initial values...
  for(i=0; i<NP; ++i) {
    for(j=0; j<3; ++j) {
      // Initial position and speed
      x[3*i+j] = drandom(L[2*j],L[2*j+1]);
      v[3*i+j] = drandom(-0.2,0.2);
    }
  }
  // Assign random mass
  for(i=0; i<NP; ++i) {
    mass[i] = drandom(low_mass,high_mass);
  }

  // Add a rotation around the z axis
  for(i=0; i<NP; ++i) {
    v[3*i+1] += x[3*i+0]/10.0;
    v[3*i+0] -= x[3*i+1]/10.0;
  }

  if (center_masses) {
    // Set the center of mass and it's speed to 0
    center_particles(x, mass);
    center_particles(v, mass);
  }

  if (bounded_state) {
    // Make sure that the total energy of the system is negative so particle don't fly in the distance
    // Set the kinetic energy to half the potential energy
    U = compute_potential_energy(x,v,mass);
    K = compute_kinetic_energy(x,v,mass);
    alpha = std::sqrt(U/(2.0*K));

    for(i=0; i<NP; ++i) {
      for(j=0; j<3; ++j) {
        v[3*i+j] *= alpha;
      }
    }
  }

  write_state(0,x);
  std::cout << "0.0  " << compute_energy(x,v,mass)/double(NP) << std::endl;

#ifdef VERLET
  clock_t start, end;
     double cpu_time_used;
     start = clock();
  int numThreadsPerBlock=256;
  int numBlocks=NP/numThreadsPerBlock;
  dim3 dimGrid(numBlocks);
  dim3 dimBlock(numThreadsPerBlock);

  // Copy data to GPU
  cudaMemcpy(d_x,&x[0],memSize,cudaMemcpyHostToDevice);
  cudaMemcpy(d_mass,&mass[0],sizeof(double)*NP,cudaMemcpyHostToDevice);
  cudaMemcpy(d_v,&v[0],memSize,cudaMemcpyHostToDevice);
  cudaMemcpy(d_L,&L[0],6*sizeof(double),cudaMemcpyHostToDevice);

  // Compute acceleration on GPU 
  compute_acceleration_gpu <<<dimGrid, dimBlock>>>(NP,epsilon,d_x,d_mass,d_acc);


  for(l=1; l<=NT; ++l) {

  // Get new coordinates done on GPU
  get_new_coordinates_gpu <<<dimGrid, dimBlock>>>(dt,d_x,d_v,d_acc,d_xnew);


    // Apply boundary conditions 
    boundary_conditions_gpu <<<dimGrid, dimBlock>>> (NP,finite_domain,d_L,d_xnew); 

    // Compute acceleration on GPU 
    compute_acceleration_gpu <<<dimGrid, dimBlock>>>(NP,epsilon,d_xnew,d_mass,d_temp);


    // Get new velocities done on GPU
    get_new_velocities_gpu <<<dimGrid, dimBlock>>> (dt,d_v,d_temp,d_acc,d_vnew);

    // Print out the system's total energy per particle (should be fairly constant)
    if (l%100 == 0) {
      compute_energy_gpu <<<dimGrid, dimBlock>>> (NP,epsilon,d_x,d_v,d_mass,d_E);
      cudaMemcpy(E,d_E,sizeof(double)*NP,cudaMemcpyDeviceToHost);
      Ef=0.0;
      for(int i=0;i<NP;i++) Ef+=E[i];
      std::cout << dt*double(l) << "  " << Ef/double(NP) << std::endl;

    }
    if ((l % write_freq) == 0){
	cudaMemcpy(&xnew[0],d_xnew,memSize,cudaMemcpyDeviceToHost);
	write_state(l,xnew);
    }

    // Now update the arrays on GPU
    update_variables_gpu <<<dimGrid, dimBlock>>> (d_xnew,d_vnew,d_temp,d_x,d_v,d_acc);

  }
  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("Exec time = %f\n",cpu_time_used);
#else
  // Fourth-order Runge-Kutta
//  double k1[6*NP],k2[6*NP],k3[6*NP],k4[6*NP];

  clock_t start, end;
  double cpu_time_used;
  start = clock();
  int numThreadsPerBlock=256;
  int numBlocks=NP/numThreadsPerBlock;
  dim3 dimGrid(numBlocks);
  dim3 dimBlock(numThreadsPerBlock);

  double *d_k1, *d_k2, *d_k3, *d_k4;
  cudaMalloc((void**)&d_k1,sizeof(double)*6*NP);
  cudaMalloc((void**)&d_k2,sizeof(double)*6*NP);
  cudaMalloc((void**)&d_k3,sizeof(double)*6*NP);
  cudaMalloc((void**)&d_k4,sizeof(double)*6*NP);

  // Copy data to GPU
  cudaMemcpy(d_x,&x[0],memSize,cudaMemcpyHostToDevice);
  cudaMemcpy(d_mass,&mass[0],sizeof(double)*NP,cudaMemcpyHostToDevice);
  cudaMemcpy(d_v,&v[0],memSize,cudaMemcpyHostToDevice);
  cudaMemcpy(d_L,&L[0],6*sizeof(double),cudaMemcpyHostToDevice);

  
  for(l=1; l<=NT; ++l) {
    // Compute acceleration on GPU 
    compute_acceleration_gpu <<<dimGrid, dimBlock>>>(NP,epsilon,d_x,d_mass,d_acc);
//    cudaMemcpy(&acc[0],d_acc,memSize,cudaMemcpyDeviceToHost);

    // Compute coefficients K1 on GPU
    kernel_RK1_gpu <<< dimGrid, dimBlock >>> (NP,dt,d_x,d_v,d_acc,d_k1,d_temp);
//    cudaMemcpy(&temp[0],d_temp,memSize,cudaMemcpyDeviceToHost);
//    cudaMemcpy(&k1[0],d_k1,sizeof(double)*6*NP,cudaMemcpyDeviceToHost);

    // Compute acceleration on GPU 
    compute_acceleration_gpu <<<dimGrid, dimBlock>>>(NP,epsilon,d_temp,d_mass,d_acc);
//    cudaMemcpy(&acc[0],d_acc,memSize,cudaMemcpyDeviceToHost);

    // Compute coefficients K2 on GPU
    kernel_RK2_gpu <<< dimGrid, dimBlock >>> (NP,dt,d_x,d_v,d_acc,d_k2,d_temp);
//    cudaMemcpy(&temp[0],d_temp,memSize,cudaMemcpyDeviceToHost);
//    cudaMemcpy(&k2[0],d_k2,sizeof(double)*6*NP,cudaMemcpyDeviceToHost);

    // Compute acceleration on GPU 
    compute_acceleration_gpu <<<dimGrid, dimBlock>>>(NP,epsilon,d_temp,d_mass,d_acc);
//    cudaMemcpy(&acc[0],d_acc,memSize,cudaMemcpyDeviceToHost);

    // Compute coefficients K3 on GPU
    kernel_RK3_gpu <<< dimGrid, dimBlock >>> (NP,dt,d_x,d_v,d_acc,d_k3,d_temp);
//    cudaMemcpy(&temp[0],d_temp,memSize,cudaMemcpyDeviceToHost);
//    cudaMemcpy(&k3[0],d_k3,sizeof(double)*6*NP,cudaMemcpyDeviceToHost);

    // Compute acceleration on GPU 
    compute_acceleration_gpu <<<dimGrid, dimBlock>>>(NP,epsilon,d_temp,d_mass,d_acc);
//    cudaMemcpy(&acc[0],d_acc,memSize,cudaMemcpyDeviceToHost);
    
    // Compute coefficients K4 on GPU
    kernel_RK4_gpu <<< dimGrid, dimBlock >>> (NP,dt,d_x,d_v,d_acc,d_k4);
//    cudaMemcpy(&k4[0],d_k4,sizeof(double)*6*NP,cudaMemcpyDeviceToHost);


    // Get new xnew and vnew on GPU
    get_newXV_RK_gpu <<< dimGrid, dimBlock >>> (NP,dt,d_x,d_v,d_k1,d_k2,d_k3,d_k4, d_xnew, d_vnew);
//    cudaMemcpy(&xnew[0],d_xnew,memSize,cudaMemcpyDeviceToHost);
//    cudaMemcpy(&vnew[0],d_vnew,memSize,cudaMemcpyDeviceToHost);


    // Apply boundary conditions 
    boundary_conditions_gpu <<<dimGrid, dimBlock>>> (NP,finite_domain,d_L,d_xnew);
//    cudaMemcpy(&xnew[0],d_xnew,memSize,cudaMemcpyDeviceToHost);


    // Print out the system's total energy per particle (should be fairly constant)
    if (l%100 == 0) {
      compute_energy_gpu <<<dimGrid, dimBlock>>> (NP,epsilon,d_x,d_v,d_mass,d_E);
      cudaMemcpy(E,d_E,sizeof(double)*NP,cudaMemcpyDeviceToHost);
      Ef=0.0;
      for(int i=0;i<NP;i++) Ef+=E[i];
      std::cout << dt*double(l) << "  " << Ef/double(NP) << std::endl;
    //  std::cout << dt*double(l) << "  " << compute_energy(xnew,vnew,mass)/double(NP) << std::endl;
    }
    if ((l % write_freq) == 0){
	cudaMemcpy(&xnew[0],d_xnew,memSize,cudaMemcpyDeviceToHost);
	write_state(l,xnew);
    }

    // Now update the arrays on GPU
    update_variables_RK_gpu <<<dimGrid, dimBlock>>> (d_xnew,d_vnew,d_x,d_v);

  }
  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("Exec time = %f\n",cpu_time_used);
#endif
  write_state(NT,x);
  cudaFree(d_x);
  cudaFree(d_v);
  cudaFree(d_acc);
  cudaFree(d_xnew);
  cudaFree(d_vnew);
  cudaFree(d_temp);
  cudaFree(d_mass);
  cudaFree(d_E);
  cudaFree(d_L);
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

