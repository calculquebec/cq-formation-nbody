This N body gravitational code is written in C++, the Makefile contains lines 
where you can specify the C++ compiler as well as any necessary compilation and 
linker flags. The resulting binary is called "nbody" and should be run with the 
command 
./nbody [parameter file]
An example of a parameter file is provided (parameters.txt); if no parameter 
file is supplied to the program, it runs with the default parameter values 
that are shown in the global.h header file. The program periodically writes 
out the particle positions to an MDL MOL file that can be visualized using a 
free tool like Jmol or PyMOL. 

The program is based on a direct (particle-particle) method and uses either 
the velocity Verlet or the fourth-order Runge-Kutta algorithm to numerically 
integrate the first-order system
x_i' = v_i
v_i' = \sum{j=1,j\ne i}^N m_j (x_i - x_j)/(eps + r_{ij}^2)^(3/2)
where x_i and v_i are the position and velocity of the i-th particle, m_i its
mass and r_ij is the L2 distance between x_i and x_j. The force softening
parameter eps is included to avoid singularities caused by near collisions of 
the particles and should be small and positive. The integration algorithm 
used is controlled by the compilation flag VERLET. The particle motions are 
confined to a toroidal box when the parameter "finite_domain" is set to true. 
The console output consists of the time value and total energy per particle, 
written out every hundred timesteps. 
      
