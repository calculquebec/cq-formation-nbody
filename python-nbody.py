import sys
import numpy as np

# Number of particles
NP = 100
# Number of time steps
NT = 10
tvalue = 0.5
# Random number seed 
seed = 0
# Frequency of writing position values to disk
write_freq = 2
# Particle mass is drawn uniformly from the  
# interval [low_mass,high_mass)
low_mass = 1.0
high_mass = 5.0
# Time increment in ODE solver
dt = 0.05
# Force softening
epsilon = 0.00000000001
# Use a finite domain with toroidal boundary 
# conditions?
finite_domain = False
# The dimensions of the finite domain
L = np.array([0.0,100.0,0.0,100.0,0.0,50.0])
# Condition for "bound state" of the particles
bounded_state = True
# Fix the center of mass at 0
center_masses = True
# Use Runge-Kutta
VERLET = False

def drandom(x, y, size):
    out = np.random.randint(256**4, dtype='<u4', size=size)
    out = out / (256.0*256*256*256)
    out = x + (y - x)*out
    return out

def write_state(timestep, x):
    f = open("nbody_%s.mol" % timestep, "w")
    f.write("nbody_%s\n" % timestep)
    f.write("  MOE2000\n")
    f.write("\n")
    f.write(f"{x.shape[1]:3d}{0:3d} 0  0  0  0  0  0  0  0   1 V2000\n")
    for row in x:
        f.write(f"{row[0]:10.4f}{row[1]:10.4f}{row[2]:10.4f} C   0  0  0  0  0  0  0  0  0  0  0  0\n")
    f.write("M  END\n")
    f.write("$$$$\n")
    f.close()

def boundary_conditions(x):
    if finite_domain:
        minimum = L[::2]
        maximum = L[1::2]
        dsize = maximum - minimum
        x[:,:] = (x - minimum) % dsize + minimum
                    
def compute_acceleration(x, mass):
    diffs = np.subtract(x[:,np.newaxis,:], x[:,:,np.newaxis])
    pfactor = mass/np.sqrt(np.sum(diffs**2, axis=0) + epsilon)**3
    np.fill_diagonal(pfactor, 0.0)
    return np.sum(diffs*pfactor, axis=2)

def compute_energy(x, v, mass):
    # First the kinetic energy...
    T = compute_kinetic_energy(x,v,mass)
    # Now the potential energy
    U = compute_potential_energy(x,v,mass)
    return T - U

def compute_kinetic_energy(x, v, mass):
    return 0.5*np.sum(np.dot(v**2, mass))

def compute_potential_energy(x, v, mass):
    U = 0.0
    for i in range(x.shape[1]):
        diff = x[:,i][:,np.newaxis] - x[:,1+i:]
        delta = np.sum(diff**2, axis=0)
        U += mass[i]*np.sum(mass[1+i:]/np.sqrt(epsilon+delta))
    return U

def compute_center_of_mass(x, mass):
    return np.dot(x, mass) / np.sum(mass)

def center_particles(x, mass):
    x -= compute_center_of_mass(x, mass)[:,np.newaxis]

def integrate():
    # Assign initial values...
    randoms = np.random.randint(256**4, dtype='<u4', size=2*3*NP) / (256.0*256*256*256)

    # Initial position and speed
    widths = (L[1::2] - L[::2])[:,np.newaxis]
    x = (L[::2][:,np.newaxis] +
         widths * randoms[::2].reshape((NP,3)).transpose().copy())
    v = -0.2 + 0.4 * randoms[1::2].reshape((NP,3)).transpose().copy()

    # Assign random mass
    mass = drandom(low_mass,high_mass,NP)

    # Add a rotation around the z axis
    v[1,:] += x[0,:]/10.0
    v[0,:] -= x[1,:]/10.0

    if center_masses:
        # Set the center of mass and it's speed to 0
        center_particles(x, mass)
        center_particles(v, mass)

    if bounded_state:
        # Make sure that the total energy of the system is negative so particle don't fly in the distance
        # Set the kinetic energy to half the potential energy
        U = compute_potential_energy(x,v,mass)
        K = compute_kinetic_energy(x,v,mass)
        alpha = np.sqrt(U/(2.0*K))
        v *= alpha

    write_state(0,x)
    print(f"0.0  {compute_energy(x,v,mass)/NP:g}")

    if VERLET:
        acc = compute_acceleration(x,mass)
        for l in range(1,NT+1):
            # Print out the system's total energy per particle (should be fairly constant)
            if l%write_freq == 0:
                print(f"{dt*float(l):g}  {compute_energy(x,v,mass)/NP:g}")

            # Now update the arrays
            x += dt*v + 0.5*dt*dt*acc
            boundary_conditions(x)
            temp = compute_acceleration(x,mass)
            v += 0.5*dt*(acc + temp)
            acc = temp
            if l%write_freq == 0:
                write_state(l,x)
    else:
        # Fourth-order Runge-Kutta
        for l in range(1,NT+1):
            acc = compute_acceleration(x, mass)
            k1 = [v, acc]
            acc = compute_acceleration(x + 0.5*dt*k1[0], mass)
            k2 = [(1.0 + 0.5*dt)*v, acc]
            acc = compute_acceleration(x + 0.5*dt*k2[0], mass)
            k3 = [(1.0 + 0.5*dt + 0.25*dt*dt)*v, acc]
            acc = compute_acceleration(x + dt*k3[0], mass)
            k4 = [(1.0 + dt + 0.5*dt*dt + 0.25*dt*dt*dt)*v, acc]
            # Now update the arrays
            x += dt*(k1[0] + 2*k2[0] + 2*k3[0] + k4[0])/6.0
            v += dt*(k1[1] + 2*k2[1] + 2*k3[1] + k4[1])/6.0

            boundary_conditions(x)

            # Print out the system's total energy per particle (should be fairly constant)
            if l%write_freq == 0:
                print(f"{dt*float(l):g}  {compute_energy(x,v,mass)/NP:g}")
                write_state(l,x)

    write_state(NT,x)

    
def read_parameters(filename):
    global NT, NP, tvalue, seed, dt, epsilon, low_mass, high_mass, write_freq
    global finite_domain, center_masses, bounded_state, L
    try:
        s = open(filename, "r")
    except:
        # If the file doesn't exist, we need to exit...
        print(f"The file {filename} cannot be found!")
        sys.exit(1)

    # Loop through all lines in the parameter file
    param = dict(
        nparticle = str(NP),
        max_time = str(tvalue),
        seed = str(seed),
        timestep = str(dt),
        epsilon = str(epsilon),
        min_mass = str(low_mass),
        max_mass = str(high_mass),
        write_frequency = str(write_freq),
        finite_domain = "yes" if finite_domain else "no",
        center_of_mass = "yes" if center_masses else "no",
        bound_state = "yes" if bounded_state else "no",
        xmin = str(L[0]),
        xmax = str(L[1]),
        ymin = str(L[2]),
        ymax = str(L[3]),
        zmin = str(L[4]),
        zmax = str(L[5]),
    )
    for line in s:
        # If it's an empty line, or if the line begins with a #, or
        # if there's no equals sign in this line, ignore it
        if line != "\n" and line[0] != '#' and '=' in line:
            # Assumes that the equals sign can only occur once in 
            # the line
            name, value = map(str.strip, line.split('=', 1))
            param[name] = value

    # Now that we have the parameter name, see if it matches
    # any of the known parameters. If so, read in the value and
    # assign it
    NP = int(param["nparticle"])
    tvalue = float(param["max_time"])
    seed = int(param["seed"])
    dt = float(param["timestep"])
    epsilon = float(param["epsilon"])
    low_mass = float(param["min_mass"])
    high_mass = float(param["max_mass"])
    write_freq = int(param["write_frequency"])
    finite_domain = param["finite_domain"] == "yes"
    center_masses = param["center_of_mass"] == "yes"
    bounded_state = param["bound_state"] == "yes"
    L = np.array([float(param[name]) for name in ["xmin", "xmax", "ymin", "ymax", "zmin", "zmax"]])
    s.close()
    # Sanity checks
    assert tvalue > np.finfo(float).eps
    assert NP > 1
    assert dt > np.finfo(float).eps
    assert epsilon > np.finfo(float).eps and epsilon < 0.1
    assert write_freq > 0
    assert low_mass > np.finfo(float).eps
    assert high_mass >= low_mass
    assert seed >= 0
    for i in range(3):
        assert L[2*i+1] > L[2*i]
    if seed == 0:
        seed = None
    np.random.seed(seed)
    NT = int(tvalue/dt)

def main():
    if len(sys.argv) > 2:
        sys.stderr.write("Usage: ./nbody parameters.txt\n")
        sys.exit(0)

    if len(sys.argv) == 2:
        read_parameters(sys.argv[1])

    integrate()

main()
