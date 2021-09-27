# Problème à N corps

Dans le fichier [`nbody.cpp`](https://github.com/calculquebec/cq-formation-nbody/blob/main/nbody.cpp),
dans la [fonction `integrate()`](https://github.com/calculquebec/cq-formation-nbody/blob/main/nbody.cpp#L177),
des points dans un espace tridimensionnel sont créés aléatoirement à partir de la
[ligne 184](https://github.com/calculquebec/cq-formation-nbody/blob/main/nbody.cpp#L184).
Différentes propriétés sont initialisées : la position, la vitesse et la masse de chaque point.

L’intégration du modèle à *N* corps se fait selon deux algorithmes choisis au moment de la compilation :
[l'intégration de Verlet](https://fr.wikipedia.org/wiki/Int%C3%A9gration_de_Verlet)
ou la [méthode de Runge-Kutta d'ordre quatre](https://fr.wikipedia.org/wiki/M%C3%A9thodes_de_Runge-Kutta).
Dans un cas comme dans l’autre, chaque point est influencé par les autres points.
Par conséquent, chaque point subit une accélération, ce qui modifie sa position et sa vitesse.
Cette mise-à-jour globale du modèle se fait à chaque itération. Voir aussi :
* [Problème à *N* corps](https://fr.wikipedia.org/wiki/Probl%C3%A8me_%C3%A0_N_corps) (Wikipédia)

L’ensemble des paramètres du modèle sont dans un fichier
[`parameters.txt`](https://github.com/calculquebec/cq-formation-nbody/blob/main/parameters.txt).
Ce nom de fichier doit être donné sur la ligne de commande. Par exemple :
```
./nbody parameters.txt
```

Il y a aussi [deux fichiers MD5](https://github.com/calculquebec/cq-formation-nbody/tree/main/solutions/md5) pour la validation :
```
md5sum -c solutions/md5/verlet_parameters_txt.md5
md5sum -c solutions/md5/rkutta_parameters_txt.md5
```

---

# The N-Body Problem

In the file [`nbody.cpp`](https://github.com/calculquebec/cq-formation-nbody/blob/main/nbody.cpp),
you will find the [function `integrate()`](https://github.com/calculquebec/cq-formation-nbody/blob/main/nbody.cpp#L177),
you will see that an array of *N* points in Euclidean 3-space are created at random locations, starting on
[line 184](https://github.com/calculquebec/cq-formation-nbody/blob/main/nbody.cpp#L184).
In addition to the position, the velocity and the mass of each body are also initialized.

The integration of the *N* body model is carried out using two different algorithms which can be chosen
at the moment of compilation: the [Verlet integration](https://en.wikipedia.org/wiki/Verlet_integration)
and the [4th order Runge-Kutta method](https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods).
In either case, each body is subject to the influence of all of the remaining *N-1* bodies
and so undergoes an acceleration, which subsequently modifies its velocity and position.
This global update of the position of the *N* bodies occurs at each iteration. See also:
* [*n*-body problem](https://en.wikipedia.org/wiki/N-body_problem) (Wikipedia)

The set of parameters for the model are stored in the file
[`parameters.txt`](https://github.com/calculquebec/cq-formation-nbody/blob/main/parameters.txt).
This filename must be given on the command line as an argument to the binary. For example:
```
./nbody parameters.txt
```

There are also [two MD5 files](https://github.com/calculquebec/cq-formation-nbody/tree/main/solutions/md5) for validation:
```
md5sum -c solutions/md5/verlet_parameters_txt.md5
md5sum -c solutions/md5/rkutta_parameters_txt.md5
```

## Technical details

This *N* body gravitational code is written in C++, the
[`Makefile`](https://github.com/calculquebec/cq-formation-nbody/blob/main/Makefile) contains lines 
where you can specify the C++ compiler as well as any necessary compilation and 
linker flags. The resulting binary is called "`nbody`" and should be run with the 
command:
```
./nbody <parameter_file>
```

An example of a parameter file is provided
([`parameters.txt`](https://github.com/calculquebec/cq-formation-nbody/blob/main/parameters.txt)); if no parameter 
file is supplied to the program, it runs with the default parameter values 
that are shown in the [`global.h` header file](https://github.com/calculquebec/cq-formation-nbody/blob/main/global.h).
The program periodically writes out the particle positions to an MDL `.mol`
file that can be visualized using a free tool like Jmol or PyMOL. 

The program is based on a direct (particle-particle) method and uses either 
the velocity Verlet or the fourth-order Runge-Kutta algorithm to numerically 
integrate the first-order system

![equation](https://latex.codecogs.com/png.latex?x_i'=v_i)

![equation](https://latex.codecogs.com/png.latex?v_i%27=\sum_{j=1,j\ne{i}}^N%20m_j%20\frac{(x_i-x_j)}{(\epsilon+r_{ij}^2)^{3/2}})

where x_i and v_i are the position and velocity of the i-th particle, m_i its
mass and r_ij is the L2 distance between x_i and x_j. The force softening
parameter eps is included to avoid singularities caused by near collisions of 
the particles and should be small and positive. The integration algorithm 
used is controlled by the compilation flag VERLET. The particle motions are 
confined to a toroidal box when the parameter "finite_domain" is set to true. 
The console output consists of the time value and total energy per particle, 
written out every hundred timesteps. 
      
