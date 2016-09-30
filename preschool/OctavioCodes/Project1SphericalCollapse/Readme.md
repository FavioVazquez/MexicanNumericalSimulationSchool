COLD COLLAPSE

1.- Generate a homogeneous sphere with N particles

- Select a random number generator (tipically the generate homogeneous numbers
between 0-1)

-  Generate coordinates x,y,z and that belong to a sphere with radius = 1 centered at the
origin. Assume vx=vy=vz=0 and with the same mass = 1/Nparticles

- Load the output into tipsy or other visualization tool. Measure the mass density profile
with the function profile if you selected tipsy.

- Evolve the evolution with nbody1.f. You have to choose softening, number of
particles,time step and total time of integration.

- Study energy conservation and evolution of the system. Discuss. (Make plots)

- Parallelize the initial condition generator with OMP.
