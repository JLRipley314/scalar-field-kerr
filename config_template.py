#!/usr/bin/env python3
##
## Launches simulation run, given parameters set in this file.
##
from sim_class import Sim
sim= Sim()
#=============================================================================
## Set simulation parameters: 

sim.output_dir = None 
sim.recompile= False

sim.run_type= 'basic_run'

sim.black_hole_mass= float(1)
sim.black_hole_spin= round(0.9999*sim.black_hole_mass,16)
sim.compactification_length= float(20)

sim.evolve_time= float(600.0) ## units of black hole mass
sim.num_saved_times= int(1200)

sim.nl= 20 ## Number of angular l values
sim.nm= 18 ## Number of angular m values
sim.nlat = 40 ## Number of theta collocation points 
sim.nphi = 40 ## Number of phi collocation points; must be a multiple of 4 

sim.nxs=   [8193]#, 33] ## Number of radial pts on each domain 
sim.rvals= [1]#,   4] ## Radial boundary of each domain (multiples of r_h) 

sim.use_cheb=False ## Radial derivatives: Chebyshev vs. finite differences
if (sim.use_cheb):
   sim.bin_name= 'default_cheb.run'
else:
   sim.cfl= 0.75 ## CFL number
   sim.bin_name= 'default_fd.run'
#-----------------------------------------------------------------------------
## Prefactor to kinetric term.
## In the action:
## (km1/psi + k0 + k1*psi + k2*psi^2)*(1/2)*(\nabla\psi)^2
sim.km1 = 0.0
sim.k0  = 1.0
sim.k1  = 0.0
sim.k2  = 0.0

## Scalar field potential
## V = (V2/2)*psi^2 + (V3/6)*psi^3 + (V4/24)*psi^4

sim.V2 = 1.0
sim.V3 = 0.0
sim.V4 = 100.0
#-----------------------------------------------------------------------------
## Initial data

sim.l_ang=  0 ## Initial data is a particular swal function.
sim.m_ang=  0 ## Must have m_ang >= 0.
sim.amp= 0.5 ## Initial amplitude of pulse 
sim.rl_0= -1.0 ## Lower(upper) bounds of initial data as a multiple
sim.ru_0= 3.0 ## of the black hole horizon.
sim.initial_data_type= 'compact_pulse'
sim.initial_data_direction= 'outgoing' 
#-----------------------------------------------------------------------------
## Details of computer setup 

sim.computer= 'della'#'home'#
sim.out_stem= '/tigress/jripley/sf/fd'

sim.time_it= True ## whether or not time run 

sim.num_omp_levels= 1 ## for nested OpenMP parallelism
sim.num_threads= 4 ## sets number of threads if using OpenMP 
#-----------------------------------------------------------------------------
## For writing slurm script

sim.walltime= '24:00:00' ## (hh:mm:ss)
sim.memory=   '8192'#'4096'#'2048' ## MB 
sim.email=    'lloydripley@gmail.com' ## For slurm notification
#=============================================================================
## Different scripts for making convergence tests, etc.

if (sim.run_type == 'basic_run'):
   sim.launch_run()
#-----------------------------------------------------------------------------
else:
   raise ValueError('run_type = '+str(sim.run_type)) 
