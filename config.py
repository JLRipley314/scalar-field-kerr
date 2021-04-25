#!/usr/bin/env python3
##
## Launches simulation run, given parameters set in this file.
##
from sim_class import Sim
sim= Sim()
#=============================================================================
## Set simulation parameters: 

sim.output_dir = None 
sim.bin_name= "default.run"
sim.recompile= False

sim.run_type= "basic_run"
sim.debug= False

sim.black_hole_mass= float(1)
sim.black_hole_spin= round(0.9*sim.black_hole_mass,16)
sim.compactification_length= float(10)

sim.evolve_time= float(5.0) ## units of black hole mass
sim.num_saved_times= int(50)

sim.nx= 64 ## Number of radial pts 
sim.nl= 16 ## Number of angular l values
sim.nm= 12 ## Number of angular m values
sim.nlat = 32 ## Number of theta collocation points 
sim.nphi = 28 ## Number of phi collocation points; must be a multiple of 4 

sim.nxs=   [64, 64, 64] ## Number of radial pts on each domain 
sim.Rvals= [0, 1, 2, 3] ## Radial Boundary of each domain 
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

sim.V2 = 0.0
sim.V3 = 0.0
sim.V4 = 0.0
#-----------------------------------------------------------------------------
## Initial data

sim.l_ang=  1 ## Initial data is a particular swal function.
sim.m_ang=  0 ## Must have m_ang >= 0.
sim.amp= 0.01 ## Initial amplitude of pulse 
sim.rl_0= -2 ## Lower(upper) bounds of initial data as a multiple
sim.ru_0=  2 ## of the black hole horizon.
sim.initial_data_type= "compact_pulse"
sim.initial_data_direction= "ingoing" 
#-----------------------------------------------------------------------------
## Details of computer setup 

sim.computer= "home"#"della"#
sim.out_stem= "/tigress/jripley/sf/highres"

sim.num_threads= 2 ## sets number of threads if using OpenMP 
#-----------------------------------------------------------------------------
## For writing slurm script

sim.walltime= "72:00:00" ## (hh:mm:ss)
sim.memory=   "2048" ## MB 
sim.email=    "lloydripley@gmail.com" ## For slurm notification
#=============================================================================
## Different scripts for making convergence tests, etc.

if (sim.run_type == "basic_run"):
   sim.launch_run()
#-----------------------------------------------------------------------------
else:
   raise ValueError("run_type = "+str(sim.run_type)) 
