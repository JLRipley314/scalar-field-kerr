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

sim.black_hole_mass= float(0.5)
sim.black_hole_spin= round(0.999999*sim.black_hole_mass,16)
sim.compactification_length= float(50)

sim.evolve_time= float(100.0) ## units of black hole mass
sim.num_saved_times= int(400)

sim.nx= 512  ## Number of radial pts 
sim.nl= 50   ## Number of angular l values
sim.nm= 40   ## Number of angular m values
sim.nlat = 104 ## Number of theta collocation points 
sim.nphi = 84 ## Number of phi collocation points; must be a multiple of 4 
#-----------------------------------------------------------------------------
## Scalar field potential
## V = (V2/2)*phi^2 + (V3/6)*phi^3 + (V4/24)*phi^4

sim.V2 = 0.5
sim.V3 = 0.0
sim.V4 = 1.0
#-----------------------------------------------------------------------------
## Initial data

sim.l_ang=  1 ## Initial data is a particular swal function.
sim.m_ang=  1 ## Must have m_ang >= 0.
sim.amp= 0.01 ## Initial amplitude of pulse 
sim.rl_0= 1.1 ## Lower(upper) bounds of initial data as a multiple
sim.ru_0= 4  ## of the black hole horizon.
sim.initial_data_type= "compact_pulse"
sim.initial_data_direction= "ingoing" 
#-----------------------------------------------------------------------------
## Details of computer setup 

sim.computer= "della"#"home"#
sim.out_stem= "/tigress/jripley/sf/highres"

sim.num_threads= 6 ## sets number of threads if using OpenMP 
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
