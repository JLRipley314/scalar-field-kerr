#!/usr/bin/env python3
#=============================================================================
## for parsing command line arguments
#=============================================================================
import argparse
parser = argparse.ArgumentParser(
      description=
      "Launches simulation run, given parameters set in this file."
      )
args = parser.parse_args()
#=============================================================================
import time
from sim_class import Sim

sim= Sim()
#=============================================================================
## Set simulation parameters here 
#=============================================================================
sim.output_dir = None 
sim.bin_name= "default.run"
sim.recompile= False

sim.run_type= "basic_run"
sim.debug= False

sim.black_hole_mass= float(0.5)
sim.black_hole_spin= round(0.9*sim.black_hole_mass,16)
sim.compactification_length= float(50)

sim.evolve_time= float(50.0) ## units of black hole mass
sim.num_saved_times= int(200)

sim.nx= 160  ## number of radial pts 
sim.nl= 20   ## number of angular values
sim.nlat = 44 ## number of theta collocation points 
sim.nphi = 44 ## number of phi collocation points; must be a multiple of 4 

sim.constraint_damping = -100 ## damping of (q - \partial_r f)==0
#-----------------------------------------------------------------------------
## scalar field potential
#-----------------------------------------------------------------------------
sim.V2 = 0.0
sim.V3 = 0.0
sim.V4 = 0.0
#-----------------------------------------------------------------------------
## Initial data:
## l_ang:                  
##    Initial data is a particular swal function
## m_ang:                  
##    Must have m_ang >= 0  
## initial_data_direction: 
##    Which way pulse is approximately "heading"
## amp:
##    Initial amplitude of real/imaginary parts of psi4
## rl(ru)_0:               
##    Lower(upper) bounds of initial data as a multiple
##    of the black hole horizon
#-----------------------------------------------------------------------------
sim.l_ang=  4 
sim.m_ang=  0 
sim.amp= 0.01 
sim.rl_0= 1.1 
sim.ru_0= 20.0
sim.initial_data_type= "compact_pulse"
sim.initial_data_direction= "ingoing"
#-----------------------------------------------------------------------------
## details of computer setup 

sim.computer= "della"#"home"#
sim.out_stem= "/tigress/jripley/sf"

sim.num_threads= 20

## for cluster/slurm script

sim.walltime= "12:00:00" ## (hh:mm:ss)
sim.memory=   "2048" ## MB 
sim.email=    "lloydripley@gmail.com" ## for slurm notification
#=============================================================================
## different scripts for making convergence tests, etc.
#=============================================================================
if (sim.run_type == "basic_run"):
   sim.launch_run()
#-----------------------------------------------------------------------------
else:
   raise ValueError("run_type = "+str(sim.run_type)) 
