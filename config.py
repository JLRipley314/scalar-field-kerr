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
sim.black_hole_spin= round(0.998*sim.black_hole_mass,16)
sim.compactification_length= float(1)

sim.evolve_time= float(250) ## units of black hole mass
sim.num_saved_times= int(1000)

sim.nx= 64 ## num radial pts 
sim.nl= 16 ## num angular values
#-----------------------------------------------------------------------------
## Initial data:
##
## l_ang:                  initial data is a particular swal function
## initial_data_direction: which way pulse is approximately "heading"
## amp:                    initial amplitude of real/imaginary parts of psi4
## rl(ru)_0:               lower(upper) bounds of initial data as a multiple
##                         of the black hole horizon
#-----------------------------------------------------------------------------
sim.l_ang=  2 
sim.amp= 0.01 
sim.rl_0= 4.2 
sim.ru_0= 5.8
sim.initial_data_direction= "i"
#-----------------------------------------------------------------------------
## details of computer setup 

sim.computer= "home"#"della"#
sim.out_stem= "/tigress/jripley/tf-out/"

sim.num_threads= 1

## for cluster/slurm script

sim.walltime= "72:00:00" ## (hh:mm:ss)
sim.memory=   "2048" ## MB 
sim.email=    "lloydripley@gmail.com" ## for slurm notification
#-----------------------------------------------------------------------------
## multiple of when can begin second order evolution 

sim.integrate_psi4_start_multiple= float(6.0)
sim.scd_order_start_multiple= float(6.0)
#=============================================================================
## different scripts for making convergence tests, etc.
#=============================================================================
if (sim.run_type == "basic_run"):
   sim.launch_run()
#-----------------------------------------------------------------------------
else:
   raise ValueError("run_type = "+str(sim.run_type)) 
