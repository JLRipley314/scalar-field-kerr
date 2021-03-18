#!/usr/bin/env python
#=============================================================================
## parameter file for evolution
## usage:
## ./setup.py [run_type] [debug]
## example:
## ./setup.py basic_run
#=============================================================================
import sys, time
from sim_class import Sim
#=============================================================================
args= sys.argv
sim= Sim(args)
#=============================================================================
sim.bin_name= 'default.run'
sim.recompile= False
#=============================================================================
sim.black_hole_mass= float(0.5)
sim.black_hole_spin= round(0.998*sim.black_hole_mass,16)
sim.compactification_length= float(1)
#=============================================================================
sim.evolve_time= float(250) ## units of black hole mass
sim.num_saved_times= int(1000)
#=============================================================================
sim.nx= 64 ## num radial pts 
sim.nl= 16 ## num l angular values
sim.nm= 16 ## num m angular values 
#=============================================================================
## whether or not only to save horizon/scriplus/norm or full field vals 
sim.sparse_save= True
#=============================================================================
sim.constrained_evo = False

sim.write_indep_res=  True
sim.write_coefs=      False
sim.write_coefs_swal= True
#=============================================================================
## details of computer setup 

sim.computer= 'della'#'home'#
sim.della_out_stem= '/tigress/jripley/tf-out/'

## for della cluster/slurm script

sim.walltime= '72:00:00' ## (hh:mm:ss)
sim.memory=   '2048' ## MB 
sim.email=    'lloydripley@gmail.com' ## for slurm notification
#=============================================================================
## Initial data
#=============================================================================
## l_ang:                  initial data is a particular swal function
## initial_data_direction: which way pulse is approximately "heading"
## amp_re(im):             initial amplitude of real/imaginary parts of psi4
## rl(ru)_0:               lower(upper) bounds of initial data as a multiple
##                         of the black hole horizon
#=============================================================================
## initial data for all linear modes 
#=============================================================================
sim.m_ang=  [-2,    2]
sim.l_ang=  [ 2,    2]
sim.amp_re= [ 0.1,  0.1]
sim.amp_im= [ 0.0,  0.0]
sim.rl_0=   [ 1.01, 1.01]
sim.ru_0=   [ 3.0 , 3.0]
sim.initial_data_direction= "ii"
#=============================================================================
## multiple of when can begin second order evolution 

sim.integrate_psi4_start_multiple= float(6.0)
#=============================================================================
if (sim.run_type == "basic_run"):
   sim.launch_run()
#=============================================================================
elif (sim.run_type == "res_study"):
   all_res= [
      [160,42],
      [172,46],
      [194,50]
   ]
  
   for i in range(len(all_res)):
      sim.nx= all_res[i][0]
      sim.nl= all_res[i][1]
      sim.launch_run()
#=============================================================================
elif (sim.run_type == "multiple"):
   all_spins= [
      round(0.998  *sim.black_hole_mass,16),
      round(0.99998*sim.black_hole_mass,16)
   ]
   all_res= [
      [160,42],
      [170,46],
      [180,50]
   ]
 
   for j in range(len(all_spins)): 
      sim.black_hole_spin= all_spins[j]
      for i in range(len(all_res)):
         sim.nx= all_res[i][0]
         sim.nl= all_res[i][1]
         sim.launch_run()
#=============================================================================
elif (sim.run_type == "amp_study"):
  
   sim.launch_run()

   for i in range(3):
      sim.amp_re= [2.0*x for x in sim.amp_re]
      sim.amp_im= [2.0*x for x in sim.amp_im]
      sim.launch_run()
#=============================================================================
elif (sim.run_type == "spin_ramp"):
   for bhs in [0,0.01,0.02,0.04,0.08,0.12,0.16,0.2,0.24,0.28,0.32]:
      sim.black_hole_spin= bhs
      sim.launch_run()
#=============================================================================
else:
   raise ValueError("run_type = "+str(sim.run_type)) 
