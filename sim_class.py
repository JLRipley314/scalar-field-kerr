#=============================================================================
import subprocess, os, sys, time, shutil
from typing import List 
from math import log
#=============================================================================
class Sim:
#=============================================================================
   def __init__(self)->None:
      self.home_dir= str(os.getcwd())
#=============================================================================
   def make_output_dir(self)->None:
      time_of_day= time.asctime().split()
      self.output_stem= str(
         time_of_day[1]
      +  '_'+time_of_day[2]
      +  '_'+time_of_day[3].replace(':','_')
      +	'_m'+str(self.black_hole_mass)
      +	'_s'+str(self.black_hole_spin)
      +	'_nx'+str(self.nx)
      +	'_nl'+str(self.nl)
      )
      if (self.computer=="home"):
         self.output_dir= "Output/"+self.output_stem
      else:
         self.output_dir=self.out_stem+self.output_stem
      os.makedirs(self.output_dir)
#=============================================================================
   def set_derived_params(self)->None:
      self.max_l = int(self.nl - 1)
#-----------------------------------------------------------------------------
## put the R_max at the location of the outer horizon.
## if near extremal limit then put R_max a m, which is in between
## the outer and inner horizons
      sqrt_term = pow(
         pow(self.black_hole_mass,2)
      -  pow(self.black_hole_spin,2)
      ,0.5)

      self.horizon= self.black_hole_mass+sqrt_term
      self.R_min = 0

      self.R_max= float(
         pow(self.compactification_length,2)
         /self.horizon
      )


      self.rl= self.horizon*self.rl_0
      self.ru= self.horizon*self.ru_0
#-----------------------------------------------------------------------------
      absa = abs(self.black_hole_spin/self.black_hole_mass)
      self.constraint_damping = abs(
         (1.0/self.black_hole_mass)*pow(abs(1.0-absa+1.0e-6),-0.5)
      )
#-----------------------------------------------------------------------------
## consider characteristic speeds: it looks like max is ~2/3 so can multiply
## the factor 6/N^2 by up to (3/2) (need to experiment)
#-----------------------------------------------------------------------------
      self.dt= float(
         9.*pow(max(self.nx,self.nl),-2)
      )
#-----------------------------------------------------------------------------
      self.nt= int(
         self.evolve_time*self.black_hole_mass/self.dt
      )
#-----------------------------------------------------------------------------
      self.t_step_save= int(
         self.nt/float(self.num_saved_times)
      )	
      if (self.t_step_save==0):
         self.t_step_save= 1
#=============================================================================
   def write_sim_params(self)->None:
      self.parameter_file = self.output_dir+'/params.txt'
      with open(self.parameter_file,'w') as f:
         attrs= vars(self)
         for param in attrs:
            f.write('{} {}\n'.format(param,attrs[param]))	
#=============================================================================
   def write_slurm_script(self, run_str:str):
      with open('{}/run.slurm'.format(self.home_dir), 'w') as f:
         f.write('#!/bin/sh\n')
         f.write('#SBATCH -J fteuk\t\t# job name\n')
         f.write('#SBATCH -t {}\t\t# walltime (dd:hh:mm:ss)\n'.format(self.walltime))
         f.write('#SBATCH -p physics\t\t# partition/queue name\n')
         f.write('#SBATCH --mem={}MB\n'.format(self.memory))
         f.write('#SBATCH --output={}\t\t# file for STDOUT\n'.format(self.output_file))
         f.write('#SBATCH --mail-user={}\t\t# Mail  id of the user\n'.format(self.email))
         #------------
         ## for openmp
         #------------
         f.write('#SBATCH --nodes=1\n')
         f.write('#SBATCH --ntasks-per-node=1\n')
         f.write('#SBATCH --cpus-per-task={}\n'.format(self.num_threads))
         f.write('\nexport OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK\n')
         #------------
         ## executable
         #------------
         f.write('\n'+run_str)

      shutil.copyfile(
      '{}/run.slurm'.format(self.home_dir),
      '{}/run.slurm'.format(self.output_dir)
      )
#=============================================================================
   def launch_run(self)->None:
      self.set_derived_params()
      
      if self.output_dir is None:
         self.make_output_dir()

      if self.recompile==True:
         self.then_recompile()

      self.write_sim_params()

      self.output_file= self.output_dir+'/out.txt'
      run_str= (
         './Bin/'+self.bin_name
         +' '+self.parameter_file
         +' '+self.output_dir
         +' | tee '+self.output_file#+' 2>&1 &'
      )
      if (self.debug):
         run_str= 'valgrind -v --track-origins=yes --leak-check=full --gen-suppressions=yes '+run_str
      if (self.computer=='home'):
         os.environ['OMP_NUM_THREADS']= str(self.num_threads)
         subprocess.call(run_str,shell=True) 
      else:
         self.write_slurm_script(run_str)
         subprocess.call('sbatch run.slurm', shell='True')		
#=============================================================================
   def then_recompile(self)->None:
      subprocess.call('make clean_obj'+self.bin_name,shell=True)
      subprocess.call('make '+self.bin_name,shell=True)
