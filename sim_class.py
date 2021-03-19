#=============================================================================
import subprocess, os, sys, time, shutil
from typing import List 
from math import log

sys.path.insert(1,os.getcwd()+'/src/tables/')
#=============================================================================
from tables_cheb     import save_cheb
from tables_legendre import save_roots_weights_Legendre
from tables_swal     import save_Gauss_quad_vals_swaL 
#=============================================================================
class Sim:
#=============================================================================
   def __init__(self,args:List[str])->None:
      self.home_dir= str(os.getcwd())

      assert len(args) > 1, (
         'argv[1] is empty-meed a run_type to run!'
      )	
      self.run_type = args[1]
      if (len(args)>2 and args[2]=='debug'):
         self.debug=True
      else:
         self.debug=False
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
         self.output_dir= "output/"+self.output_stem
      elif (self.computer=="della"):
         self.output_dir=self.della_out_stem+self.output_stem
      else:
         raise ValueError("self.computer="+self.computer+" not supported")
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

      self.R_max= float(
         pow(self.compactification_length,2)
         /self.horizon
      )

      self.rl= [self.horizon*x for x in self.rl_0]
      self.ru= [self.horizon*x for x in self.ru_0]
#-----------------------------------------------------------------------------
      absa = abs(self.black_hole_spin/self.black_hole_mass)
      self.constraint_damping = abs(
         (1.0/self.black_hole_mass)*pow(abs(1.0-absa+1.0e-6),-0.5)
      )
#-----------------------------------------------------------------------------
## when to begin metric reconstruction
      self.integrate_psi4_start_time = self.black_hole_mass*max(
         [  self.integrate_psi4_start_multiple*(
               (2.0/self.black_hole_mass)*(x - self.horizon)
            +  4.0*log(max(x,self.horizon)/self.horizon)
            ) for x in self.ru
         ]
      )
      self.scd_order_start_time = self.black_hole_mass*max(
         [  self.scd_order_start_multiple*(
               (2.0/self.black_hole_mass)*(x - self.horizon)
            +  4.0*log(max(x,self.horizon)/self.horizon)
            ) for x in self.ru
         ]
      )
#-----------------------------------------------------------------------------
      self.lin_pos_m = [ m for m in self.lin_m if m>=0]
      self.scd_pos_m = [ m for m in self.scd_m if m>=0]
#-----------------------------------------------------------------------------
## for parallelism (if applied)
      self.num_threads= max(len(self.lin_m),len(self.scd_m))
#-----------------------------------------------------------------------------
## Gauss points for integration
## want to exactly integrate polynomials of order
## 2l + 2m(i.e. lmin) + alpha + beta (so being a bit conservative here) 
      lmin= max(max([abs(m) for m in self.lin_m]),abs(self.psi_spin))
      self.ny= (self.nl
      +	int(abs(2*lmin))
      + int(abs(2*max([abs(m) for m in self.lin_m])+self.psi_spin))
      +	int(abs(2*max([abs(m) for m in self.lin_m])-self.psi_spin))
      )
      if (self.ny%2!=0):
         self.ny+= 1
#-----------------------------------------------------------------------------
## consider characteristic speeds: it looks like max is ~2/3 so can multiply
## the factor 6/N^2 by up to (3/2) (need to experiment)
#-----------------------------------------------------------------------------
      self.dt= float(
         9.*pow(max(self.nx,self.ny),-2)
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
#-----------------------------------------------------------------------------
      self.max_m= max(
         max(self.lin_m),
         max(self.scd_m),
         1)
      self.min_m= min(
         min(self.lin_m),
         min(self.scd_m),
         1)
#-----------------------------------------------------------------------------
      self.max_s=  3
      self.min_s= -3
#-----------------------------------------------------------------------------
      assert(len(self.lin_m)==len(set(self.lin_m)))
      assert(len(self.scd_m)==len(set(self.scd_m)))
      assert(len(self.write_lin_m)==len(set(self.write_lin_m)))
      assert(len(self.write_scd_m)==len(set(self.write_scd_m)))
#=============================================================================
   def make_tables_dir(self)->None:
      self.tables_dir= self.output_dir+"/tables"
      os.makedirs(self.tables_dir)
#=============================================================================
   def write_sim_params(self)->None:
      with open(self.output_dir+'/sim_params.txt','w') as f:
         attrs= vars(self)
         for param in attrs:
            if type(attrs[param])!=list:
               f.write('{} {}\n'.format(param,attrs[param]))	
            else:
               f.write('len_'+param+' '+str(len(attrs[param]))+'\n')

               write_str= param
               for val in attrs[param]:
                  write_str+= ' '+str(val)
               f.write(write_str+'\n')
#=============================================================================
   def write_slurm_script(self):
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
         run_str= './bin/{} {}\n\n'.format(self.bin_name, self.output_dir)
         if (self.debug):
            run_str= 'valgrind -v --track-origins=yes --leak-check=full '+run_str
         f.write('\n'+run_str)

      shutil.copyfile(
      '{}/run.slurm'.format(self.home_dir),
      '{}/run.slurm'.format(self.output_dir)
      )
#=============================================================================
   def compile(self)->None:
      subprocess.call('make '+self.bin_name,shell=True)
#=============================================================================
   def launch_run(self)->None:
      self.set_derived_params()

      self.make_output_dir()

      self.make_tables_dir()
      save_cheb(self.tables_dir,self.nx)
      save_roots_weights_Legendre(self.tables_dir,self.ny)
      self.make_Gauss_pts()

      self.write_sim_params()
      if self.recompile==True:
         self.then_recompile()

      self.output_file= self.output_dir+'/output.txt'
      if (self.computer=='home'):
         run_str= (
            'time ./bin/'+self.bin_name+' '+self.output_dir+' > '+self.output_file+' 2>&1 &'
         )
         if (self.debug):
            run_str= 'valgrind -v --track-origins=yes --leak-check=full '+run_str
         os.environ['OMP_NUM_THREADS']= str(self.num_threads)
         subprocess.call(run_str,shell=True) 
      elif (self.computer=='della'):
         self.write_slurm_script()
         subprocess.call('sbatch run.slurm', shell='True')		
      else:
         raise ValueError('computer= '+self.computer+' not yet supported')
#=============================================================================
   def make_Gauss_pts(self)->None:
      for spin in [-3,-2,-1,0,1,2,3]:
         for m_ang in range(self.min_m,self.max_m+1):
            save_Gauss_quad_vals_swaL(self.tables_dir,spin,m_ang,self.nl,self.ny) 
#=============================================================================
   def then_recompile(self)->None:
      subprocess.call('make clean_obj'+self.bin_name,shell=True)
      subprocess.call('make '+self.bin_name,shell=True)
