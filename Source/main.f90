!
! Evolve linear Teukolsky field, reconstruct metric, and evolve
! second order Teukolsky field
!
!=============================================================================
program main
!=============================================================================
   use, intrinsic :: iso_fortran_env, only: stdout=>output_unit

   use mod_prec
   use mod_params, only: &
      read_params, &
      sparse_save, &
      nt, dt, t_step_save, black_hole_mass, &
      lin_m, len_lin_m, &
      psi_spin, psi_boost, &
      integrate_psi4_start_time

   use mod_field,        only: set_field, shift_time_step, time_integrate_field
   use mod_cheb,         only: cheb_init, cheb_filter, cheb_test
   use mod_swal,         only: swal_init, swal_filter, swal_test_orthonormal
   use mod_ghp,          only: ghp_init
   use mod_teuk,         only: teuk_init, teuk_time_step
   use mod_initial_data, only: set_initial_data
   use mod_bkgrd_np,     only: bkgrd_np_init
   use mod_write_level,  only: write_level, write_diagnostics

   use mod_fields_list, only: &
      psi4_lin_p, psi4_lin_q, psi4_lin_f, &
      res_lin_q, & 

      psi4_integrated_lin_f, &

      psi4_twice_integrated_lin_f

   implicit none
!=============================================================================
! Put everything in a block so valgrind doesn't get confused about 
! automatically deallocated memory
!=============================================================================
clean_memory: block
!=============================================================================
! declare and initialize variables, fields, etc.
!=============================================================================
   integer(ip) :: m_i, t_step
   real(rp)    :: time
!=============================================================================
   write (*,*) "Reading in params"
   call read_params()
!=============================================================================
   write (*,*) "Initializing fields"   
!-----------------------------------------------------------------------------
! first order metric field
!-----------------------------------------------------------------------------
   call set_field(fname="lin_p",spin=psi_spin,boost=psi_boost,falloff=1_ip,f=psi4_lin_p)
   call set_field(fname="lin_q",spin=psi_spin,boost=psi_boost,falloff=2_ip,f=psi4_lin_q)
   call set_field(fname="lin_f",spin=psi_spin,boost=psi_boost,falloff=1_ip,f=psi4_lin_f)

   call set_field(fname="integrated_lin_f", &
      spin=psi_spin,boost=psi_boost,falloff=1_ip,f=psi4_integrated_lin_f &
   )
   call set_field(fname="twice_integrated_lin_f", &
      spin=psi_spin,boost=psi_boost,falloff=1_ip,f=psi4_twice_integrated_lin_f &
   )
!-----------------------------------------------------------------------------
! independent residual fields
!-----------------------------------------------------------------------------
   call set_field(fname="res_lin_q",spin=-2_ip,boost=-2_ip,falloff=2_ip,f=res_lin_q)
!-----------------------------------------------------------------------------
! initialize chebyshev diff matrices, swal matrices, etc.
!-----------------------------------------------------------------------------
   call cheb_init()
   call swal_init()
   call ghp_init()
   call bkgrd_np_init()
   call teuk_init()
!=============================================================================
! initial data 
!=============================================================================
   write (stdout,*) "Setting up initial data"
!-----------------------------------------------------------------------------
   time = 0.0_rp

   do m_i=1,len_lin_m
      call set_initial_data(m_i, psi4_lin_p, psi4_lin_q, psi4_lin_f)
   end do

   if (sparse_save) then
      call write_diagnostics(time / black_hole_mass)
   else
      call write_level(      time / black_hole_mass)
   end if
!=============================================================================
! integrate in time 
!=============================================================================
   write (stdout,*) "Beginning time evolution"
!-----------------------------------------------------------------------------
   time_evolve: do t_step=1,nt
      time = t_step*dt
      !--------------------------------------------------------------------
      ! \Psi_4^{(1)} evolution 
      !--------------------------------------------------------------------
      !$OMP PARALLEL DO NUM_THREADS(len_lin_m) IF(len_lin_m>1)
      preserve_m_evo: do m_i=1,len_lin_m
         call teuk_time_step(lin_m(m_i),psi4_lin_p,psi4_lin_q,psi4_lin_f)
         !------------------------------------
         ! low pass filter (in spectral space)
         !------------------------------------
         call cheb_filter(lin_m(m_i),psi4_lin_p)
         call cheb_filter(lin_m(m_i),psi4_lin_q)
         call cheb_filter(lin_m(m_i),psi4_lin_f)

         call swal_filter(lin_m(m_i),psi4_lin_p)
         call swal_filter(lin_m(m_i),psi4_lin_q)
         call swal_filter(lin_m(m_i),psi4_lin_f)
         !------------------------------------
         if (time > integrate_psi4_start_time) then
            call time_integrate_field(lin_m(m_i),psi4_lin_f,           psi4_integrated_lin_f)
            call time_integrate_field(lin_m(m_i),psi4_integrated_lin_f,psi4_twice_integrated_lin_f)
         end if
      end do preserve_m_evo
      !$OMP END PARALLEL DO
      !-----------------------------------------------------------------------
      ! save to file 
      !-----------------------------------------------------------------------
      if (mod(t_step,t_step_save)==0) then
         if (sparse_save) then
            call write_diagnostics(time / black_hole_mass)
         else
            call write_level(      time / black_hole_mass)
         end if
      end if
      !-----------------------------------------------------------------------
      ! shift time steps
      !-----------------------------------------------------------------------
      do m_i=1,len_lin_m
         call shift_time_step(lin_m(m_i),psi4_lin_p)
         call shift_time_step(lin_m(m_i),psi4_lin_q)
         call shift_time_step(lin_m(m_i),psi4_lin_f)

         call shift_time_step(lin_m(m_i),psi4_integrated_lin_f)
         call shift_time_step(lin_m(m_i),psi4_twice_integrated_lin_f)
      end do

   end do time_evolve
!=============================================================================
   write (*,*) "Finished evolution"
!=============================================================================
end block clean_memory
!=============================================================================
end program main
