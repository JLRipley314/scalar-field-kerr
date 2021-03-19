!
! module for writing fields to file
!
!=============================================================================
module mod_write_level
!=============================================================================
   use, intrinsic :: iso_fortran_env, only: &
      stdout=>output_unit, stdin=>input_unit, stderr=>error_unit

   use mod_prec

   use mod_params, only: nx, ny 
   use mod_field,  only: field

   use mod_io, only: write_csv

   use mod_cheb, only: cheb_real_to_coef
   use mod_swal, only: swal_real_to_coef

   use mod_teuk, only: compute_res_q

   use mod_params, only: &
      write_lin_m, len_write_lin_m, &
      write_indep_res, &
      write_coefs, write_coefs_swal, constrained_evo

   use mod_fields_list, only: &
      psi4_lin_p, psi4_lin_q, psi4_lin_f, &
      res_lin_q, & 

      psi4_integrated_lin_f, &

      psi4_twice_integrated_lin_f
!=============================================================================
   implicit none
   private

   public :: write_out, write_diagnostics, write_level 

   interface write_horizon_or_scriplus
      module procedure &
         write_field_horizon_or_scriplus
   end interface
!=============================================================================
contains
!=============================================================================
   subroutine write_field_horizon_or_scriplus(time,location,m_ang,f)
      real(rp),     intent(in) :: time
      character(*), intent(in) :: location
      integer(ip),  intent(in) :: m_ang
      type(field),  intent(inout) :: f
      !----------------------------------------------
      character(:), allocatable :: fn
      complex(rp),  allocatable :: vals(:)

      fn = location//"_"//f%fname 

      !------------------------------------------------
      !------------------------------------------------
      select case (location)
      !------------------------------------------------
      case ("horizon")
         allocate(vals(lbound(f%np1,2):ubound(f%np1,2)))
         vals = f%np1(nx,:,m_ang) 
         call write_csv(fn, time, m_ang, vals)
      !------------------------------------------------
      case ("scriplus")
         allocate(vals(lbound(f%np1,2):ubound(f%np1,2)))
         vals = f%np1(1,:,m_ang) 
         call write_csv(fn, time, m_ang, vals)
      !------------------------------------------------
      case ("coef_horizon")
         allocate(vals(lbound(f%coefs_swal,2):ubound(f%coefs_swal,2)))
         call swal_real_to_coef( &
            f%spin, &
            m_ang, &
            f%np1, &
            f%coefs_swal &
         )
         vals = f%coefs_swal(nx,:,m_ang) 
         call write_csv(fn, time, m_ang, vals)
      !------------------------------------------------
      case ("coef_scriplus")
         allocate(vals(lbound(f%coefs_swal,2):ubound(f%coefs_swal,2)))
         call swal_real_to_coef( &
            f%spin, &
            m_ang, &
            f%np1, &
            f%coefs_swal &
         )
         vals = f%coefs_swal(1,:,m_ang) 
         call write_csv(fn, time, m_ang, vals)
      !------------------------------------------------
      case default
         ! do nothing
      !------------------------------------------------
      end select 
      !------------------------------------------------
      !----------------------------------------------
   end subroutine write_field_horizon_or_scriplus
!=============================================================================
! computes two norm
!-----------------------------------------------------------------------------
   subroutine write_norm(time,m_ang,f)
      real(rp),     intent(in) :: time
      integer(ip),  intent(in) :: m_ang
      type(field),  intent(in) :: f
      !----------------------------------------------
      character(:), allocatable :: fn
      real(rp)                  :: norm

      fn = "norm_"//f%fname 

      norm = sum(abs(conjg(f%np1(:,:,m_ang))*f%np1(:,:,m_ang)))

      norm = sqrt(norm / real(size(f%np1(:,:,m_ang)),kind=rp))

      call write_csv(fn, time, m_ang, norm)

   end subroutine write_norm
!=============================================================================
   subroutine write_out(time)
      real(rp), intent(in) :: time

      write(stdout,'(e14.6)') time

   end subroutine write_out
!=============================================================================
   subroutine write_diagnostics(time)
      real(rp), intent(in) :: time

      integer(ip) :: i 
      !-----------------------------------------------------------------------
      write(stdout,'(e14.6)') time
      !-----------------------------------------------------------------------
      ! field values at future null infinity and horizon
      !-----------------------------------------------------------------------
      do i=1,len_write_lin_m
         call write_horizon_or_scriplus(time,"horizon", write_lin_m(i),psi4_lin_f)
         call write_horizon_or_scriplus(time,"scriplus",write_lin_m(i),psi4_lin_f)
         call write_horizon_or_scriplus(time,"scriplus",write_lin_m(i),psi4_integrated_lin_f)
         call write_horizon_or_scriplus(time,"scriplus",write_lin_m(i),psi4_twice_integrated_lin_f)
      end do
      !-----------------------------------------------------------------------
      if (  (write_indep_res) &
      .and. (.not. constrained_evo) &
      ) then
         do i=1,len_write_lin_m
            call compute_res_q( write_lin_m(i),psi4_lin_q,psi4_lin_f,res_lin_q)
            call write_norm(time,write_lin_m(i),res_lin_q)
         end do
      end if
      !-----------------------------------------------------------------------
      if (write_coefs_swal) then
         !--------------------------------------------------------------------
         do i=1,len_write_lin_m
            call write_horizon_or_scriplus(time,"coef_horizon", write_lin_m(i),psi4_lin_f)
            call write_horizon_or_scriplus(time,"coef_scriplus",write_lin_m(i),psi4_lin_f)
            call write_horizon_or_scriplus(time,"coef_scriplus",write_lin_m(i),psi4_integrated_lin_f)
            call write_horizon_or_scriplus(time,"coef_scriplus",write_lin_m(i),psi4_twice_integrated_lin_f)
         end do 
         !--------------------------------------------------------------------
      end if

   end subroutine write_diagnostics
!=============================================================================
   subroutine write_level(time)
      real(rp), intent(in) :: time

      integer(ip) :: i
      !-----------------------------------------------------------------------
      ! \Psi_4^{(1)}
      !-----------------------------------------------------------------------
      do i=1,len_write_lin_m
         call write_csv(time,write_lin_m(i),psi4_lin_f)
      end do 
      !-----------------------------------------------------------------------
      if (write_indep_res) then
         if (.not. constrained_evo) then
            do i=1,len_write_lin_m
               call compute_res_q( write_lin_m(i),psi4_lin_q,psi4_lin_f,res_lin_q)
               call write_csv(time,write_lin_m(i),res_lin_q)
            end do
         end if
      end if
      !-----------------------------------------------------------------------
      if (write_coefs) then
         !--------------------------------------------------------------------
         do i=1,len_write_lin_m
            call cheb_real_to_coef( &
               write_lin_m(i), &
               psi4_lin_f%np1, &
               psi4_lin_f%coefs_cheb, &
               psi4_lin_f%re, &
               psi4_lin_f%im, &
               psi4_lin_f%coefs_cheb_re, &
               psi4_lin_f%coefs_cheb_im  &
            )
            call swal_real_to_coef( &
               psi4_lin_f%spin, &
               write_lin_m(i), &
               psi4_lin_f%coefs_cheb, &
               psi4_lin_f%coefs_both &
            )
            call write_csv( &
               "coefs_"//psi4_lin_f%fname, &
               time, &
               write_lin_m(i), &
               psi4_lin_f%coefs_both(:,:,write_lin_m(i)) &
            )
         end do 
      end if
      !-----------------------------------------------------------------------
      if (write_coefs_swal) then
         !--------------------------------------------------------------------
         do i=1,len_write_lin_m
            call swal_real_to_coef( &
               psi4_lin_f%spin, &
               write_lin_m(i), &
               psi4_lin_f%np1, &
               psi4_lin_f%coefs_swal &
            )
            call write_csv( &
               "swal_coefs_"//psi4_lin_f%fname, &
               time, &
               write_lin_m(i), &
               psi4_lin_f%coefs_swal(:,:,write_lin_m(i)) &
            )
         end do 
         !--------------------------------------------------------------------
      end if
      !-------------------------------------------------------------------- 
   end subroutine write_level
!=============================================================================
end module mod_write_level
