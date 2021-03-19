!
! Newman-Penrose fields of definite spin and boost weight 
! for the background Kerr spacetime.
! The fields have been rescaled by the highest power of R that allows
! them to remain regular and nonzero from the black hole horizon to
! future null infinity.
!
! rh_0:   rescaled by R
! mu_0:   rescaled by R
! ta_0:   rescaled by R^2
! pi_0:   rescaled by R^2
! ep_0:   rescaled by R^2
! psi2_0: rescaled by R^3
!
!=============================================================================
module mod_bkgrd_np
!=============================================================================
   use mod_prec
   use mod_params, only: &
      ZI, nx, ny, & 
      cl=>compactification_length, &
      bhm=>black_hole_mass, &
      bhs=>black_hole_spin

   use mod_cheb, only: R=>Rarr
   use mod_swal, only: cy=>cyarr, sy=>syarr
!=============================================================================
   implicit none
   private

   public :: bkgrd_np_init

   complex(rp), allocatable, public, protected :: &
      mu_0(:,:), ta_0(:,:), pi_0(:,:), rh_0(:,:), ep_0(:,:), psi2_0(:,:)  
!=============================================================================
contains
!=============================================================================
   subroutine bkgrd_np_init()
      allocate(mu_0(nx,ny))
      allocate(ta_0(nx,ny))
      allocate(pi_0(nx,ny))
      allocate(rh_0(nx,ny))
      allocate(ep_0(nx,ny))
      allocate(psi2_0(nx,ny))
      !----------------------------
      mu_0 = &
         1.0_rp / (-(cl**2) + ZI*bhs*R*cy)
      !----------------------------
      ta_0 = &
         (ZI*bhs*sy/sqrt(2.0_rp)) / ((cl**2 - ZI*bhs*R*cy)**2)
      !----------------------------
      pi_0 = &
      -  (ZI*bhs*sy/sqrt(2.0_rp)) / (cl**4 + (bhs*cy*R)**2)
      !----------------------------
      rh_0 = &
      -  0.5_rp*( &
            cl**4 - 2.0_rp*(cl**2)*bhm*R + (bhs*R)**2 &
         )/( &
            ((cl**2 - ZI*bhs*R*cy)**2)*(cl**2 + ZI*bhs*R*cy) &
         )
      !----------------------------
      ep_0 = &
         0.5_rp*( &
            (cl**2)*bhm - (bhs**2)*R - ZI*bhs*(cl**2-bhm*R)*cy &
         )/( &
            ((cl**2 - ZI*bhs*R*cy)**2)*(cl**2 + ZI*bhs*R*cy) &
         )
      !----------------------------
      psi2_0 = &
      -  bhm / ((cl**2 - ZI*bhs*R*cy)**3)
      !----------------------------
   end subroutine bkgrd_np_init
!=============================================================================
end module mod_bkgrd_np
