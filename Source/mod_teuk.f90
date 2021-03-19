!
! Teukolsky field solver
!
module mod_teuk
!=============================================================================
   use mod_prec
   use mod_field, only: field, set_level
   use mod_params, only: &
      dt, nx, ny, max_l, &
      min_m, max_m, &
      spin=>psi_spin, & 
      cl=>compactification_length, &
      bhm=>black_hole_mass, &
      bhs=>black_hole_spin, &
      constraint_damping, &
      constrained_evo

   use mod_cheb, only: R=>Rarr, compute_DR
   use mod_swal, only: cy=>cyarr, sy=>syarr, compute_swal_laplacian
!=============================================================================
   implicit none
   private
   public :: teuk_init, teuk_time_step, compute_res_q 

   interface teuk_time_step 
      module procedure teuk_lin_time_step
   end interface
!=============================================================================
   real(rp), allocatable :: &
      A_pp(:,:,:), A_pq(:,:,:)

   complex(rp), allocatable :: &
      B_pp(:,:,:), B_pq(:,:,:), B_pf(:,:,:)

   real(rp), allocatable :: prefactor_p(:,:)

   complex(rp), parameter :: ZI = (0.0_rp, 1.0_rp)
!=============================================================================
contains
!=============================================================================
   subroutine teuk_init()
      integer(ip) :: m_ang

      allocate(A_pp(nx,ny,min_m:max_m))
      allocate(A_pq(nx,ny,min_m:max_m))

      allocate(B_pp(nx,ny,min_m:max_m))
      allocate(B_pq(nx,ny,min_m:max_m))
      allocate(B_pf(nx,ny,min_m:max_m))

      allocate(prefactor_p(nx,ny))
      !----------------------------
      m_loop: do m_ang=min_m,max_m
      !----------------------------
         A_pp(:,:,m_ang) = (2.0_rp / (cl**4)) * ( &
            (cl**6) &
         +  (cl**2)*(bhs**2 - 8*(bhm**2))*(R**2) &
         +  4*(bhs**2)*bhm*(R**3) &
         )

         A_pq(:,:,m_ang) = ((R**2) / (cl**4)) * ( &
            (cl**4) &
         -  2*(cl**2)*bhm*R &
         +  (bhs**2)*(R**2) &
         )
      !----------------------------
         B_pp(:,:,m_ang) = - ( &
            2*ZI*bhs*m_ang*((cl**2) + 4*bhm*R) / (cl**2) &
         -  2*(bhs**2)*(cl**2 + 6*bhm*R)*R / (cl**4)   &
         -  4*bhm*spin &
         +  8*(bhm**2)*(2 + spin)*R / (cl**2) &
         +  2*ZI*bhs*spin*cy &
         )

         B_pq(:,:,m_ang) = (2.0_rp * R / (cl**4)) * ( &
            2*(bhs**2)*(R**2) &
         +  (1 + spin)*(cl**4)  &
         -  (ZI*bhs*m_ang + (3 + spin)*bhm) * ((cl**2) * R) &
         )

         B_pf(:,:,m_ang) = - (2.0_rp * R / (cl**4)) * ( &
            ZI*bhs*(cl**2)*m_ang &
         -  (bhs**2)*R &
         +  (cl**2)*bhm*(1 + spin) &
         )
      end do m_loop
      !------------------------------------------------------------------
      prefactor_p = 1.0_rp / ( &
         8*bhm*(2*(cl**2)*bhm - (bhs**2)*R)*((cl**2) + 2*bhm*R) / (cl**4) &
      -  (bhs*sy)**2 &
      )

      do m_ang=min_m,max_m
         A_pp(:,:,m_ang) = A_pp(:,:,m_ang) * prefactor_p
         A_pq(:,:,m_ang) = A_pq(:,:,m_ang) * prefactor_p 

         B_pp(:,:,m_ang) = B_pp(:,:,m_ang) * prefactor_p 
         B_pq(:,:,m_ang) = B_pq(:,:,m_ang) * prefactor_p 
         B_pf(:,:,m_ang) = B_pf(:,:,m_ang) * prefactor_p 
      end do
      !------------------------------------------------------------------
   end subroutine teuk_init
!=============================================================================
   subroutine set_k(step, m_ang, & 
         p, q, f, &
         kp, kq, kf) 

      integer(ip), intent(in) :: step, m_ang 
      type(field), intent(inout) :: p, q, f
      complex(rp), dimension(nx,ny,min_m:max_m), intent(inout) :: kp, kq, kf 

      call set_level(step, m_ang, p)
      call set_level(step, m_ang, q)
      call set_level(step, m_ang, f)

      call compute_DR(step, m_ang, p)
      call compute_DR(step, m_ang, q)
      call compute_DR(step, m_ang, f)

      call compute_swal_laplacian(step, m_ang, f)

      !-------------------------------------
      kp(:,:,m_ang) = &
         A_pp(:,:,m_ang) * p%DR(:,:,m_ang) &
      +  A_pq(:,:,m_ang) * q%DR(:,:,m_ang) &

      +  B_pp(:,:,m_ang) * p%level(:,:,m_ang) &
      +  B_pq(:,:,m_ang) * q%level(:,:,m_ang) &
      +  B_pf(:,:,m_ang) * f%level(:,:,m_ang) &

      +  prefactor_p     * f%swal_lap(:,:,m_ang)
      !-------------------------------------
      if (.not. constrained_evo) then
         kq(:,:,m_ang) = &
            p%DR(:,:,m_ang) &
         -  constraint_damping * (q%level(:,:,m_ang) - f%DR(:,:,m_ang))
      end if
      !-------------------------------------
      kf(:,:,m_ang) = &
         p%level(:,:,m_ang)

   end subroutine set_k
!=============================================================================
! RK4 time integrator: linear teukolsky wave
!=============================================================================
   subroutine teuk_lin_time_step(m_ang, p, q, f) 
      integer(ip), intent(in)    :: m_ang
      type(field), intent(inout) :: p, q, f 
   !--------------------------------------------------------
   ! if first time then k1 has not been set from k5
   !--------------------------------------------------------
      if (f%first_time) then
         call set_k(1_ip, m_ang, p, q, f, p%k1, q%k1, f%k1) 

         p%first_time = .false.
         q%first_time = .false.
         f%first_time = .false.

         if (constrained_evo) then
            call compute_DR(1_ip, m_ang, f)
            q%n(:,:,m_ang) = f%DR(:,:,m_ang)
         end if
      end if

      p%l2(:,:,m_ang)= p%n(:,:,m_ang)+0.5_rp*dt*p%k1(:,:,m_ang)
      f%l2(:,:,m_ang)= f%n(:,:,m_ang)+0.5_rp*dt*f%k1(:,:,m_ang)

      if (constrained_evo) then
         call compute_DR(2_ip, m_ang, f)
         q%l2(:,:,m_ang) = f%DR(:,:,m_ang)
      else
         q%l2(:,:,m_ang)= q%n(:,:,m_ang)+0.5_rp*dt*q%k1(:,:,m_ang)
      end if
   !--------------------------------------------------------
      call set_k(2_ip, m_ang, p, q, f, p%k2, q%k2, f%k2) 

      p%l3(:,:,m_ang)= p%n(:,:,m_ang)+0.5_rp*dt*p%k2(:,:,m_ang)
      f%l3(:,:,m_ang)= f%n(:,:,m_ang)+0.5_rp*dt*f%k2(:,:,m_ang)

      if (constrained_evo) then
         call compute_DR(3_ip, m_ang, f)
         q%l3(:,:,m_ang) = f%DR(:,:,m_ang)
      else
         q%l3(:,:,m_ang)= q%n(:,:,m_ang)+0.5_rp*dt*q%k2(:,:,m_ang)
      end if
   !--------------------------------------------------------
      call set_k(3_ip, m_ang, p, q, f, p%k3, q%k3, f%k3) 

      p%l4(:,:,m_ang)= p%n(:,:,m_ang)+dt*p%k3(:,:,m_ang)
      f%l4(:,:,m_ang)= f%n(:,:,m_ang)+dt*f%k3(:,:,m_ang)

      if (constrained_evo) then
         call compute_DR(4_ip, m_ang, f)
         q%l4(:,:,m_ang) = f%DR(:,:,m_ang)
      else
         q%l4(:,:,m_ang)= q%n(:,:,m_ang)+dt*q%k3(:,:,m_ang)
      end if
   !--------------------------------------------------------
      call set_k(4_ip, m_ang, p, q, f, p%k4, q%k4, f%k4) 

      p%np1(:,:,m_ang)= p%n(:,:,m_ang) &
      +  (dt/6.0_rp)*(p%k1(:,:,m_ang)+2.0_rp*p%k2(:,:,m_ang)+2.0_rp*p%k3(:,:,m_ang)+p%k4(:,:,m_ang))

      f%np1(:,:,m_ang)= f%n(:,:,m_ang) &
      +  (dt/6.0_rp)*(f%k1(:,:,m_ang)+2.0_rp*f%k2(:,:,m_ang)+2.0_rp*f%k3(:,:,m_ang)+f%k4(:,:,m_ang))

      if (constrained_evo) then
         call compute_DR(5_ip, m_ang, f)
         q%np1(:,:,m_ang) = f%DR(:,:,m_ang)
      else
         q%np1(:,:,m_ang)= q%n(:,:,m_ang) &
         +  (dt/6.0_rp)*(q%k1(:,:,m_ang)+2.0_rp*q%k2(:,:,m_ang)+2.0_rp*q%k3(:,:,m_ang)+q%k4(:,:,m_ang))
      end if
   !------------------------------------------------------------
   ! want k5 for computing independent residuals
   !------------------------------------------------------------
      call set_k(5_ip, m_ang, p, q, f, p%k5, q%k5, f%k5) 
   end subroutine teuk_lin_time_step
!=============================================================================
! independent residula: q - \partial_R f
!=============================================================================
   subroutine compute_res_q(m_ang, q, f, res) 
      integer(ip), intent(in)    :: m_ang
      type(field), intent(in)    :: q 
      type(field), intent(inout) :: f
      type(field), intent(inout) :: res

      integer(ip), parameter :: step = 5_ip

      call compute_DR(step, m_ang, f)
      
      res%np1(:,:,m_ang) = f%DR(:,:,m_ang) - q%np1(:,:,m_ang)

   end subroutine compute_res_q
!=============================================================================
end module mod_teuk

