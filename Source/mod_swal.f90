!
! Spin-weighted associated Legendre function routines, along with
! Legendre polynomial roots and weights for Gaussian qudrature to
! go to/from coefficient and real space
!
module mod_swal 
!=============================================================================
   use mod_prec
   use mod_field,  only: field, set_level
   use mod_io,     only: set_arr
   use mod_params, only: &
      nx, ny, nl, max_l, &
      min_m, max_m, min_s, max_s

   implicit none
!=============================================================================
   private

   ! gauss points y, cos(y), sin(y)
   real(rp), allocatable, protected, public :: Yvec(:),   cyvec(:),   syvec(:) 
   real(rp), allocatable, protected, public :: Yarr(:,:), cyarr(:,:), syarr(:,:)

   ! subroutines 
   public :: swal_init, compute_swal_laplacian, swal_lower, swal_raise

   public :: swal_filter

   public :: swal_test_orthonormal, swal_test_to_from

   public :: swal_real_to_coef

   ! weights for Gaussian integration 
   real(rp), allocatable :: weights(:)

   ! going to/from swal sapce
   real(rp), allocatable, protected, public :: swal(:,:,:,:)

   ! evaluate in real space 
   real(rp), allocatable :: &
      laplacian(:,:,:,:), lower(:,:,:,:), raise(:,:,:,:), low_pass(:,:,:,:) 
!=============================================================================
contains
!=============================================================================
   integer(ip) function compute_lmin(spin, m_ang) result(lmin)
      integer(ip), intent(in)  :: spin, m_ang

      lmin = max(abs(spin),abs(m_ang))
   end function compute_lmin
!=============================================================================
   subroutine swal_init()
      character(:), allocatable :: mstr, sstr
      integer(ip) :: m_ang, spin, i, j, k, min_l
      !-----------------------------------------------------------------
      allocate(Yvec(ny))
      allocate(cyvec(ny))
      allocate(syvec(ny))
      allocate(Yarr(nx,ny))
      allocate(cyarr(nx,ny)) 
      allocate(syarr(nx,ny)) 

      allocate(weights(ny))

      allocate(swal(ny,0:max_l,min_m:max_m,min_s:max_s))
      allocate(laplacian(ny,ny,min_m:max_m,min_s:max_s)) 
      allocate(lower(    ny,ny,min_m:max_m,min_s:max_s))
      allocate(raise(    ny,ny,min_m:max_m,min_s:max_s))
      allocate(low_pass( ny,ny,min_m:max_m,min_s:max_s))
      !-----------------------------------------------------------------
      call set_arr('roots_legendre.txt', ny, Yvec)

      call set_arr('cos.txt', ny, cyvec)
      call set_arr('sin.txt', ny, syvec)

      do i=1,nx
         Yarr(i,:)  = Yvec
         cyarr(i,:) = cyvec
         syarr(i,:) = syvec
      end do
      !-----------------------------------------------------------------
      ! weights for Gaussian quadrature
      call set_arr('weights_legendre.txt', ny, weights)
      !-----------------------------------------------------------------
      ! spin-weighted spherical associated Legendre polynomials 
      do spin=min_s,max_s

         sstr = '     '
         write (sstr,'(i5)') spin
         sstr = trim(adjustl(sstr))

         do m_ang=min_m,max_m

            mstr = '     '
            write (mstr,'(i5)') m_ang
            mstr = trim(adjustl(mstr))

            call set_arr('s_'//sstr//'_m_' //mstr//'.txt', ny,nl,swal(:,:,m_ang,spin))
         end do
      end do
      !-----------------------------------------------------------------
      ! precompute swal laplacian array
      !-----------------------------------------------------------------
      laplacian = 0
      do spin =min_s,max_s
      do m_ang=min_m,max_m
         do i=1,ny
         do j=1,ny
            min_l = compute_lmin(spin, m_ang)

            do k=min_l,max_l
               laplacian(j,i,m_ang,spin) = laplacian(j,i,m_ang,spin) &
               -  real(k-spin,rp)*real(k+spin+1,rp)   &
                  *weights(j)*swal(j,k,m_ang,spin)  &
                  *swal(i,k,m_ang,spin)
            end do
         end do
         end do
      end do
      end do
      !-----------------------------------------------------------------
      ! precompute raising array
      !-----------------------------------------------------------------
      raise = 0
      do spin =min_s,max_s-1
      do m_ang=min_m,max_m
         do i=1,ny
         do j=1,ny
            min_l = max(abs(spin), abs(spin+1), abs(m_ang))

            do k=min_l,max_l
               raise(j,i,m_ang,spin) = raise(j,i,m_ang,spin) &
               +  sqrt(real(k-spin,rp)*real(k+spin+1.0_rp,rp)) &
                  *weights(j)*swal(j,k,m_ang,spin)  &
                  *swal(i,k,m_ang,spin+1)
            end do
         end do
         end do
      end do
      end do
      !-----------------------------------------------------------------
      ! precompute lowering array
      !-----------------------------------------------------------------
      lower = 0
      do spin =min_s+1,max_s
      do m_ang=min_m,  max_m
         do i=1,ny
         do j=1,ny
            min_l = max(abs(spin), abs(spin-1), abs(m_ang))

            do k=min_l,max_l
               lower(j,i,m_ang,spin) = lower(j,i,m_ang,spin) &
               -  sqrt(real(k+spin,rp)*real(k-spin+1,rp))   &
                  *weights(j)*swal(j,k,m_ang,spin)  &
                  *swal(i,k,m_ang,spin-1)
            end do
         end do
         end do
      end do
      end do
      !-----------------------------------------------------------------
      ! precompute swal low-pass filter array
      !-----------------------------------------------------------------
      low_pass = 0
      do spin =min_s,max_s
      do m_ang=min_m,max_m
         do i=1,ny
         do j=1,ny
            do k=0,max_l
               low_pass(j,i,m_ang,spin) = low_pass(j,i,m_ang,spin) &
               +  exp(-40.0_rp*(real(k,rp)/real(max_l,rp))**16) &
                  *weights(j)*swal(j,k,m_ang,spin)  &
                  *swal(i,k,m_ang,spin)
            end do
         end do
         end do
      end do
      end do
      !-----------------------------------------------------------------
      !-----------------------------------------------------------------
   end subroutine swal_init
!=============================================================================
! Gaussian quadrature.
! Integrate against complex conjugate cc[s_Y^m_l]=(-1)^{m+s} {-s}_Y^{-m}_l
!=============================================================================
   subroutine swal_real_to_coef(spin,m_ang,vals,coefs)
      integer(ip), intent(in) :: spin
      integer(ip), intent(in) :: m_ang
      complex(rp), dimension(nx,     ny,min_m:max_m), intent(in)  :: vals
      complex(rp), dimension(nx,0:max_l,min_m:max_m), intent(out) :: coefs

      integer(ip) :: lmin, j, k

      lmin  = compute_lmin(spin,m_ang)

      coefs(:,:,m_ang) = 0.0_rp

      do k=lmin,max_l
      do j=1,ny
         coefs(:,k,m_ang) =  &
            coefs(:,k,m_ang) &
         +  (vals(:,j,m_ang) * weights(j) * swal(j,k,m_ang,spin))
      end do
      end do
   end subroutine swal_real_to_coef
!=============================================================================
! coefficient synthesis
!=============================================================================
   subroutine swal_coef_to_real(spin,m_ang,coefs,vals)
      integer(ip), intent(in) :: spin
      integer(ip), intent(in) :: m_ang
      complex(rp), dimension(nx,0:max_l,min_m:max_m), intent(in)  :: coefs
      complex(rp), dimension(nx,     ny,min_m:max_m), intent(out) :: vals

      integer(ip) :: lmin, j, k

      lmin = compute_lmin(spin,m_ang)
      vals(:,:,m_ang) = 0.0_rp

      do j=1,ny
      do k=lmin,max_l
         vals(:,j,m_ang) = &
            vals(:,j,m_ang) &
         + (coefs(:,k,m_ang) * swal(j,k,m_ang,spin))
      end do
      end do
   end subroutine swal_coef_to_real
!=============================================================================
   subroutine swal_laplacian(spin,m_ang,vals,vals_lap)
      integer(ip), intent(in) :: spin
      integer(ip), intent(in) :: m_ang
      complex(rp), dimension(nx,ny,min_m:max_m), intent(in)  :: vals
      complex(rp), dimension(nx,ny,min_m:max_m), intent(out) :: vals_lap

      integer(ip) :: j, jp

      vals_lap(:,:,m_ang) = 0

      do j =1,ny
      do jp=1,ny
         vals_lap(:,j,m_ang) = vals_lap(:,j,m_ang) &
         +  vals(:,jp,m_ang)*laplacian(jp,j,m_ang,spin)
      end do
      end do

   end subroutine swal_laplacian
!=============================================================================
   subroutine compute_swal_laplacian(step,m_ang,f)
      integer(ip), intent(in) :: step, m_ang
      type(field), intent(inout) :: f

      call set_level(step,m_ang,f)

      call swal_laplacian(f%spin,m_ang,f%level,f%swal_lap)

   end subroutine compute_swal_laplacian
!=============================================================================
   subroutine swal_lower(spin,m_ang,vals,vals_lowered)
      integer(ip), intent(in) :: spin
      integer(ip), intent(in) :: m_ang
      complex(rp), dimension(nx,ny,min_m:max_m), intent(in)  :: vals
      complex(rp), dimension(nx,ny,min_m:max_m), intent(out) :: vals_lowered 

      integer(ip) :: j, jp

      vals_lowered(:,:,m_ang) = 0

      do j =1,ny
      do jp=1,ny
         vals_lowered(:,j,m_ang) = vals_lowered(:,j,m_ang) &
         +  vals(:,jp,m_ang)*lower(jp,j,m_ang,spin)
      end do
      end do

   end subroutine swal_lower
!=============================================================================
   subroutine swal_raise(spin,m_ang,vals,vals_raised)
      integer(ip), intent(in) :: spin
      integer(ip), intent(in) :: m_ang
      complex(rp), dimension(nx,ny,min_m:max_m), intent(in)  :: vals
      complex(rp), dimension(nx,ny,min_m:max_m), intent(out) :: vals_raised

      integer(ip) :: j, jp

      vals_raised(:,:,m_ang) = 0

      do j =1,ny
      do jp=1,ny
         vals_raised(:,j,m_ang) = vals_raised(:,j,m_ang) &
         +  vals(:,jp,m_ang)*raise(jp,j,m_ang,spin)
      end do
      end do

   end subroutine swal_raise
!=============================================================================
! Low pass filter. A smooth filter appears to help prevent Gibbs-like ringing
!=============================================================================
   subroutine swal_filter(m_ang,f)
      integer(ip), intent(in) :: m_ang
      type(field), intent(inout) :: f

      integer(ip) :: j, jp

      f%level(:,:,m_ang) = 0

      do j =1,ny
      do jp=1,ny
         f%level(:,j,m_ang) = f%level(:,j,m_ang) &
         +  f%np1(:,jp,m_ang)*low_pass(jp,j,m_ang,f%spin)
      end do
      end do

      f%np1(:,:,m_ang) = f%level(:,:,m_ang)

   end subroutine swal_filter
!=============================================================================
! test that the swal are orthogonal
!=============================================================================
   subroutine swal_test_orthonormal()

      integer(ip) :: s, m, l1, l2, j
      real(rp)    :: integral

      write (*,*) 'swal_test_orthonormal'

      do s=min_s,max_s
      do m=min_m,max_m
      do l1= 0,max_l
      do l2=l1,max_l
         integral = 0.0_rp
         do j=1,ny 
            integral = integral + weights(j)*swal(j,l1,m,s)*swal(j,l2,m,s)!*(((-1.0_rp)**(m+s))*swal(j,l2,-m,-s))
         end do
         if (abs(integral)>1e-14) then
            write (*,*) s,m,l1,l2,integral
         end if
      end do 
      end do 
      end do 
      end do 

   end subroutine swal_test_orthonormal
!=============================================================================
   subroutine swal_test_to_from(spin,m_ang,vals_in,coefs,vals_out)
      integer(ip), intent(in) :: spin, m_ang
      complex(rp), dimension(nx,ny,     min_m:max_m), intent(in)    :: vals_in
      complex(rp), dimension(nx,0:max_l,min_m:max_m), intent(inout) :: coefs
      complex(rp), dimension(nx,ny,     min_m:max_m), intent(out)   :: vals_out

      write(*,*) "swal_test_to_from. spin: ", spin, ", m_ang: ", m_ang

      call swal_real_to_coef(spin,m_ang,vals_in,coefs)
      call swal_coef_to_real(spin,m_ang,coefs,vals_out)

   end subroutine swal_test_to_from
!=============================================================================
end module mod_swal
