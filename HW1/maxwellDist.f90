module maxwellDist

   implicit none

   type :: particle
      real :: mass
      real :: x(3)
      real :: c(3)
   end type particle

   interface
      module function velDist(m,T,V,N,p)
         type(particle), intent(in) :: p
         type(real), intent(in) :: m
         type(integer), intent(in) :: N
         type(real), intent(in) :: T
         type(real), intent(in) :: V
      end function velDist
   end interface
end module maxwellDist

submodule (maxwellDist) genVel
   contains
      module procedure velDist
      integer :: i
      real :: k, pi, R1, R2

      k = 1.380649E-23
      pi = 3.14159265359

      call random_number(R1)
      call random_number(R2)

      beta = (m/2/k/T)

      do i = 1, N
         p(i)%c(1) = V(1) + SIN(2*pi*R1)*(-log(R2))**0.5
         p(i)%c(2) = V(2) + SIN(2*pi*R1)*(-log(R2))**0.5
         p(i)%c(3) = V(3) + SIN(2*pi*R1)*(-log(R2))**0.5
      end do
      end procedure velDist
end submodule genVel
