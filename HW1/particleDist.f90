module maxwellDist
   implicit none

   type :: particle
      real :: mass
      real :: x(3)
      real :: c(3)
   end type particle
   contains
      subroutine particleDist(m,V,T,p)
         type(particle), intent(inout) :: p
         real, intent(in) :: m
         real, intent(in) :: T
         real, intent(in) :: V(3)

         integer :: i
         real :: k, pi, beta, R(6)

         k = 1.380649E-23
         pi = 3.14159265359

         call random_number(R)

         beta = (m/2/k/T)

         p%c(1) = V(1) + SIN(2*pi*R(1))*(-log(R(2)))**0.5
         p%c(2) = V(2) + SIN(2*pi*R(4))*(-log(R(3)))**0.5
         p%c(3) = V(3) + SIN(2*pi*R(5))*(-log(R(6)))**0.5
   end subroutine particleDist
end module maxwellDist
