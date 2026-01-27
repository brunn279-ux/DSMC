module maxwellDist
   implicit none

   type :: particle
      real :: mass
      real :: x(3)
      real :: c(3)
   end type particle
   contains
      subroutine particleDist(n,m,V,T,p)
         type(particle), intent(in) :: p
         type(real), intent(in) :: m
         type(integer), intent(in) :: N
         type(real), intent(in) :: T
         type(real), intent(in) :: V

         integer :: i
         real :: k, pi, beta, R(6)

         k = 1.380649E-23
         pi = 3.14159265359

         call random_number(R)

         beta = (m/2/k/T)

      do i = 1, N
         p(i)%c(1) = V(1) + SIN(2*pi*R(1))*(-log(R(2)))**0.5
         p(i)%c(2) = V(2) + SIN(2*pi*R(4))*(-log(R(3)))**0.5
         p(i)%c(3) = V(3) + SIN(2*pi*R(5))*(-log(R(6)))**0.5
      end do
   end subroutine particleDist
end module maxwellDist
