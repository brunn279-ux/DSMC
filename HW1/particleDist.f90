program makeVelocities

   !input of number of particles bulk velocity, etc
   integer :: N
   real(8) :: mass, T, cx, cy, cz
   character(len=48) :: args
   type(particle) :: p(N)
   real(8) :: V(3)
   integer :: i
   integer :: io

   N = 1
   mass = 1.0
   T = 0.0


   call get_command_argument(1,args)
   read(args, *) N
   read(args, *) mass
   read(args, *) T
   read(args, *) cx
   read(args, *) cy
   read(args, *) cz

   !actual physics start

   V(1) = cx
   V(2) = cy
   V(3) = cz

   p(:)%mass = mass

   do i = 1,N
      call particleDist(V, T, p(i))
   end do

   open(newunit=io, file="vel.txt")
   do i = 1, N
      print *, p(i)%c
   end do
   close(io)
end program

module maxwellDist
   implicit none

   type :: particle
      real :: mass
      real :: x(3)
      real :: c(3)
   end type particle
   contains
      subroutine particleDist(V,T,p)
         type(particle), intent(inout) :: p
         real, intent(in) :: T
         real, intent(in) :: V(3)

         integer :: i
         real :: k, pi, beta, R(6)

         k = 1.380649E-23
         pi = 3.14159265359

         call random_number(R)

         beta = (p%mass/2/k/T)

         p%c(1) = V(1) + SIN(2*pi*R(1))*(-log(R(2)))**0.5
         p%c(2) = V(2) + SIN(2*pi*R(4))*(-log(R(3)))**0.5
         p%c(3) = V(3) + SIN(2*pi*R(5))*(-log(R(6)))**0.5
   end subroutine particleDist
end module maxwellDist


