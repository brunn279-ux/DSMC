program makeVelocities

   use maxwellDist

   !input of number of particles bulk velocity, etc
   integer :: N
   real(8) :: mass, T, cx, cy, cz
   character(len=48) :: args
   type(particle), allocatable :: p(:)
   real(8) :: V(3)
   integer :: i
   integer :: io


   N = 1
   mass = 1.0
   T = 0.0


   call get_command_argument(1,args)
   read(args, *) N
   call get_command_argument(2,args)
   read(args, *) mass
   call get_command_argument(3,args)
   read(args, *) T
   call get_command_argument(4,args)
   read(args, *) cx
   call get_command_argument(5,args)
   read(args, *) cy
   call get_command_argument(6,args)
   read(args, *) cz

   allocate(p(N))

   !actual physics start

   V(1) = cx
   V(2) = cy
   V(3) = cz

   !print *, "size(p) =", size(p)
   p(:)%mass = mass/6.02e23/1000

   do i = 1,N
      call particleDist(V, T, p(i))
   end do

   open(newunit=io, file="vel.txt")
   write(io, *)  "cx    cy    cz"
   write(io, *)  "--------------"
   do i = 1, N
      write(io, *) p(i)%c(1),",",p(i)%c(2),",",p(i)%c(3)
   end do
   close(io)

   print *, "N: ",N, "Mass: ", mass, "T: ", T,"cx: ", cx, "cy: ", cy, "cz: ", cz
   !print *, p(:)

   deallocate(p)
end program

module maxwellDist
   implicit none

   type :: particle
      real(8) :: mass
      real(8) :: x(3)
      real(8) :: c(3)
   end type particle
   contains
      subroutine particleDist(V,T,p)
         type(particle), intent(inout) :: p
         real(8), intent(in) :: T
         real(8), intent(in) :: V(3)

         integer :: i
         real(8) :: k, pi, beta, R(6)

         k = 1.380649E-23
         pi = 3.14159265359

         call random_number(R)

         !print *, R

         beta = (p%mass/2/k/T)**0.5
         !print *, beta, " ", p%mass, " ", k, " ", T

         p%c(1) = V(1) + 1/beta*SIN(2*pi*R(1))*(-log(R(2)))**0.5
         p%c(2) = V(2) + 1/beta*SIN(2*pi*R(4))*(-log(R(3)))**0.5
         p%c(3) = V(3) + 1/beta*SIN(2*pi*R(5))*(-log(R(6)))**0.5
   end subroutine particleDist
end module maxwellDist


