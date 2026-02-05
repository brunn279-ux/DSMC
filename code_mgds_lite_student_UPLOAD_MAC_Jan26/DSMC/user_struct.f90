!=======================================================================
!  Module : userdata
!-----------------------------------------------------------------------
!  Purpose : Module to store the variables associated with user routines
!
!  Remarks :
!    -To access variables in a subroutine, add 'use userdata' after
!     the subroutine name but before 'Implicit None'
!=======================================================================
  module userdata

     Implicit None

     ! Modules in fortran are very useful. For example, this is a good
     ! place to store parameters e.g. pi (parameter just means you 
     ! can't change the value later)
     real*8,parameter :: pi = 4.0d0*atan(1.0d0)
     real*8, parameter :: k_b = 1.380649d-23 !m^2 kg s^-2 K^-1
     real*8, parameter :: M_ar = 39.948d0 ! g/mol
     real*8, parameter :: M_N2 = 28.014d0 ! g/mol
     real*8, parameter :: bulkT = 200d0 !K
     real*8, parameter :: bulkrho = 1.6d-5 !kg/m^3
     real*8, parameter :: d_N2 = 4.0d-10
     real*8, parameter :: bulkV(3) = [50.0d0, 0.0d0, 0.0d0]
     real*8, parameter :: bulkMag = (bulkV(1)*bulkV(1)+bulkV(2)*bulkV(2)+bulkV(3)*bulkV(3))**0.5
     real*8, parameter :: lam = M_N2/1000/6.022e23/bulkrho/2**0.5/pi/d_N2**2
     real*8, parameter :: tau_c = lam/bulkMag
     ! Cell where our point of interest lies in
     integer :: mycell_id
  end module userdata
