!=======================================================================
!  Subroutine : mgds_USER_compute_grid_scale
!-----------------------------------------------------------------------
!  Purpose : The user must set the grid scale here. 
! 
!  Arguments : 
!    grid_scale  [out] reference scale grid that initial grid will 
!                      be sized to. Note that scale_factors from the 
!                      input deck are applied outside this routine.
!
!  Remarks :
!    -This scale applies to all Cartesian directions for 3D and to the 
!     x and y direction for 2D (this is handled elsewhere)
!    -We need an integer # of cells, so the grid will not necessarily 
!     be sized exactly to your scale. The grid scale you provide is the
!     maximum allowable size. Here is an example:
!       For x direction, domain size is: 3.0 m. Your grid scale is 1.1.
!       MGDS will find the minimum # of cells that satisfies your grid
!       scales, # of cells in x = 4 and the x grid scale = 0.75. 
!    -This integer sizing happens to all three grid levels (L1,L2,L3)
!     starting from L1 and the scaling factor provided in the input
!     deck in L1_scale_factor (so this factor X grid_scale). Therefore,
!     you should not expect a specific grid sizing, just expect that 
!     at maximum, your grid will be smaller than grid_scale X scale_
!     factor
!    -The factor from the inputdeck is applied outside this routine, so
!     you can just modify L3_scale_factor instead of recompiling with 
!     different grid_scale.
!=======================================================================
  subroutine mgds_USER_compute_grid_scale(grid_scale)

     use userdata

     Implicit None

     real*8,intent(out) :: grid_scale

     grid_scale = lam

  end subroutine mgds_USER_compute_grid_scale

!=======================================================================
!  Subroutine : mgds_USER_compute_dT_Wp
!-----------------------------------------------------------------------
!  Purpose : The user must set the timestep and particle weight (i.e. 
!     how many real particles each simulated particle represents)
! 
!  Arguments : 
!     No arguments 
!
!  Remarks :
!    -The outputs of this routine are stored in the mgds_inputvars and 
!     therefore there are no inputs or outputs.
!    -There are redundant variables for the timestep that are set 
!     outside this routine.
!       dTRef -> accessed through mgds_dsmc module
!    -NOTE: There are variables in the inputdeck that you may find 
!           useful here, but you DO NOT have to use them. 
!           The variables are:
!               InputDeck%dT_fact
!               InputDeck%NumPartCell
!=======================================================================
subroutine mgds_USER_compute_dT_Wp
  use mgds_inputvars, only : runparams,InputDeck
  use userdata
  implicit none
 
  ! InputDeck%dT_fact and InputDeck%NumPartCell are not used by default!! 
  runparams %WpRef = 1.0d0
  runparams %dTRef = tau_c/5

end subroutine mgds_USER_compute_dT_Wp
