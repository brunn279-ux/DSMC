!=======================================================================
!  Subroutine : mgds_USER_MeshCube
!-----------------------------------------------------------------------
!  Purpose : Subroutine to override the default L1,L2,L3 spacings
! 
!  Remarks : 
!    -DO NOT USE THIS IF YOU DON'T KNOW WHAT YOU ARE DOING
!=======================================================================
  subroutine mgds_USER_MeshCube()
    Implicit None
  
  end subroutine mgds_USER_MeshCube

!=======================================================================
!  Subroutine : mgds_USER_Init
!-----------------------------------------------------------------------
!  Purpose : Subroutine to initialize user functions and data structures
! 
!  Remarks : 
!    -You may have to add specific modules to access certain built-in 
!     variables or data structures
!=======================================================================
  subroutine mgds_USER_Init
  
     use userdata, only : mycell_id
     ! We need to grab the # of cells and the celld array (array
     ! of mgds cell types -> see main_struct definitions)
     use mgds_geometry, only : num_cd,celld
  
     Implicit None

     integer :: icell
     real*8 :: x_point,y_point
     real*8 :: xo,xe,yo,ye

     ! Example: How to find which cell a point is in and then flag it

     ! For Homework 1, you will need to output some data from a single 
     ! cell. This is is how would do this (this assumes 2D!)

     ! Point of interest
     x_point = lam*50
     y_point = lam*25
     
     ! Initial value to check if we didn't flag a cell
     mycell_id = -10
     ! Loop over the cells and check if this point lies within 
     do icell = 1,num_cd

        ! Get origin and termination point of mgds cell in 2D
        xo = celld(icell)%x0
        yo = celld(icell)%y0
        xe = celld(icell)%xe
        ye = celld(icell)%ye

        ! Check if we are in this cell (only inclusive at lower bounds,
        ! otherwise multiple cells can trigger right a cell boundaries)
        if( (x_point .ge. xo) .and. &
            (x_point .lt. xe) .and. &
            (y_point .ge. yo) .and. &
            (y_point .lt. ye) ) then

           ! We have found our cell!
           mycell_id = icell
        endif
     enddo

     if(mycell_id .eq. -10)then
        write(*,*)'Point not found in any cell!'
     endif

  end subroutine mgds_USER_Init

!=======================================================================
!  Subroutine : mgds_USER_ScreenOutput
!-----------------------------------------------------------------------
!  Purpose : Subroutine to print to the screen selected vqariables
! 
!  Remarks : 
!    -You may have to add specific modules to access certain built-in 
!     variables or data structures
!    -You may not want to output at every timestep; that is what the
!     if(mod(nit,1) statement is for. Change 1 to how often you want to
!     output to the screen.
!=======================================================================
  subroutine mgds_USER_ScreenOutput()
  
     use userdata
     use mgds_dsmc
   
     Implicit None
   
     ! Write iteration, # of particles, # of collisions, and # of bad
     ! collisions to screen every timestep
     if(mod(nit,1) .eq. 0)write(*,*)nit,num_ptcl,numtotcoll,numbadcoll
  
  end subroutine mgds_USER_ScreenOutput

!=======================================================================
!  Subroutine : mgds_USER_PostWriteRestart
!-----------------------------------------------------------------------
!  Purpose : Subroutine to perform operations right after a restartfile 
!    was created
! 
!  Remarks : 
!    -This is a good place to write user data to file 
!=======================================================================
  subroutine mgds_USER_PostWriteRestart()
  
     use userdata

     Implicit None
     ! For now, does nothing
  
  end subroutine mgds_USER_PostWriteRestart

!=======================================================================
!  Subroutine : mgds_USER_PreTimestep
!-----------------------------------------------------------------------
!  Purpose : Subroutine to perform operations right before a timestep is
!    started
!=======================================================================
  subroutine mgds_USER_PreTimestep()   

     use userdata

     Implicit None
     ! For now, does nothing
  end subroutine mgds_USER_PreTimestep
!=======================================================================
!  Subroutine : mgds_USER_PostTimestep
!-----------------------------------------------------------------------
!  Purpose : Subroutine to perform operations right after a timestep is
!    completed
!=======================================================================
  subroutine mgds_USER_PostTimestep()
  
     use mgds_dsmc    ! you know timestep number from here
     use mgds_inputvars 
     use mgds_geometry
   
     Implicit None
   
     integer :: file_id,npart,num_pcell,ncell
     type(particle_t),Pointer  :: pCurrent
   
     ! EXAMPLE: Write particle location at every timestep to visualize 
     !          location.
     ! NOTE ! You uncomment this section to use it

     ! Loop over all cells and grab particles (ParIniArr is only 
     ! 1 million elements big, but you shouldn't be writing more than
     ! this to a file anyways...

     !num_pcell = 0
     !do ncell = 1,num_cd
     !   ! First, grab particles from cell list and store in ParIniArr
     !   if(associated( celld(ncell)%head )) then
     !      num_pcell = num_pcell + 1
     !      pCurrent => celld(ncell) %head
     !      parIniArr(num_pcell) = pCurrent    ! copies properties
   
     !      do while(associated(pCurrent%next))
     !         num_pcell = num_pcell + 1
     !         pCurrent => pCurrent%next
     !         parIniArr(num_pcell) = pCurrent  ! copies properties
     !      enddo
     !   endif
     !enddo
   
     !! Write to a tecplot file
   
     !file_id = 100
     !if(nit .eq. 1)then
     !   open(unit = file_id,file='part_loc.dat')
     !   write(file_id,*)'VARIABLES = "x" "y"'
     !endif
   
     !write(file_id,*)'Zone T="t=',REAL(dTref*real(nit-1), KIND=4),&
     !   '", STRANDID=',nit,', SOLUTIONTIME=',&
     !   REAL(dTref*real(nit-1), KIND=4),', I=',num_pcell,', F=POINT'
   
     !do npart = 1,num_pcell
     !   write(file_id,*)parIniArr(npart)%px,parIniArr(npart)%py
     !enddo
   
     !if(nit .eq. nstep)then
     !   close(file_id)
     !endif
  end subroutine mgds_USER_PostTimestep

