!=======================================================================
!  Program : mgdsviz_lite
!-----------------------------------------------------------------------
!  Purpose :
!    Main driver for the DSMC course code post-processor.
!
!    1) First read in the species information from the species data
!       file. 
!    2) Then, read in restart file info: inpuctdeck, DSMC run parameters,
!       geometry, surfaces and sampling info (surfaces and grid)
!    3) Output routines for flowfield 
!    4) Output routines for surfaces 
!    5) An Adaptive Mesh Routine (AMR) can also be called, which 
!       resolves the grid to a user specified local mean-free path 
!       (+ a few other options)
!=======================================================================

  program mgdsviz_lite
  
     use flow_field_output_data
     use mgds_dsmc, only : nsp
  
     Implicit None

     logical :: AMR_switch

     ! User needs to set set some things here (we can make another 
     ! input deck for this if we want, but I think this will be fine)

     ! Whether or not to perform an AMR
     AMR_switch = .false.

     ! The number of flow variables to output
     NumFlowVar = 7
     ! The number of surface variables to output
     NumSurfVar = 1

     ! All relevant info (inputs, geometry, surfaces, and sampling)
     ! are read here
     call mgdsviz_Read()

     ! Set the "total" number of species (+1 on number of species for
     ! averages)
     totspec = nsp + 1
  
     ! flow field output   
     call output_flow_field()
     write(*,*)"finished plotting flow field"
         
     ! Now, need to read in the surface sampling data and plot
     call output_surface()
     write(*,*)"finished plotting surface sampling"
  
     ! AMR routines 
     if(AMR_switch)then
        write(*,*)'[================== ENTERING AMR ===================]'
        call SYSTEM('sleep 1')
        ! Perform AMR
        call mgdsviz_AMR
        ! Update AMR grid in restart
        call mgdsviz_UpdateGeometry
     endif
  
     write(*,*)'Post-processor completed.'
  
  end program mgdsviz_lite

!=======================================================================
!  Subroutine : mgds_TerminateError
!-----------------------------------------------------------------------
!  Purpose :
!    Routine to throw an error message and end the program execution.
!=======================================================================

  subroutine mgds_TerminateError()

     Implicit None

     Character*(50) :: date_string

     write(*,*)'The program terminated with an error'
  
     call fdate( date_string )
     write(*,*)'Terminating on: ',trim(date_string)
     STOP

  end subroutine

