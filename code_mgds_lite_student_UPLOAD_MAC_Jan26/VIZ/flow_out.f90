!=======================================================================
!  Subroutine : output_flow_field
!-----------------------------------------------------------------------
!  Purpose :
!    Calculates and outputs flow field quantities of interest
!=======================================================================
  subroutine output_flow_field()
  
     use flow_field_output_data
     use mgds_geometry
   
     Implicit None
   
     ! Storage for computed properties (local)
     real*8,Dimension(ncL3,NumFlowVar) :: flow_prop
     character*30,dimension(NumFlowVar) :: out_name

     integer :: iout
     character*30 :: var_string

     ! Set default output names here (Do not touch!!)
     do iout = 1,NumFlowVar
        write(var_string,'(I0)')iout
        var_string = 'Var'//trim(var_string)
        out_name(iout) = var_string
     enddo
    
     ! Calculation of flow field properties
     call flow_calc(flow_prop,out_name)
     write(*,*)"Finished calculating species and total flow field quantities"
     
     ! Output will be written in 2D if simulation is 2D (this is read in
     ! from restart)
     call flow_plot(flow_prop,out_name)
     write(*,*)"finished plotting flow field quantities"  
     
     ! Other output options can be added (i.e. line extraction)

  end subroutine output_flow_field

!=======================================================================
!  Subroutine : flow_calc
!-----------------------------------------------------------------------
!  Purpose :
!    Calculates  flow field quantities of interest (in flow_prop)
!  
!  Arguments : 
!    flow_prop [inout] array to store computed variables for output
!    out_name  [inout] storage for variable names
!=======================================================================
  subroutine flow_calc(flow_prop,out_name)

     use flow_field_output_data
     use mgds_inputvars
     use mgds_geometry
     use mgds_AMR

     Implicit None
    
     ! Storage per L3 cell (we store outputs here)
     real*8,dimension(ncL3,NumFlowVar),intent(inout) :: flow_prop
     ! Storage for output names
     character*30,dimension(NumFlowVar),intent(inout) :: out_name
   
     ! Geometric variables for each cell
     real*8 :: L3geom(1:6)
     
     integer :: icell,iNumSample
     real*8 :: npart
     real*8  :: Wp_ref,dvol_inv
     real*8 :: sumCx, sumCy, sumCz
     real*8 :: sumCx2, sumCy2, sumCz2
     real*8 :: bulkCx, bulkCy, bulkCz

     ! NOTE ! If you want to perform an AMR, you MUST fill in this array
     ! NOTE ! for each L3 cell. The first element should be the length
     ! NOTE ! you want to resolve to and the second is a local time 
     ! NOTE ! scale. The time scale doesn't do anything, but is useful
     ! NOTE ! for setting the timestep for the run after an AMR
     allocate(mfp_store(ncL3,2))
     ! This array is initialized to something crazy so that an AMR won't 
     ! be performed on accident if this array isn't filled correctly
     mfp_store(:,:) = -99.99d99

     ! If you want to change the output names from the default (Var1,
     ! Var2, ...) do so here.
     ! NOTE ! Geometric variables (x, y, and z) are handled
     ! NOTE ! automatically 

     out_name(1) = 'NParticles'
     out_name(2) = 'Vx'
     out_name(3) = 'Vy'
     out_name(4) = 'Vz'
     out_name(5) = 'compT'
     out_name(6) = 'compRho'
     out_name(7) = 'compP'

     ! get the number of samples
     iNumSample = runparams%nsam  
     ! Get reference particle weight (from parameters)
     Wp_ref =  runparams%WpRef
     do icell = 1,ncL3
   
        L3geom(1) = celld(icell) %x0
        L3geom(2) = celld(icell) %y0
        L3geom(3) = celld(icell) %z0
        L3geom(4) = celld(icell) %xe-celld(icell) %x0
        L3geom(5) = celld(icell) %ye-celld(icell) %y0
        L3geom(6) = celld(icell) %ze-celld(icell) %z0


        ! Get inverse volume (is correct even for cut-cells)
        dvol_inv = celld(icell) %dVolInv
   
        ! We only compute the average # of simulated particles per cell (NOT the
        ! number density) for a single species. If more output is desired, then
        ! the user must increase the size of NumFlowVar and add the required 
        ! computations here. 
        ! See the collide.f90 routine for more details on what sums(:,:) is 
        ! sampling.
        npart = celld(icell)%sums(1,1)/dble(iNumSample)
        if(iNumSample .gt. 0)then
           flow_prop(icell,1) = npart
        else
           flow_prop(icell,1) = 0.0d0
        endif

        sumCx = celld(icell)%sums(2,1)/dble(iNumSample)
        sumCy = celld(icell)%sums(3,1)/dble(iNumSample)
        sumCz = celld(icell)%sums(4,1)/dble(iNumSample)
        sumCx2 = celld(icell)%sums(5,1)/dble(iNumSample)
        sumCy2 = celld(icell)%sums(6,1)/dble(iNumSample)
        sumCz2 = celld(icell)%sums(7,1)/dble(iNumSample)
        print *, sumcx, sumcy, sumcz, sumcx2, sumcy2, sumcz2

        flow_prop(icell, 2) = sumCx/npart
        flow_prop(icell, 3) = sumCy/npart
        flow_prop(icell, 4) = sumCz/npart

        if(npart .gt. 0)then
           bulkCx = sumCx/npart
           bulkCy = sumCy/npart
           bulkCz = sumCz/npart
           flow_prop(icell, 5) = 1/3/1.380649d-23*28.02/1000/6.022e23* &
              ((sumCx2+sumCy2+sumCz2)/npart - bulkCx*bulkCx - bulkCy*bulkCy - bulkCz*bulkCz)
           print *, 1/3/1.380649d-23*28.02/1000/6.022e23*&
              ((sumCx2+sumCy2+sumCz2)/npart - bulkCx*bulkCx - bulkCy*bulkCy - bulkCz*bulkCz)
        else
           flow_prop(icell,5) = 0.0d0
        endif

        ! This is also where we would set our AMR factors (if we wanted)
        ! mfp_store(icell,1) = local_length_scale
        ! mfp_store(icell,2) = local_time_scale
     enddo
  
  end subroutine flow_calc

  !subroutine thermo_calc(bulkCx2, bulkCy2, bulkCz2, cellCx2, cellCy2, cellCz2, cellT, numParts)
  !   real*8, intent(in) :: bulkCx2, bulkCy2, bulkCz2, cellCx2, cellCy2, cellCz2, numParts
  !   real*8, intent(out) :: cellT, cellP, cellRho

  !   cellT = 1/3*1/1.380649d-23*28.02/1000/6.022*10d23*numParts(-bulkCx2 - bulkCy2 - bulkCz2 + cellCx2 + cellCy2 + cellCz2)
  !end subroutine thermo_calc



!=======================================================================
!  Subroutine : flow_plot
!-----------------------------------------------------------------------
!  Purpose :
!    Writes out computed data in tecplot format
!
!  Arguments :
!    flow_prop [in] array with stored output variables
!    out_name  [in] array with variable names
!  
!  Remarks :
!    Writes 2D or 3D data depending on flag_2D (read in from restart)
!=======================================================================
  subroutine flow_plot(flow_prop,out_name)

     use mgds_inputvars
     use flow_field_output_data
     use mgds_geometry

     Implicit None

     real*8,dimension(ncL3,NumFlowVar),intent(in) :: flow_prop
     character*30,dimension(NumFlowVar),intent(in) :: out_name

     integer :: npl,iNode,iout,n,nn
     real*8 :: xn, yn, zn
     character*20  :: int_string
     character*120 :: variable_list

     ! Temporary vector for writing in BLOCK tecplot format
     real*8,dimension(:),allocatable :: temp

     ! Retained from when a group of num_spec files were written
     integer :: file_id
     character*40 :: filename

     ! Set filename
     filename = "flow_field.dat"
     file_id = 101
     ! Set number of data points per line in block format
     npl = 5
     ! Check to make sure we don't have less than npl outputs
     if(npl .gt. ncL3) npl = ncL3
     
     if(flag_2D)then

       ! 4 nodes per element for FEQUAD
       nn = 4*ncL3

       ! Build variable list off of out_name
       variable_list = 'variables = x y'
       do iout = 1,NumFlowVar
          variable_list = trim(variable_list)//' '//trim(out_name(iout))
       enddo

       ! Write tecplot header
       open(unit=file_id,file=filename,status='unknown')
       write(file_id,*)trim(variable_list)
       if(NumFlowVar .eq. 1)then
           write(file_id,*)'zone T="flow2d", n=',nn,', e=', &
               ncL3,', DATAPACKING=BLOCK, ZONETYPE=FEQUADRILATERAL',&
               ', VARLOCATION=([3]=CELLCENTERED)'
       else
          write(int_string,'(I0)')NumFlowVar+2
          write(file_id,*)'zone T="flow2d", n=',nn,', e=', &
                  ncL3,', DATAPACKING=BLOCK, ZONETYPE=FEQUADRILATERAL',&
                  ', VARLOCATION=([3-',trim(int_string),']=CELLCENTERED)'
       endif
     else

        ! 8 nodes per element for FEBRICK
        nn = 8*ncL3

        ! Build variable list off of out_name
        variable_list = 'variables = x y z'
        do iout = 1,NumFlowVar
           variable_list = trim(variable_list)//' '//trim(out_name(iout))
        enddo

        ! Write tecplot header
        open(unit=file_id,file=filename,status='unknown')
        write(file_id,*)trim(variable_list)
        if(NumFlowVar .eq. 1)then
           write(file_id,*)'zone T="flow3d", n=',nn,', e=', &
                ncL3,', DATAPACKING=BLOCK, ZONETYPE=FEBRICK',&
                ', VARLOCATION=([4]=CELLCENTERED)'

        else
           write(int_string,'(I0)')NumFlowVar+3
           write(file_id,*)'zone T="flow3d", n=',nn,', e=', &
                   ncL3,', DATAPACKING=BLOCK, ZONETYPE=FEBRICK',&
                   ', VARLOCATION=([4-',trim(int_string),']=CELLCENTERED)'
        endif
     endif


     ! Allocate temporary array to store data
     allocate(temp(nn))

     ! Store x location data
     if(flag_2D)then
        do n=1,ncL3
           temp(4*(n-1)+1) = celld(n)%x0
           temp(4*(n-1)+2) = celld(n)%xe
           temp(4*(n-1)+3) = celld(n)%xe
           temp(4*(n-1)+4) = celld(n)%x0
        enddo
     else
        do n=1,ncL3
           temp(8*(n-1)+1) = celld(n)%x0
           temp(8*(n-1)+2) = celld(n)%xe
           temp(8*(n-1)+3) = celld(n)%xe
           temp(8*(n-1)+4) = celld(n)%x0
           temp(8*(n-1)+5) = celld(n)%x0
           temp(8*(n-1)+6) = celld(n)%xe
           temp(8*(n-1)+7) = celld(n)%xe
           temp(8*(n-1)+8) = celld(n)%x0
        enddo
     endif

     ! Write x location data (nodal)
     do n = 1,nn/npl
        write(file_id,100)temp(npl*(n-1)+1:npl*(n-1)+npl)
     enddo
     if(mod(nn,npl) .ne. 0)then
        write(file_id,100)temp(nn-mod(nn,npl)+1:nn)
     endif

     ! Store y location data
     if(flag_2D)then
        do n=1,ncL3
           temp(4*(n-1)+1) = celld(n)%y0
           temp(4*(n-1)+2) = celld(n)%y0
           temp(4*(n-1)+3) = celld(n)%ye
           temp(4*(n-1)+4) = celld(n)%ye
        enddo
     else
        do n=1,ncL3
           temp(8*(n-1)+1) = celld(n)%y0
           temp(8*(n-1)+2) = celld(n)%y0
           temp(8*(n-1)+3) = celld(n)%ye
           temp(8*(n-1)+4) = celld(n)%ye
           temp(8*(n-1)+5) = celld(n)%y0
           temp(8*(n-1)+6) = celld(n)%y0
           temp(8*(n-1)+7) = celld(n)%ye
           temp(8*(n-1)+8) = celld(n)%ye
        enddo
     endif

     ! Write y location data (nodal)
     do n = 1,nn/npl
        write(file_id,100)temp(npl*(n-1)+1:npl*(n-1)+npl)
     enddo
     if(mod(nn,npl) .ne. 0)then
        write(file_id,100)temp(nn-mod(nn,npl)+1:nn)
     endif

     ! Only write z location for 3D
     if(flag_2D .eqv. .false.)then
        ! Store z location data
        do n=1,ncL3
           temp(8*(n-1)+1) = celld(n)%z0
           temp(8*(n-1)+2) = celld(n)%z0
           temp(8*(n-1)+3) = celld(n)%z0
           temp(8*(n-1)+4) = celld(n)%z0
           temp(8*(n-1)+5) = celld(n)%ze
           temp(8*(n-1)+6) = celld(n)%ze
           temp(8*(n-1)+7) = celld(n)%ze
           temp(8*(n-1)+8) = celld(n)%ze
        enddo

        ! Write y location data (nodal)
        do n = 1,nn/npl
           write(file_id,100)temp(npl*(n-1)+1:npl*(n-1)+npl)
        enddo
        if(mod(nn,npl) .ne. 0)then
           write(file_id,100)temp(nn-mod(nn,npl)+1:nn)
        endif
     endif

     ! We no longer need temp
     deallocate(temp)

     ! Loop over output variables 
     do iout = 1,NumFlowVar
        ! Write data for iout
        do n = 1,ncL3/npl
           write(file_id,100)flow_prop(npl*(n-1)+1:npl*(n-1)+npl,iout)
        enddo
        if(mod(ncL3,npl) .ne. 0)then
           write(file_id,100)flow_prop(ncL3-mod(ncL3,npl)+1:ncL3,iout)
        endif
     enddo

     ! Lastly, write connectivity 
     if(flag_2D)then
        ! Write nodal connectivity
        do n = 1,ncL3
          iNode = (n - 1) * 4
          write(file_id,110)iNode+1,iNode+2,iNode+3,iNode+4
        end do
     else
        do n = 1, ncL3
          iNode = (n - 1) * 8
          write(file_id,110)iNode+1,iNode+2,iNode+3,iNode+4,&
                            iNode+5,iNode+6,iNode+7,iNode+8
        end do
     endif

     ! Close file
     close(file_id)

100  format(100(1x,e13.6))
110  format(100(1x,i9))

  end subroutine flow_plot
