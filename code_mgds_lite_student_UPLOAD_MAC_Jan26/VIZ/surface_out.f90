!=======================================================================
!  Subroutine : output_surface
!-----------------------------------------------------------------------
!  Purpose :
!    Subroutine calls routines to compute desired output on surfaces
!    and then write data to tecplot format
!=======================================================================
  subroutine output_surface

     use flow_field_output_data
     use mgds_surfaces
  
     Implicit none
  
     ! Storage for computed surface properties (local storage)
     real*8,dimension(ntri,NumSurfVar) :: surf_prop
     character*30,dimension(NumFlowVar) :: out_name

     integer :: iout
     character*30 :: var_string

     ! Set default output names here (Do not touch!!)
     do iout = 1,NumFlowVar
        write(var_string,'(I0)')iout
        var_string = 'Var'//trim(var_string)
        out_name(iout) = var_string
     enddo
  
     call surf_calc(surf_prop,out_name)
     write(*,*)"finished calculating species and total surface quantities"
  
     call surf_plot(surf_prop,out_name)
     write(*,*)"finished plotting species and total surface quantities"
  
  end subroutine output_surface

!=======================================================================
!  Subroutine : surf_quantities_calc
!-----------------------------------------------------------------------
!  Purpose :
!    Subroutine calculates and stores per surface quantities 
! 
!  Arguments : 
!    surf_prop  [inout] storage for computed surface properties
!    out_name   [inout] storage for variable names
!=======================================================================
  subroutine surf_calc(surf_prop,out_name)

     use flow_field_output_data
     use mgds_inputvars
     use mgds_surfaces
  
     Implicit None
  
     ! Storage per triangle (we store outputs here)
     real*8,dimension(ntri,NumSurfVar),intent(inout) :: surf_prop
     ! Storage for output names
     character*30,dimension(NumFlowVar),intent(inout) :: out_name
  
     integer :: iNumSample,itri,iout
     real*8 :: dt,area,Wp_ref
     
     ! If you want to change the output names from the default (Var1,
     ! Var2, ...) do so here.
     ! NOTE ! Geometric variables (x, y, and z) are handled
     ! NOTE ! automatically
     out_name(1) = 'nflux'

     ! Get the number of samples
     iNumSample = runparams%nsam
     Wp_ref =  runparams%WpRef
     ! Set the timestep size
     dt = runparams%dtref
  
     ! Initialize the triangle information before any surface info calculated
     do iout = 1,NumSurfVar
        do itri = 1,ntri
           surf_prop(itri,iout) = 0.0d0
        enddo
     enddo


     ! Loop through all triangles
     do itri = 1,ntri
        ! Set the triangle area (ss2 is 2*Surface area)
        area = triang(itri)%ss2 /2.0d0
  
        if(dt*area*real(iNumSample) .gt. 0.0d0)then
           surf_prop(itri,1) =  sdata(itri)%sums(1,1)*Wp_ref/iNumSample/dt/area
           surf_prop(itri,1) = area
        else
           surf_prop(itri,1) = 0.0d0
        endif
     enddo
  
  end subroutine surf_calc

!=======================================================================
!  Subroutine : plot_surface_quantities
!-----------------------------------------------------------------------
!  Purpose :
!    Plot per species and total surface quantities
! 
!  Arguments : 
!    surf_prop  [in] array with computed surface properties
!    out_name   [in] array with variable names
!=======================================================================
  subroutine surf_plot(surf_prop,out_name)

     use flow_field_output_data
     use mgds_surfaces
  
     Implicit None
  
     real*8,dimension(ntri,NumSurfVar),intent(in) :: surf_prop
     character*30,dimension(NumFlowVar),intent(in) :: out_name

     real*8 :: xn, yn, zn
     integer :: file_id,itri,iout,iNode,n,nn
     character*20 :: int_string
     character*40 :: filename
     character*120 :: variable_list

     ! Temporary vector for writing in BLOCK tecplot format
     real*8,dimension(:),allocatable :: temp
     integer :: npl
     
     ! Setup the file output
     file_id = 100
     filename = "surface.dat"
     
     ! Set number of writes per line
     npl = 5

     ! Build variable list off of out_name
     variable_list = 'variables = x y z'
     do iout = 1,NumFlowVar
        variable_list = trim(variable_list)//' '//trim(out_name(iout))
     enddo

     ! Write header
     open(unit=file_id,file=filename,status='unknown')
     write(file_id,*)trim(variable_list)
     if(NumSurfVar .eq. 1)then
        write(file_id,*)'zone T="surf", n=',ntri*3,', e=', &
           ntri,', et=triangle, f=feblock',&
           ', VARLOCATION=([4]=CELLCENTERED)'
     else
        write(int_string,'(I0)')NumSurfVar+3
        write(file_id,*)'zone T="surf", n=',ntri*3,', e=', &
           ntri,', et=triangle, f=feblock',&
           ', VARLOCATION=([4',trim(int_string),']=CELLCENTERED)'
     endif
     
     nn = 3*ntri
     allocate(temp(nn))
     ! Store x location
     do itri=1,ntri
        temp(3*(itri-1)+1) = triang(itri)%x1
        temp(3*(itri-1)+2) = triang(itri)%x2
        temp(3*(itri-1)+3) = triang(itri)%x3
     enddo

     ! Write x location data (nodal)
     do n = 1,nn/npl
        write(file_id,100)temp(npl*(n-1)+1:npl*(n-1)+npl)
     enddo
     if(mod(nn,npl) .ne. 0)then
        write(file_id,100)temp(nn-mod(nn,npl)+1:nn)
     endif

     ! Store y location
     do itri=1,ntri
        temp(3*(itri-1)+1) = triang(itri)%y1
        temp(3*(itri-1)+2) = triang(itri)%y2
        temp(3*(itri-1)+3) = triang(itri)%y3
     enddo

     ! Write y location data (nodal)
     do n = 1,nn/npl
        write(file_id,100)temp(npl*(n-1)+1:npl*(n-1)+npl)
     enddo
     if(mod(nn,npl) .ne. 0)then
        write(file_id,100)temp(nn-mod(nn,npl)+1:nn)
     endif

     ! Store z location
     do itri=1,ntri
        temp(3*(itri-1)+1) = triang(itri)%z1
        temp(3*(itri-1)+2) = triang(itri)%z2
        temp(3*(itri-1)+3) = triang(itri)%z3
     enddo

     ! Write z location data (nodal)
     do n = 1,nn/npl
        write(file_id,100)temp(npl*(n-1)+1:npl*(n-1)+npl)
     enddo
     if(mod(nn,npl) .ne. 0)then
        write(file_id,100)temp(nn-mod(nn,npl)+1:nn)
     endif

     ! We not longer need temp
     deallocate(temp)

     ! Loop over surface variables and write
     do iout = 1,NumSurfVar
        do n = 1,ntri/npl
           write(file_id,100)surf_prop(npl*(n-1)+1:npl*(n-1)+npl,iout)
        enddo
        if(mod(ntri,npl) .ne. 0)then
           write(file_id,100)surf_prop(ntri-mod(ntri,npl)+1:ntri,iout)
        endif
     enddo
     
     ! Write node index
     do itri = 1,ntri
        iNode = (itri-1)*3
        write(file_id,110)iNode+1,iNode+2,iNode+3
     enddo
  
     ! Close file
     close(file_id)
  
100  Format(100(1x,e16.9))
110  Format(100(1x,i9))
  
  end subroutine surf_plot
