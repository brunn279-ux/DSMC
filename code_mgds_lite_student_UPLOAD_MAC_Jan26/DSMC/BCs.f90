!=======================================================================
!  Subroutine : user_Volume_BoundaryConditions
!-----------------------------------------------------------------------
!  Purpose :
!    Subroutine to put particles in inflow cells using a volume
!    particle generation method
!
!  Arguments :
!    itimestep  [in]  current timestep
!=======================================================================
  subroutine mgds_USER_Volume_BoundaryConditions(itimestep)

    use mgds_dsmc      ! ParIniArr from here
    use mgds_geometry  ! celld structure from here
    use userdata 
    use mgds_random, only : mgds_randu
  
    Implicit None

    integer,intent(IN) :: itimestep

    integer :: icell,nptcl_add,iPart,ngen
    real*8  :: rRand,rand_angle
    real*8  :: rx,ry,rz,vx,vy,vz
    real*8  :: erot,evib 
  
    real*8 :: xo
  
    do icell=1,num_cd
  
       ! Generate particles based on x-position
       xo = celld(icell)%x0
  
       ! Decides which cells generate particles
       ! NOTE ! This can be anything! The user is responsible for
       ! NOTE ! deciding how to do this.
       if(xo <= lam*5 .or. xo >= lam*95)then

          ! Before generating particles, we remove the current
          ! particles from the linked list (DO NOT REMOVE)
          call mgds_Particles_RemoveAllFromCell(icell)
  
          ! Choose how many particles to generate in this cell
          ! (hard-coded for test)
          ngen = 20 

          ! Next, loop over particles and generate properties
          ! Particle properties are temporarily stored in ParIniArr
          ! before being added to the linked list
          nptcl_add = 0        ! number of particles in local cell
          do iPart = 1,ngen
             nptcl_add = nptcl_add + 1  ! local cell particle count
          
             ! Need to give the particle a position in the cell and 
             ! a velocity.

             ! Position: Random in x, center of cell in y and z
             ! Random number between 0 and 1
             rRand = mgds_randu()
             rx = celld(icell)%x0+&
                ( celld(icell)%xe - celld(icell)%x0 )*rRand
             ry = celld(icell)%y0+&
                ( celld(icell)%ye - celld(icell)%y0 )/2.0d0
             rz = celld(icell)%z0+&
                ( celld(icell)%ze - celld(icell)%z0 )/2.0d0

             ! Velocity: no z velocity (2D), but randomly 
             ! pick a velocity angle between -45 and 45. 
             ! Velocity magnitude is 10 m/s
             !vmag = 10.0d0
             !V = [50.0d0, 0.0d0, 0.0d0]
             ! Random number between 0 and 1
             rRand = mgds_randu()
             ! Random angle between -45 and 45
             !rand_angle = -45.0d0+90.0d0*rRand
             ! Convert to radians (pi from userdata!)
             !rand_angle = rand_angle*pi/180.0d0
             !vx = vmag*cos(rand_angle)
             !vy = vmag*sin(rand_angle)
             !vz = 0.0d0
             call velocityGeneration(bulkV(1), vx)
             call velocityGeneration(bulkV(2), vy)
             call velocityGeneration(bulkV(3), vz)

             ! We stage the particle info in parIniArr, which is 
             ! an array of user-defined particle types. parIniArr is 
             ! large (1 million), so it is very unlikely you will 
             ! exceed the max size. The particles stored in 
             ! parIniArr are added to the linked list of the current
             ! cell in mgds_Particles_AddCelldataParticles
             parIniArr(nptcl_add)%Px = rx
             parIniArr(nptcl_add)%Py = ry
             parIniArr(nptcl_add)%Pz = rz
             parIniArr(nptcl_add)%Vx = vx
             parIniArr(nptcl_add)%Vy = vy
             parIniArr(nptcl_add)%Vz = vz
             ! No internal energy
             parIniArr(nptcl_add)%Erot = 0.0d0
             parIniArr(nptcl_add)%Evib = 0.0d0

             ! This is the particle species index. For a single species,
             ! this is not important, but will be for nsp > 1. Note that 
             ! the # of species must be consistent with the inputdeck 
             ! values (so mgds_lite can handle io for you)
             parIniArr(nptcl_add)%spectype = 1
  
             ! These are user variables that can be used to store
             ! info per particle:
             ! 1) rfoo -> a real*8 
             ! 2) ifoo -> an integer
             ! For now, they are now used by the underlying code 
             ! in any way.
             parIniArr(nptcl_add)%rfoo = 1.0d0
             parIniArr(nptcl_add)%ifoo = 1
          enddo
  
          ! Add nptcl_add particles to the linked list in cell
          ! icell
          call mgds_Particles_AddCelldataParticles(icell,nptcl_add)
       endif
    enddo
  
  end Subroutine mgds_USER_Volume_BoundaryConditions

!=======================================================================
!  Subroutine : mgds_USER_surface
!-----------------------------------------------------------------------
!  Purpose :
!    This routines allows the user to alter particle properties after
!    a collision with a surface is detected. 
!  
!  Remarks
!    -Note that why the particle position is available, it is 'read-
!     only' -> that is, you cannot change the particle position from 
!     this routine.
!
!  Arguments :
!    itri        [in]   index of triangle that particle collides with
!    ispec      [inout] species index of particle
!    vx,vy,vz   [inout] particle velocity in x, y, and z
!    erot,evib  [inout] internal energy for rotation and vibration
!    ifoo,rfoo  [inout] user integer and real variables 
!    delete     [inout] flag to delete particle (if desired)
!=======================================================================
  subroutine mgds_USER_surface(itri,ispec,vx,vy,vz,erot,evib,&
       rfoo,ifoo,px,py,pz,delete)
  
     use userdata
     use mgds_dsmc
     use mgds_surfaces  
     use mgds_random, only : mgds_randu
   
     Implicit None

     integer,intent(in) :: itri
     integer,intent(inout) :: ispec
     real*8, intent(inout) :: vx,vy,vz
     real*8, intent(inout) :: erot,evib
     integer,intent(inout) :: ifoo
     real*8, intent(inout) :: rfoo
     real*8, intent(in)    :: px,py,pz
     logical,intent(inout) :: delete

     integer :: izone
     real*8 :: nx,ny,nz

     real*8 :: vn_mag
     real*8,dimension(3) :: vn,vt
   
     ! Get some info for this triangle (triang is user-defined triangle
     ! type -> see data structures)
     ! Surface normal
     nx = triang(itri)%nx
     ny = triang(itri)%ny
     nz = triang(itri)%nz
     ! Triangle zone (if multiple surfaces are loaded, this is the index
     ! of the surface this triangle is apart of. Useful if different 
     ! behavior is desired for different surfaces.
     izone = triang(itri)%izone
   
     ! Flag to determine whether particle is deleted or not.
     ! For now, the default behavior is to NOT delete particles.
     delete = .false.
   
     ! As an example, if we wanted to delete all particles from the 2nd
     ! zone, we would do:
     ! if(izone .eq. 2) then
     !    delete = .true.
     ! endif
   
     ! No need to set particle velocities if we are deleting particle
     ! anyway.
     if(delete)then
        return
     endif
     
     ! By default, specular reflections are done on surfaces. 
     ! For specular reflections, the normal component is simply reflected.
     vn_mag = nx*vx+vy*ny+vz*nz
     vn(1) = vn_mag*nx
     vn(2) = vn_mag*ny
     vn(3) = vn_mag*nz

     vt(1) = vx-vn(1)
     vt(2) = vy-vn(2)
     vt(3) = vz-vn(3)

     ! New velocity just has the normal component direction
     ! reversed
     vx = vt(1)-vn(1)
     vy = vt(2)-vn(2)
     vz = vt(3)-vn(3)

     ! No internal energy
     erot = 0.0d0
     evib = 0.0d0

  end subroutine mgds_USER_surface

   subroutine velocityGeneration(V, c)
      use userdata
      use mgds_random, only : mgds_randu

      implicit none

      real*8, intent(in) :: V
      real*8, intent(out) :: c
      real*8 :: R1, R2, indMass, beta
      R1 = mgds_randu()
      R2 = mgds_randu()
      indMass = M_ar/1000/6.022e23

      beta = (indMass/2/k_b/bulkT)**0.5

      c = V + 1/beta*SIN(2*pi*R1)*(-log(R2))**0.5
   end subroutine
