!=======================================================================
!  Subroutine : mgds_collide
!-----------------------------------------------------------------------
!  Purpose : Subroutine to perform collisions on the particles stored 
!    in ParIniArr (up to nLocalParticles; all other info is junk).
!
!  Arguments :
!    icell           [in]    current mgds cell index
!    nLocalParticles [inout] # of particles in cell (in ParIniArr)
!    nBadColl        [out]   # of bad collisions
!    nColl           [out]   # of collisions
!    iNumThirdBody   [inout] # of third bodies reserved
!    timestep        [in]    current timestep
! 
!  Remarks : 
!    -This doesn't change behavior as compared to the original MGDS 
!     code, but allows for a simpler alteration when the MGDS collision
!     routines are removed for the student version.
!=======================================================================
  subroutine mgds_collide(icell,           &
                          nLocalParticles, &
                          nPartSpec,       &
                          nBadColl,        &
                          ncoll,           &
                          iNumThirdBody,   &
                          timestep)


     use mgds_geometry
     use mgds_dsmc
     use mgds_random, only : mgds_randu
     use userdata, only : mycell_id

     Implicit none

     integer,intent(in)     :: icell,timestep
     integer,intent(inout)  :: nLocalParticles
     integer,dimension(nsp) :: nPartSpec
     integer,intent(inout)  :: nColl,nBadColl,iNumThirdBody

     integer :: iPart,idist,Rdmpos,temp
  
     ! Counters for local cell on this timestep. Added to global counters
     ! outside this routine.
     ! # of collisions
     nColl = 0
     ! # of bad collisions
     nBadcoll = 0
  
     ! RdmCollPairs simply holds the ordered indices of
     ! all particles in a cell (needed for randomization
     ! of pairs)
     ! NOTE ! RdmCollPairs is actually much larger than nLocalParticles 
     ! NOTE ! (currently hard-coded to 1 million) to avoid re-allocating
     ! NOTE ! memory constantly. The same is true for ParIniArr.
     do iPart=1, nLocalParticles
        RdmCollPairs(iPart) = iPart
     enddo

     ! FREE-BEE -> Fisher-Yates shuffle
     ! One efficient way to randomly shuffle our list of indices. If
     ! our shuffling is truly random, then any two consecutive pairs in the 
     ! list are also random. So you only have to grab two consecutive
     ! elements of RdmCollPairs (which refer to elements of ParIniArr) to 
     ! test random pairs.

     ! This is the Fisher-Yates algorithm, which randomly shuffles the 
     ! indices of RdmCollPairs.
     ! NOTE ! Calling mgds_randu() gives a uniformly random number from 
     ! NOTE ! 0 to 1
     do iPart = 1, nLocalParticles
        idist = nLocalParticles-(iPart-1)
        Rdmpos = iPart + int(mgds_randu() * idist )
        temp = RdmCollPairs(iPart)
        RdmCollPairs(iPart) = RdmCollPairs(Rdmpos)
        RdmCollPairs(Rdmpos) = temp
     enddo

     ! Collide particles here!
     ! This collision loop will look like:
     ! do iPair = 1,2*max_pairs,2
     !    ! Get inddices of pairs
     !    P1 =  RdmCollPairs(ipair)
     !    P2 =  RdmCollPairs(ipair+1)
     !    Then, you can access the particles properties of the 
     !    shuffled indices.
     !    vx1 = parIniArr(Pindex1)%Vx
     !    vy1 = parIniArr(Pindex1)%Vy
     !    vz1 = parIniArr(Pindex1)%Vz
     !    vx2 = parIniArr(Pindex2)%Vx
     !    vy2 = parIniArr(Pindex2)%Vy
     !    vz2 = parIniArr(Pindex2)%Vz
     !  
     !    Then, test for a collision (somehow....)
     !    If we collide, Ncoll = Ncoll+1
     ! end do
  
     ! Sample properties based on sampling parameters and current
     ! time step. 
     ! isam = 0 -> no sampling
     ! isam = 1 -> continue sampling from last run
     ! NOTE ! Regardless of issm, with isam = 1 we start sampling 
     ! NOTE ! immediately. issm = 0 (forced) regardless of input deck
     ! isam = 2 -> Kill old sampling, and start sampling after issm
     !             iterations.
     if(isam.gt.0 .and. nci.gt.issm) then
        call mgds_dsmc_AddSamples(icell,nLocalParticles)
     endif


     ! If we wanted to do something with the particles in the cell we 
     ! flagged the user initialization (see mgds_USER_Init in extra_utils.f90), 
     ! uncomment this section and do so here.
     !if(icell .eq. mycell_id)then
     !   ! Do stuff
     !endif
  
  end subroutine mgds_collide

!=======================================================================
!  Subroutine : mgds_dsmc_AddSamples
!-----------------------------------------------------------------------
!  Purpose :
!    Subroutine to sample particle data in a cell
!
!  Arguments :
!    icell     [in]  cell number
!    num_pcell [in]  number of dsmc particles in this cell
!=======================================================================
  subroutine mgds_dsmc_AddSamples(icell,num_pcell)

     use mgds_geometry
     use mgds_particles
     use mgds_dsmc

     Implicit None

     integer :: icell,num_pcell
     integer :: isp,n

     do n = 1,num_pcell

        isp = parIniArr(n) % spectype
        if(isp .gt. 0) then

           celld(icell) %sums( 1,isp) = &
              celld(icell) %sums( 1,isp) + 1.0d0
           celld(icell) %sums( 2,isp) = &
              celld(icell) %sums( 2,isp) + parIniArr(n) %vx
           celld(icell) %sums( 3,isp) = &
              celld(icell) %sums( 3,isp) + parIniArr(n) %vy 
           celld(icell) %sums( 4,isp) = &
              celld(icell) %sums( 4,isp) + parIniArr(n) %vz
           celld(icell) %sums( 5,isp) = &
              celld(icell) %sums( 5,isp) + (parIniArr(n) %vx)**2
           celld(icell) %sums( 6,isp) = &
              celld(icell) %sums( 6,isp) + (parIniArr(n) %vy)**2
           celld(icell) %sums( 7,isp) = &
              celld(icell) %sums( 7,isp) + (parIniArr(n) %vz)**2
           celld(icell) %sums( 8,isp) = &
              celld(icell) %sums( 8,isp) + parIniArr(n) %Erot
           celld(icell) %sums( 9,isp) = &
              celld(icell) %sums( 9,isp) + parIniArr(n) %Evib
        endif
     enddo
  end subroutine mgds_dsmc_AddSamples
