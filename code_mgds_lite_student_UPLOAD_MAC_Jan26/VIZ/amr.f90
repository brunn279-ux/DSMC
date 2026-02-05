!=======================================================================
!  Subroutine : mgdsviz_AMR
!-----------------------------------------------------------------------
!  Purpose :
!    Subroutine to perform Adaptive Mesh Refinement (AMR)
! 
!  Remarks :
!   This methodology should be consistent with the routine
!   "mgds_Geom_MeshCube()
!=======================================================================
  subroutine mgdsviz_AMR()

     use mgds_inputvars
     use flow_field_output_data
     use mgds_AMR

     Implicit None

     integer :: i,i1,j1,k1, i2,j2,k2, i3,j3,k3, n,n1,n2,n3,mn2
     integer :: im2,jm2,km2, im3,jm3,km3
     integer :: ni2,nj2,nk2, nn2,nn3, no2
     integer :: myi3,myj3,myk3
     integer :: ierr

     integer,dimension(3) ::  ips,jps

     real*8 :: rlam, dref, gsp, rn, rlnew
     real*8 :: rlmx,rlmn,mfp_max,mfp_min
     real*8 :: x0,y0,z0,x1,y1,z1
     real*8 :: dxl,dyl,dzl, dxs,dys,dzs
     real*8 :: dx2,dy2,dz2, dx3,dy3,dz3, tmp
     real*8 :: min_dt,mfp_fraction
     real*8,dimension(3) :: xno,xne,xxo,xxe, xxxo,xxxe

     real*8,parameter  :: eps = 1d-9

     ! It is useful need to find the smallest mean-collision time
     min_dt = 9.99d99

     ! 3-level Cartesian hierarchical loop
     x0 = 0.0d0
     y0 = 0.0d0
     z0 = 0.0d0
     n1 = 0
     n2 = 0
     n3 = 0
     do k1 = 1,kmL1
     do j1 = 1,jmL1
     do i1 = 1,imL1
        n1 = n1 + 1
        do k2 = 1,kmL2(i1,j1,k1)
        do j2 = 1,jmL2(i1,j1,k1)
        do i2 = 1,imL2(i1,j1,k1)
           n2 = n2 + 1
           do k3 = 1,kmL3(n2)
           do j3 = 1,jmL3(n2)
           do i3 = 1,imL3(n2)
              n3 = n3 + 1

              if(mfp_store(n3,2).lt.min_dt) min_dt = mfp_store(n3,2)
           enddo
           enddo
           enddo
        enddo
        enddo
        enddo
     enddo
     enddo
     enddo

     write(*,*) 'Smallest mean collision time: ',min_dt

     !
     ! ******** Create a new grid *********
     ! First, collect largest lambda and build L2 cells
     ! Second, collect smallest lambda and build L3 cells
      
     ! For now, L1 cells do not change.
     WRITE(*,*) 'Level one cells:', imL1,jmL1,kmL1

     allocate(newL1(imL1,jmL1,kmL1),&
             nimL2(imL1,jmL1,kmL1),&
             njmL2(imL1,jmL1,kmL1),&
             nkmL2(imL1,jmL1,kmL1),&
             nidL2(imL1,jmL1,kmL1), STAT=ierr)
     if(ierr.ne.0) then
        write(*,*) 'Could not allocate new L1 grid structures'
        STOP
     endif

     WRITE(*,*) 'Done allocating new L2 sizing arrays'

     x0 = 0.0d0
     y0 = 0.0d0
     z0 = 0.0d0

     ! Start by making new L2 cells in each L1 cell
     n1 = 0
     n2 = 0
     nn2 = 0
     do k1 = 1,kmL1
     do j1 = 1,jmL1
     do i1 = 1,imL1
        n1 = n1 + 1
        x1 = x0 + dble(i1-1)*bdx1
        y1 = y0 + dble(j1-1)*bdy1
        z1 = z0 + dble(k1-1)*bdz1

        ! Calculate new sizes of the L2 cells

        ! Find min and max mean-free path in L1 cell of original grid
        rlmn = 9.9d99
        rlmx = -9.99d99
        do k2 = 1,kmL2(i1,j1,k1)
        do j2 = 1,jmL2(i1,j1,k1)
        do i2 = 1,imL2(i1,j1,k1)
           n2 = n2 + 1
           do k3 = 1,kmL3(n2)
           do j3 = 1,jmL3(n2)
           do i3 = 1,imL3(n2)
              n3 = cellsL1(i1,j1,k1) %PLevelTwoCellInfo(i2,j2,k2)&
                                     %PLevelThreeCellInfo(i3,j3,k3)&
                                     %icelld
              rlam = mfp_store(n3,1)
              if(rlam.lt.rlmn) then
                 rlmn = rlam
              endif
              if(rlam .gt. rlmx)then
                 rlmx = rlam
              endif
           enddo
           enddo
           enddo
        enddo
        enddo
        enddo

        if(rlmn .lt. mfp_min) rlmn = mfp_min
        if(rlmn .gt. mfp_max) rlmn = mfp_max

        ! For now, L2 sizing is a hard-coded choice
        dxl = sqrt(bdx1*rlmn)
        dyl = sqrt(bdy1*rlmn)
        dzl = sqrt(bdz1*rlmn)

        tmp = bdx1/dxl
        im2 = ceiling(tmp+eps)
        dx2 = bdx1/dble(im2)

        tmp = bdy1/dyl
        jm2 = ceiling(tmp+eps)
        dy2 = bdy1/dble(jm2)

        tmp = bdz1/dzl
        km2 = ceiling(tmp+eps)
        dz2 = bdz1/dble(km2)

        if(flag_2D)then
           ! If 2-D, only allow one cell in z direction
           km2 = 1
           dz2 = bdz1/dble(km2)
        endif

        !--- sizes of L2 cells have been calculated above

        allocate(newL1(i1,j1,k1) %PLevelTwoCellInfo(im2,jm2,km2), STAT=ierr)
        if(ierr.ne.0) then
           write(*,*) 'Could not allocate new L2 grid structures'
           write(*,*) 'i1,j1,k1=',i1,j1,k1
           STOP
        endif
        nidL2(i1,j1,k1) = nn2     ! offset of L2 cell in L2 size arrays
        nimL2(i1,j1,k1) = im2
        njmL2(i1,j1,k1) = jm2
        nkmL2(i1,j1,k1) = km2
        newL1(i1,j1,k1) %dChildCellX = dx2
        newL1(i1,j1,k1) %dChildCellY = dy2
        newL1(i1,j1,k1) %dChildCellZ = dz2

        do nk2 = 1,nkmL2(i1,j1,k1)
        do nj2 = 1,njmL2(i1,j1,k1)
        do ni2 = 1,nimL2(i1,j1,k1)
           newL1(i1,j1,k1) %PLevelTwoCellInfo(ni2,nj2,nk2) %DeltaX2 = dx2
           newL1(i1,j1,k1) %PLevelTwoCellInfo(ni2,nj2,nk2) %DeltaY2 = dy2
           newL1(i1,j1,k1) %PLevelTwoCellInfo(ni2,nj2,nk2) %DeltaZ2 = dz2
        enddo
        enddo
        enddo

        nn2 = nn2 + nimL2(i1,j1,k1)*njmL2(i1,j1,k1)*nkmL2(i1,j1,k1)
     enddo
     enddo
     enddo

     write(*,*) 'Created new L2 cells:',nn2
     nncL2 = nn2

!    Make new L3 cells in every L2 cell of each L1 cell
     allocate(nimL3(nncL2),&
              njmL3(nncL2),&
              nkmL3(nncL2),&
              nidL3(nncL2), STAT=ierr)
     if(ierr.ne.0) then
        write(*,*) 'Could not allocate new L3 integer arrays'
        STOP
     endif
     nimL3(:) = -1111111
     njmL3(:) = -2222222
     nkmL3(:) = -3333333
     nidl3(:) = -4444444

     myi3 = 0  !--- maximum L2 subdivision to L3
     myj3 = 0
     myk3 = 0

     n1 = 0
     n2 = 0
     nn3 = 0
     do k1 = 1,kmL1
     do j1 = 1,jmL1
     do i1 = 1,imL1
        n1 = n1 + 1
        x1 = x0 + dble(i1-1)*bdx1
        y1 = y0 + dble(j1-1)*bdy1
        z1 = z0 + dble(k1-1)*bdz1
        nn2 = nidL2(i1,j1,k1)

        ! Loop over new L2 cells to pick up min lamda for each one
        do nk2 = 1,nkmL2(i1,j1,k1)
        do nj2 = 1,njmL2(i1,j1,k1)
        do ni2 = 1,nimL2(i1,j1,k1)
           nn2 = nn2 + 1

           ! We need to calculate L3 cell sizes here
     
           dx2 = newL1(i1,j1,k1) %PLevelTwoCellInfo(ni2,nj2,nk2) %DeltaX2
           dy2 = newL1(i1,j1,k1) %PLevelTwoCellInfo(ni2,nj2,nk2) %DeltaY2
           dz2 = newL1(i1,j1,k1) %PLevelTwoCellInfo(ni2,nj2,nk2) %DeltaZ2
     
     !     Bounds of new L2 cell
           xno(1) = x1 + dble(ni2-1)*dx2
           xno(2) = y1 + dble(nj2-1)*dy2
           xno(3) = z1 + dble(nk2-1)*dz2
           xne(1) = xno(1) + dx2
           xne(2) = xno(2) + dy2
           xne(3) = xno(3) + dz2
     
     !     Find the smallest MFP in old L2 cells
           rlmn = 9.9d99
     
     !     Loop over old L2 cells in L1 cell
           mn2 = idL2(i1,j1,k1)
           do k2 = 1,kmL2(i1,j1,k1)
           do j2 = 1,jmL2(i1,j1,k1)
           do i2 = 1,imL2(i1,j1,k1)
              mn2 = mn2 + 1
              xxo(1) = x1 + dble(i2-1)*cellsL1(i1,j1,k1) %dChildCellX
              xxo(2) = y1 + dble(j2-1)*cellsL1(i1,j1,k1) %dChildCellY
              xxo(3) = z1 + dble(k2-1)*cellsL1(i1,j1,k1) %dChildCellZ
              xxe(1) = xxo(1) + cellsL1(i1,j1,k1) %dChildCellX
              xxe(2) = xxo(2) + cellsL1(i1,j1,k1) %dChildCellY
              xxe(3) = xxo(3) + cellsL1(i1,j1,k1) %dChildCellZ
     
              ips(:) = 0
              ! Segment test
              if( ((xno(1)+eps.ge.xxo(1)) .and. (xno(1)-eps.le.xxe(1))) .or.&
                  ((xne(1)+eps.ge.xxo(1)) .and. (xne(1)-eps.le.xxe(1)))) ips(1) = 1
              if( ((xno(2)+eps.ge.xxo(2)) .and. (xno(2)-eps.le.xxe(2))) .or.&
                  ((xne(2)+eps.ge.xxo(2)) .and. (xne(2)-eps.le.xxe(2)))) ips(2) = 1
              if( ((xno(3)+eps.ge.xxo(3)) .and. (xno(3)-eps.le.xxe(3))) .or.&
                  ((xne(3)+eps.ge.xxo(3)) .and. (xne(3)-eps.le.xxe(3)))) ips(3) = 1
              ! Complete encapsulation test
              if((xno(1)-eps.le.xxo(1)) .and. (xne(1)+eps.ge.xxe(1))) ips(1) = 1
              if((xno(2)-eps.le.xxo(2)) .and. (xne(2)+eps.ge.xxe(2))) ips(2) = 1
              if((xno(3)-eps.le.xxo(3)) .and. (xne(3)+eps.ge.xxe(3))) ips(3) = 1

              ! Check if new cell is inside an old cell
              if((xno(1)+eps.ge.xxo(1)) .and. (xne(1)-eps.le.xxe(1))) ips(1) = 1
              if((xno(2)+eps.ge.xxo(2)) .and. (xne(2)-eps.le.xxe(2))) ips(2) = 1
              if((xno(3)+eps.ge.xxo(3)) .and. (xne(3)-eps.le.xxe(3))) ips(3) = 1

              ! All directions pass test; passing means overlapping cells
              if((ips(1) .ne.0) .and. (ips(2).ne.0) .and. (ips(3).ne.0) ) then
                 n3 = idL3(mn2)
                 do k3 = 1,kmL3(mn2)
                 do j3 = 1,jmL3(mn2)
                 do i3 = 1,imL3(mn2)
                    n3 = n3 + 1
        
                    xxxo(1) = xxo(1) + dble(i3-1)*&
                              cellsL1(i1,j1,k1) %PLevelTwoCellInfo(i2,j2,k2)&
                                                %dChildCellX
                    xxxo(2) = xxo(2) + dble(j3-1)*&
                              cellsL1(i1,j1,k1) %PLevelTwoCellInfo(i2,j2,k2)&
                                                %dChildCellY
                    xxxo(3) = xxo(3) + dble(k3-1)*&
                              cellsL1(i1,j1,k1) %PLevelTwoCellInfo(i2,j2,k2)&
                                                %dChildCellZ
                    xxxe(1) = xxxo(1) +&
                              cellsL1(i1,j1,k1) %PLevelTwoCellInfo(i2,j2,k2)&
                                                %dChildCellX
                    xxxe(2) = xxxo(2) +&
                              cellsL1(i1,j1,k1) %PLevelTwoCellInfo(i2,j2,k2)&
                                                %dChildCellY
                    xxxe(3) = xxxo(3) +&
                              cellsL1(i1,j1,k1) %PLevelTwoCellInfo(i2,j2,k2)&
                                                %dChildCellZ
        
                    jps(:) = 0
                    ! Segment test
                    if( ((xno(1)+eps.ge.xxxo(1)) .and. (xno(1)-eps.le.xxxe(1))) .or.&
                        ((xne(1)+eps.ge.xxxo(1)) .and. (xne(1)-eps.le.xxxe(1))) ) jps(1)=1
                    if( ((xno(2)+eps.ge.xxxo(2)) .and. (xno(2)-eps.le.xxxe(2))) .or.&
                        ((xne(2)+eps.ge.xxxo(2)) .and. (xne(2)-eps.le.xxxe(2))) ) jps(2)=1
                    if( ((xno(3)+eps.ge.xxxo(3)) .and. (xno(3)-eps.le.xxxe(3))) .or.&
                        ((xne(3)+eps.ge.xxxo(3)) .and. (xne(3)-eps.le.xxxe(3))) ) jps(3)=1
                    ! Complete encapsulation test
                    if((xno(1)-eps.le.xxxo(1)) .and. (xne(1)+eps.ge.xxxe(1))) jps(1) = 1
                    if((xno(2)-eps.le.xxxo(2)) .and. (xne(2)+eps.ge.xxxe(2))) jps(2) = 1
                    if((xno(3)-eps.le.xxxo(3)) .and. (xne(3)+eps.ge.xxxe(3))) jps(3) = 1
                    
                    if((xno(1)+eps.ge.xxxo(1)) .and. (xne(1)-eps.le.xxxe(1))) jps(1) = 1
                    if((xno(2)+eps.ge.xxxo(2)) .and. (xne(2)-eps.le.xxxe(2))) jps(2) = 1
                    if((xno(3)+eps.ge.xxxo(3)) .and. (xne(3)-eps.le.xxxe(3))) jps(3) = 1
        
                    ! All directions test; passing means overlapping cells
                    if((jps(1) .ne.0) .and. (jps(2).ne.0) .and. (jps(3).ne.0)) then
        
                       n = cellsL1(i1,j1,k1) %PLevelTwoCellInfo(i2,j2,k2)&
                          %PLevelThreeCellInfo(i3,j3,k3)&
                          %icelld
                       rlam = mfp_store(n,1)
                       if(rlmn.gt.rlam) rlmn = rlam
                    endif
                 enddo   ! old L3 cells in old L2
                 enddo
                 enddo
              endif
           enddo   !--- sweeps of old L2 cells in L1
           enddo
           enddo

           ! User defined bounds of minimum mean-free path
           if(rlmn .lt. mfp_min) rlmn = mfp_min
           if(rlmn .gt. mfp_max) rlmn = mfp_max

           ! We know the smallest lamda in this new L2 cell
           ! Need to multiply times fraction of mfp we want
           
           dxs = rlmn
           dys = rlmn
           dzs = rlmn
           if(dxs.gt.dx2) dxs = dx2
           if(dys.gt.dy2) dys = dy2
           if(dzs.gt.dz2) dzs = dz2
     
     !     Create divisions for new L3 cells within this new L2 cell
           tmp = dx2/dxs
           im3 = ceiling(tmp + eps)
           dx3 = dx2/dble(im3)
     
           tmp = dy2/dys
           jm3 = ceiling(tmp + eps)
           dy3 = dy2/dble(jm3)
     
           tmp = dz2/dzs
           km3 = ceiling(tmp + eps)
           dz3 = dz2/dble(km3)
     
           if(flag_2D)then
              ! Force L3 cell size to 1
              km3 = 1
              dz3 = dz2/dble(km3)
           endif
           !--- sizes of L3 cells have been calculated above
     
           allocate(newL1(i1,j1,k1) %PLevelTwoCellInfo(ni2,nj2,nk2)&
              %PLevelThreeCellInfo(im3,jm3,km3), STAT=ierr)
           if(ierr.ne.0) then
               write(*,*) 'Could not allocate new L3 grid structures'
               write(*,*) 'i1,j1,k1=',i1,j1,k1
               write(*,*) 'ni2,nj2,nk2=',ni2,nj2,nk2
               STOP
           endif
           nidL3(nn2) = nn3    !! offset of L3 cells in the L3 size array
           nimL3(nn2) = im3
           njmL3(nn2) = jm3
           nkmL3(nn2) = km3
     
           newL1(i1,j1,k1) %PLevelTwoCellInfo(ni2,nj2,nk2) %dChildCellX = dx3
           newL1(i1,j1,k1) %PLevelTwoCellInfo(ni2,nj2,nk2) %dChildCellY = dy3
           newL1(i1,j1,k1) %PLevelTwoCellInfo(ni2,nj2,nk2) %dChildCellZ = dz3
     
           if(myi3.lt.im3) myi3 = im3   !-- keeping track of subdivision
           if(myj3.lt.jm3) myj3 = jm3
           if(myk3.lt.km3) myk3 = km3
           do k3 = 1,km3
           do j3 = 1,jm3
           do i3 = 1,im3
              nn3 = nn3 + 1
     
              newL1(i1,j1,k1) %PLevelTwoCellInfo(ni2,nj2,nk2)&
                              %PLevelThreeCellInfo(i3,j3,k3) %DeltaX3 = dx3
              newL1(i1,j1,k1) %PLevelTwoCellInfo(ni2,nj2,nk2)& 
                              %PLevelThreeCellInfo(i3,j3,k3) %DeltaY3 = dy3
              newL1(i1,j1,k1) %PLevelTwoCellInfo(ni2,nj2,nk2)& 
                              %PLevelThreeCellInfo(i3,j3,k3) %DeltaZ3 = dz3
              newL1(i1,j1,k1) %PLevelTwoCellInfo(ni2,nj2,nk2)& 
                              %PLevelThreeCellInfo(i3,j3,k3) %GlobalIndex =nn3
     
              ! Reset SigmaCrMax
              newL1(i1,j1,k1) %PLevelTwoCellInfo(ni2,nj2,nk2)& 
                                     %PLevelThreeCellInfo(i3,j3,k3)&
                                      %sigcrmax3 = 1.0d-16
           enddo
           enddo
           enddo
        enddo    ! new L2 cells in this L1
        enddo
        enddo
     nn2 = nn2 + nimL2(i1,j1,k1)*njmL2(i1,j1,k1)*nkmL2(i1,j1,k1)
     enddo    ! L1 cells
     enddo
     enddo

     write(*,*) 'Created new L3 cells:',nn3
     write(*,*) 'Maximum L3 divisions:',myi3,myj3,myk3 

     nncL3 = nn3
  end subroutine 

