      PROGRAM LinPeriodicityTest
        IMPLICIT NONE
        INTEGER*4 :: i,j,k,l,ix,iy,iz,ip,ie,iee,PerCells,BB,iBin(3),
     >               x_per_flag,y_per_flag, z_per_flag,ppiclf_neltbbb,
     >               nfx,nfy,nfz, fcount,ii ,jj, kk, nBin(3)
 
        REAL*8    :: d2l,d2i,shift(3),shifti(3),centeri(3,3390),
     >               centerPer(4,3390), xcoord, ycoord, zcoord,
     >               Max_EleLen(3),binlength(3),BoundDist(3),
     >               CurrentbinBMin(3), binbound(3),
     >               CurrentbinBMax(3), ppiclf_binb(6),ppiclf_bins_dx(3)

        ppiclf_neltbbb = 3390
        ppiclf_binb(1) = -20.0D0
        ppiclf_binb(2) = 25.0D0
        ppiclf_binb(3) = 0.0D0
        ppiclf_binb(4) = 45.0D0
        ppiclf_binb(5) = 20.0D0
        ppiclf_binb(6) = 65.0D0
        nBin(1) = 1
        nBin(2) = 1
        nBin(3) = 1
        nfx = 15
        nfy = 15
        nfz = 15
        x_per_flag = 1
        y_per_flag = 1
        z_per_flag = 1

        DO l = 1,3
          ppiclf_bins_dx(l) = (ppiclf_binb(2*l) - ppiclf_binb(2*l-1))
     >                        /nBin(l)
          Max_EleLen(l) = 1.0D0
        END DO
        fcount = 1
        DO i = 1,14
          IF(i<3) THEN
            xcoord = ppiclf_binb(2) - 0.5*Max_EleLen(1) !Periodic
          ELSE
            xcoord = ppiclf_binb(1) !in current bin
          END IF
          IF(i < 10) THEN
            ycoord = ppiclf_binb(3) !in current bin
          ELSE
            ycoord = ppiclf_binb(4) !Periodic
          END IF
          IF(i<9) THEN
            zcoord = ppiclf_binb(5) + i/10*Max_EleLen(3) !in current bin
          ELSE
            zcoord = ppiclf_binb(6) + i/20*Max_EleLen(3) !Periodic
          END IF
          centeri(1,fcount) = xcoord 
          centeri(2,fcount) = ycoord
          centeri(3,fcount) = zcoord
          fcount = fcount + 1
        END DO

        ! Cells in ppiclf domain when nBins = 9
        DO i = 15,nfx+15
          xcoord = ppiclf_binb(1) + 10*Max_EleLen(1) + 
     >              (i-1)*Max_EleLen(1)
          DO j = 1,nfy
            ycoord = ppiclf_binb(3) +10*Max_EleLen(2) + 
     >               (j-1)*Max_EleLen(2)
            DO k = 1,nfz 
              zcoord = ppiclf_binb(5) +10*Max_EleLen(3) + 
     >               (k-1)*Max_EleLen(3)
              centeri(1,fcount) = xcoord
              centeri(2,fcount) = ycoord
              centeri(3,fcount) = zcoord
              fcount = fcount + 1
            END DO
          END DO
        END DO 
        ! *** Linear Periodicity Start ***
        ! Copy cell in centerPer array if linear periodicity is turned on
        ! and fluid cell is one of two layers of cells near boundary.
        ! centerPer(1-x;2-y;3-z;4-FluidCellNumber , Number of Periodic Cells)

        ! Set multiple of maximum element length in each direction
        DO l = 1,3
          binlength(l) = ppiclf_binb(2*l)-ppiclf_binb(2*l-1)
        END DO
        
        ! Find current bin boundaries
        DO l = 1,3
          iBin(l) = 0 
          CurrentbinBMin(l) = ppiclf_bins_dx(l)*iBin(l)+
     >                        ppiclf_binb(2*l-1)
          CurrentbinBMax(l) = ppiclf_bins_dx(l)*
     >                        (iBin(l)+1)+ppiclf_binb(2*l-1)
        END DO
        WRITE(*,*) 'BinbMin x:', CurrentbinBMin(1)
        WRITE(*,*) 'BinbMax x:', CurrentbinBMax(1)
        WRITE(*,*) 'x periodicity:', x_per_flag
        WRITE(*,*) 'y periodicity:', y_per_flag
        WRITE(*,*) 'z periodicity:', z_per_flag
        WRITE(*,*) 'Number of Bins (x,y,z):', nBin(1),nBin(2),nBin(3)
        ! if ix/y/z ==2, then periodicity on in that dimension
        ix = 1
        iy = 1
        iz = 1
        IF (x_per_flag .EQ. 1) ix = 2
        IF (y_per_flag .EQ. 1) iy = 2
        IF (z_per_flag .EQ. 1) iz = 2
        ! Number of copied fluid cells due to periodicity 
        PerCells = 0
        IF(ix.EQ.2 .OR. iy.EQ.2 .OR. iz.EQ.2) THEN
          ! At least one dimension is periodic
          ! Loop through twice to look at max and min bin boundaries
          DO BB = 1,2
            DO l = 1,3
              IF(BB .EQ. 1) THEN
                ! Comparing cell with minimum bin domain boundary
                shift(l)    =  binlength(l)
                binbound(l) =  ppiclf_binb((2*l)-1)
              ELSE
                ! Comparing cell with maximum bin domain boundary
                shift(l)    = -binlength(l)
                binbound(l) =  ppiclf_binb(2*l)
              END IF 
            END DO !l
            ! Loop through every cell mapped to bin
            DO ie = 1,ppiclf_neltbbb
              ! x periodicity loop
              DO ii = 1,ix
                IF(ii.EQ.2) THEN
                  BoundDist(1) = ABS(centeri(1,ie)-binbound(1))
                  ! 1.51 since comparing to cell center location
                  IF(BoundDist(1) < 1.51*Max_EleLen(1)) THEN
                    shifti(1) = shift(1)
                  ELSE
                    CYCLE
                  END IF
                ELSE
                  shifti(1) = 0.0D0
                END IF
                ! y periodicity loop
                DO jj = 1,iy
                  IF(jj.EQ.2) THEN
                    BoundDist(2) = ABS(centeri(2,ie)-binbound(2))
                    IF(BoundDist(2) < 1.51*Max_EleLen(2)) THEN
                      shifti(2) = shift(2)
                    ELSE
                      CYCLE 
                    END IF
                  ELSE
                    shifti(2) = 0.0D0
                  END IF
                  ! z periodicity loop
                  DO kk = 1,iz
                    IF(kk.EQ.2) THEN
                      BoundDist(3) = ABS(centeri(3,ie)-binbound(3))
                      IF(BoundDist(3) < 1.51*Max_EleLen(3)) THEN
                        shifti(3) = shift(3)
                      ELSE
                        CYCLE
                      END IF
                    ELSE
                      shifti(3) = 0.0D0
                    END IF
                    !Making sure we only duplicate periodic cells
                    IF(ii.EQ.1 .AND. jj.EQ.1 .AND. kk.EQ.1)CYCLE
                    ! Save periodic cell copy 
                    PerCells = PerCells + 1
                    ! Save fluid cell number
                    centerPer(4,PerCells) = ie
                    DO l = 1,3
                      ! shifti will be nonzero if dimension is periodic
                      ! and close to ppiclf domain boundary
                      centerPer(l,PerCells) = centeri(l,ie) + shifti(l)

                      ! Check if copied cell is outside of the
                      ! bin fluid mapping bounds in any dimension.
                      ! If it,then rewrite over this cell index in next loop cycle.
                      IF(centerPer(l,PerCells) .LT.
     >                (CurrentbinBMin(l) - 1.52*Max_EleLen(l))
     >                .OR. centerPer(l,PerCells) .GT.
     >                (CurrentbinBMax(l) + 1.52*Max_EleLen(l))) THEN 
                         PerCells = PerCells -1
                         EXIT
                      END IF
                    END DO !l
                  END DO !iz
                END DO !iy
              END DO !ix
            END DO !ie
          END DO !BB
        END IF !lin periodicity on
        ! *** Linear Periodicity End ***
        WRITE(*,*) 'Periodic cell copied information:'
        WRITE(*,*) 'Periodic cells:',PerCells
        DO i = 1,PerCells
          WRITE(*,*) centerPer(1,i), centerPer(2,i), centerPer(3,i), 
     >               centerPer(4,i)
        END DO
      END PROGRAM 
