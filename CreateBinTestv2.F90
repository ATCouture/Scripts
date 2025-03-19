      PROGRAM CreateBin
      IMPLICIT NONE
  
      INTEGER*4 :: nBin(3), nProc, l, totBin, nBinMax, ix, iy, iz, iBin(3), iBinTot, &
                   idealBin(3), k, m
                
      REAL*8    :: lBin(3), dxBin(3), filter
      LOGICAL :: check, changed
      CHARACTER(2) :: strg
      nProc   = 50
      lBin(1) = 6.0D0
      lBin(2) = 6.2D0
      lBin(3) = 6.1D0
      filter  = 0.1D0
      totBin  = 1
      WRITE(*,*) '********************************************'
      WRITE(*,*) ' User Input Values for bin test:'
      WRITE(*,*) 'Neighbor/filter width:', filter
      WRITE(*,*) 'Number of Processors:',nProc
      WRITE(*,*) 'Length of bin in x: ', lBin(1)
      WRITE(*,*) 'Length of bin in y: ', lBin(2)
      WRITE(*,*) 'Length of bin in z: ', lBin(3)
      WRITE(*,*) '********************************************'
      WRITE(*,*) ''

      ! Number of bins based on least bin surface area criteria
      nBin(1) = INT((nProc**(1.0D0/3.0D0))*(lBin(1)**(2.0D0/3.0D0))/ &
                ((lBin(2)**(1.0D0/3.0D0))*(lBin(3))**(1.0D0/3.0D0)))
    
      nBin(2) = INT((nProc**(1.0D0/3.0D0))*(lBin(2)**(2.0D0/3.0D0))/ &
                ((lBin(1)**(1.0D0/3.0D0))*(lBin(3))**(1.0D0/3.0D0)))
     
      nBin(3) = INT((nProc**(1.0D0/3.0D0))*(lBin(3)**(2.0D0/3.0D0))/ &
                ((lBin(2)**(1.0D0/3.0D0))*(lBin(1))**(1.0D0/3.0D0))) 
      totBin  = 0
      iBinTot = 0

      ! Filterwidth criteria check
      DO l = 1,3
        DO
          ! Checking with nBin + 1 will make sure all 3 iBin combinations 
          ! in DO loop below will meet filterwidth criteria. 
          IF ((nBin(l)+1)/lBin(l) > filter) EXIT
          nBin(l) = nBin(l) - 1 
          IF (nBin(l) < 0) THEN ! < 0 instead of <1 since we added 1 in if statement
            WRITE(*,*) 'Filterwidth criteria is violated. Particle domain is less than filterwidth'
            RETURN !will call ppiclf_exit subroutine in actual code
          END IF
        END DO
      END DO


      ! Since bin must be an integer, check -1, +0, +1 number of bins for each bin dimension
      ! ideal number of bins will be max value while less than number of processors.
      ! Will not check total bin value (cycle do loop) if filterwidth criteria is violated.
      DO ix = 1,3
        iBin(1) = nBin(1) + (ix-2)
        dxBin(1) = lBin(1)/iBin(1)
        IF(dxBin(1) < filter .OR. iBin(1) < 1) CYCLE
        DO iy = 1,3
          iBin(2) = nBin(2) + (iy-2)
          dxBin(2) = lBin(2)/iBin(2)
          IF(dxBin(2) < filter .OR. iBin(2) < 1) CYCLE
          DO iz = 1,3
            iBin(3) = nBin(3) + (iz-2)
            dxBin(3) = lBin(3)/iBin(3)
            IF(dxBin(3) < filter .OR. iBin(3) < 1) CYCLE
            iBinTot = iBin(1)*iBin(2)*iBin(3)
            IF(iBinTot > totBin .AND. iBinTot <= nProc) THEN
              totBin = 1
              DO l = 1,3
                idealBin(l) = iBin(l)
                totBin = totBin*idealBin(l)
              END DO
              ! This loop is to make sure the dimension with the longest
              ! dxBin gets more bins in the case where two dimensions are within 
              ! 1 bin division.
              ! The code may be hard to follow, but has been tested to confirm it works
              DO l = 1,2
                ! l = 1 compares x with y and x with z.
                ! l = 2 compares y with z
                k  = l + 1
                IF(l==2) k = 2
                IF(ABS(idealBin(l)-idealBin(k)) == 1) THEN
                  IF(lBin(l)/(iBin(l)-1) < lBin(k)/iBin(k)) THEN
                    idealBin(l) = idealBin(l) - 1
                    idealBin(k) = idealBin(k) + 1
                  END IF
                  IF(lBin(l)/iBin(l) > lBin(k)/(iBin(k)-1)) THEN
                    idealBin(l) = idealBin(l) + 1
                    idealBin(k) = idealBin(k) - 1
                  END IF
                END IF
                m = l + 2
                IF(l==2) m = 3
                IF(ABS(idealBin(l)-idealBin(m)) == 1) THEN
                  IF(lBin(l)/(idealBin(l)-1) < lBin(m)/idealBin(m)) THEN
                    idealBin(l) = idealBin(l) - 1
                    idealBin(m) = idealBin(m) + 1
                  END IF
                  IF(lBin(l)/idealBin(l) > lBin(m)/(idealBin(m)-1)) THEN
                    idealBin(l) = idealBin(l) + 1
                    idealBin(m) = idealBin(m) - 1
                  END IF
                END IF
              END DO !l
            END IF
          END DO !iz
        END DO !iy
      END DO !ix
            
      DO l = 1,3
        IF (l .EQ. 1) strg = 'x'
        IF (l .EQ. 2) strg = 'y'
        IF (l .EQ. 3) strg = 'z'
        nBin(l) = idealBin(l)
        dxBin(l) = lBin(l)/nBin(l)
        WRITE(*,*) 'Number of bins in ',strg,': ',nBin(l)
      END DO
      WRITE(*,*) 'Total Bins:',totBin
      WRITE(*,*) 'Number of Processors:',nProc
      changed = .FALSE.

      ! Loop ensuring number of bins less than number of processors.
      ! Shouldn't ever be violated, but will subtract from dimension
      ! with most bins to decrease totBin by smallest value if violated.
      DO
        check = .TRUE. 
        IF (totBin > nProc) THEN
          check = .FALSE.
          l = MAXLOC(nBin(:),1)
          nBin(l) = nBin(l) - 1
          dxBin(l) = lBin(l)/nBin(l)
          changed = .TRUE.
          totBin = 1
          DO l = 1,3
            totBin = totBin*nBin(l)
          END DO
        END IF
        IF (check) EXIT
      END DO

      ! Loop to see if we can add one to dimension with largest number of bins
      ! Shouldn't ever be needed, but good test for now.
      nBinMax = MAXLOC(nBin(:),1)
      DO
        IF ((totbin/nBin(nBinMax))*(nBin(nBinMax)+1) < nProc) THEN
          nBin(nBinMax) = nBin(nBinMax)+1
          dxBin(nBinMax) = lBin(nBinMax)/nBin(nBinMax)
          totBin = 1
          DO l = 1,3
            totBin = totBin*nBin(l)
          END DO
          changed = .TRUE.
        ELSE
          EXIT
        END IF
      END DO

      WRITE(*,*) ''
      IF (changed) THEN
        WRITE(*,*)  '********************************************'
        WRITE(*,*)  'Number of bins Changed from inital formula'
        WRITE(*,*)  '********************************************'
        DO l = 1,3
          IF (l .EQ. 1) strg = 'x'
          IF (l .EQ. 2) strg = 'y'
          IF (l .EQ. 3) strg = 'z'
          WRITE(*,*) 'Number of bins in ',strg,': ',nBin(l)
          WRITE(*,*) 'Length of bin dx in ',strg,': ',dxBin(l)
        END DO  
        WRITE(*,*) 'Total Bins:',totBin
        WRITE(*,*) 'Number of Processors:',nProc
      END IF
      END PROGRAM
