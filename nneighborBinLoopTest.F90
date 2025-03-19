PROGRAM nneighborTest 
  IMPLICIT NONE
  
  INTEGER*4 :: part_num, i, j, k, l, overlap, x_b, y_b, z_b, nx_b, ny_b, nz_b,&
               found, tempbin, jpart, ibin, jbin, kbin, binLoop, partCount, totbin
  REAL*8    :: part_vol_fract, part_dia, part_vol, x_min, x_max, y_min, y_max, &
               rand, z_min, z_max, x_range, y_range, z_range, total_vol, pi, delta,&
               start, finish, temp
  REAL*8,    DIMENSION(:,:), ALLOCATABLE :: part_pos
  INTEGER*4, DIMENSION(:,:), ALLOCATABLE :: near
  INTEGER*4, DIMENSION(:,:), ALLOCATABLE :: bin
  INTEGER*4, DIMENSION(:),   ALLOCATABLE :: bincounter

!*******************************************************************
!*******************************
!***Set dimensional Parameter***
!*******************************
    pi             = 3.14159265359d0
    part_vol_fract = 0.10d0 
    part_dia       = 0.0001d0
    part_vol       = pi/6.0d0*part_dia**3

    x_min = -0.00300d0 + part_dia/2.0d0 
    x_max =  0.00300d0 - part_dia/2.0d0
    y_min = -0.00300d0 + part_dia/2.0d0
    y_max =  0.00300d0 - part_dia/2.0d0
    z_min = -0.00300d0 + part_dia/2.0d0
    z_max =  0.00300d0 - part_dia/2.0d0

    x_range = x_max - x_min
    y_range = y_max - y_min
    z_range = z_max - z_min
    total_vol = x_range*y_range*z_range
    part_num = INT(part_vol_fract*total_vol/part_vol)  
    WRITE(*,'("Total Number of Particles: ",i6)') part_num 
    ALLOCATE (part_pos(part_num,3)) 
    ALLOCATE (near(part_num,2))
!*******************************************************************
!***************************
!***GENERATE PARTICLE BED***
!***************************
      i = 1 
      main: DO 
        IF (i == part_num) EXIT main
        CALL RANDOM_NUMBER(rand)
        part_pos(i,1) = rand*x_range+x_min
        CALL RANDOM_NUMBER(rand)
        part_pos(i,2) = rand*y_range+y_min
        CALL RANDOM_NUMBER(rand)
        part_pos(i,3) = rand*z_range+z_min
        overlap = 0
        IF (i > 1) THEN
          checker: DO j = 1,i
            delta = ((part_pos(i,1) - part_pos(j,1))**2 + (part_pos(i,2) - part_pos(j,2))**2 + (part_pos(i,3) - part_pos(j,3))**2)
            IF (delta <= (part_dia)**2 .AND. i /= j) THEN
              overlap = 1
              EXIT checker
            END IF
          END DO checker
        END IF
        IF (overlap == 0) THEN
          i = i + 1
        END IF
      END DO main
!*******************************************************************
!***************************
!***BRUTE FORCE SEARCH***
!***************************!*******************************************************************
      ! Time to do a nearest neighbor search!
      CALL CPU_TIME(start)
      found = 0
      WRITE(*,'("Brute Force: Particle # and within x dps of Particle i")')
      DO i = 1,part_num
        DO j = 1,part_num
          near(j,1) = j 
          IF (i == j) THEN
            near(j,2) = -99 !This means it is over 5*dp away from ith particle
            CYCLE
          END IF
          ! Returns -99 if far away, or # of dp away.  Ex: near(j,2) = 2 means 
          ! jth particle is greater than 1*dp and less than 2*dp away from ith particle
          CALL BruteNN(i,j,part_dia,part_pos(i,1:3),part_pos(j,1:3),near(j,2))
          IF (i == part_num .AND. near(j,2) > 0) WRITE(*,*) near(j,1),near(j,2)
          IF (i == part_num .AND. near(j,2) > 0) found = found + 1 
        END DO
      END DO
      Call CPU_TIME(finish)
      WRITE(*,'("Number of particles within 5dp: ",i4)') found
      WRITE(*,'("Time for brute force method: ", f8.3)') (finish-start)
!***************************************************************************
! ****** Subbin Method
      ! Time to do a nearest neighbor search
      CALL CPU_TIME(start)
      WRITE(*,'("Subbin Method:")')
      ! l is the ratio of subbin lengths to particle diameter.
      ! l should be set to largest radius desired to search for neighbors
      l = 5    
      WRITE(*,'("************************************")')
      ! Find total number of subbins
      nx_b = FLOOR((x_range)/(l*part_dia)) 
      ny_b = FLOOR((y_range)/(l*part_dia)) 
      nz_b = FLOOR((z_range)/(l*part_dia))
      !total number of bins
      totbin = (nz_b+1)*(nx_b+1)*(ny_b+1)
      ALLOCATE (bincounter(0:(totbin-1)))
      ALLOCATE (bin(0:(totbin-1),part_num))
      near = 0
      tempbin = 0
      bincounter = 0
      WRITE(*,'("Number of x, y, and z subbins: ", i3,", ", i3,", ", i3)') nx_b+1,ny_b+1,nz_b+1
      WRITE(*,'("Subbin lengths equal ",i2," x particle diameter")') l
      
      ! Store each particle in the bin array
      DO i = 1,part_num
        ! find the subbin associated with the ith particle:
        x_b = FLOOR((part_pos(i,1) - x_min)/(l*part_dia))
        y_b = FLOOR((part_pos(i,2) - y_min)/(l*part_dia))
        z_b = FLOOR((part_pos(i,3) - z_min)/(l*part_dia))
        tempbin = x_b + y_b*nx_b + z_b*nx_b*ny_b
        !array tracking current index of each bin as it is built  
        bincounter(tempbin) = bincounter(tempbin) + 1
        ! adds particles to a subbin 
        bin(tempbin,bincounter(tempbin)) = i
      END DO
      found = 0
      WRITE(*,'("Subbin: Particle # and within x dps of Particle i")')
      DO i = 1,part_num
        ! find the subbin associated with the ith particle
        x_b = FLOOR((part_pos(i,1) - x_min)/(l*part_dia))
        y_b = FLOOR((part_pos(i,2) - y_min)/(l*part_dia))
        z_b = FLOOR((part_pos(i,3) - z_min)/(l*part_dia))
        tempbin = x_b + y_b*nx_b + z_b*nx_b*ny_b !subbin of ith particle
        ! partCount counts the number of particles being compared to ith particle
        partCount = 0
        DO ibin = 1,3     ! to look at -1,0,+1 x-dir subbins 
          DO jbin = 1,3   ! to look at -1,0,+1 y-dir subbins
            DO kbin = 1,3 ! to look at -1,0,+1 z-dir subbins
              !loops through 27 adjacent subbins
              binLoop = tempbin + (-2+ibin) + (-2+jbin)*(nx_b) + (-2+kbin)*(nx_b)*(ny_b) 
              ! Ensures binLoop index isn't outside of bin array
              if (binLoop > -1 .AND. binLoop <= totbin) THEN
                DO j = 1,bincounter(binLoop)
                  jpart = bin(binLoop,j)    ! particle index being compared to ith particle
                  partCount = partCount + 1 ! Counter for near array, to account for jparticles from other subbins
                  near(partCount,1) = jpart
                  ! Returns -99 if far away, or # of dp away.  Ex: near(j,2) = 2 means 
                  ! jth particle is greater than 1*dp and less than 2*dp away from ith particle
                  CALL BruteNN(i,jpart,part_dia,part_pos(i,1:3),part_pos(jpart,1:3),near(partCount,2))
                  IF (i == part_num .AND. near(partCount,2) > 0) found = found + 1
                  IF (i == part_num .AND. near(partCount,2) > 0) WRITE(*,*) near(partCount,1), near(partCount,2)
                END DO 
              END IF
            END DO
          END DO
        END DO
        IF (i ==part_num) WRITE(*,'("Number of particles within 5dp: ",i4)') found
      END DO
      CALL CPU_TIME(finish)
      WRITE(*,'("Total time for subbin method: ", f8.3)') (finish-start) 
    
    DEALLOCATE (bin)
    DEALLOCATE (bincounter)
    DEALLOCATE (near)
    DEALLOCATE (part_pos)
END PROGRAM nneighborTest



SUBROUTINE BruteNN(i,j,dia,ipos,jpos,iNear)
IMPLICIT NONE
! Input
  INTEGER*4, INTENT(IN)  :: i,j
  REAL*8,    INTENT(IN)  :: dia, ipos(3), jpos(3)    
  INTEGER*4, INTENT(OUT) :: iNear
! Internal
  REAL*8    :: deltai      = 0.0d0
  INTEGER*4 :: delta_ratio = - 99
! delta ratio is the distance between centers normalized by particle diameter  
    deltai = ((ipos(1) - jpos(1))**2 + (ipos(2) - jpos(2))**2 + (ipos(3) - jpos(3))**2)
    delta_ratio = CEILING(deltai/(dia)**2)
    IF (delta_ratio > 5) THEN
      iNear = -99 
    ELSE
      iNear  = delta_ratio
    END IF
END SUBROUTINE

