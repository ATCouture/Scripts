PROGRAM RandomPart
IMPLICIT NONE

INTEGER*4         :: part_num, i, j, overlap
CHARACTER(len=46) :: fname        
CHARACTER(len=2)  :: part_mat      
REAL*8 :: part_vol_fract, part_dia, part_vol, x_min, x_max, y_min, y_max, &
        rand, z_min, z_max, x_range, y_range, z_range, total_vol, pi, delta
REAL*8, DIMENSION(:), allocatable :: x_pos, y_pos, z_pos


REAL :: start, finish

CALL CPU_TIME(start)

  pi             = 3.14159265359d0
  part_vol_fract = 0.35d0 
  part_dia       = 0.0001d0
  part_vol       = pi/6.0d0*part_dia**3
  part_mat       = 'St'

  x_min = -0.006d0 + part_dia/2.0d0 
  x_max =  0.006d0 - part_dia/2.0d0
  y_min = -0.006d0 + part_dia/2.0d0
  y_max =  0.006d0 - part_dia/2.0d0
  z_min =  0.034d0 + part_dia/2.0d0
  z_max =  0.040d0 - part_dia/2.0d0

  x_range = x_max - x_min
  y_range = y_max - y_min
  z_range = z_max - z_min

  total_vol = pi*((x_max + part_dia/2.0d0)**2)*(z_range + part_dia)
  part_num = INT(part_vol_fract*total_vol/part_vol)  
  
  allocate (x_pos(1:part_num))
  allocate (y_pos(1:part_num))
  allocate (z_pos(1:part_num))
  
  i = 1 
  main: DO 
    IF (i == part_num) EXIT main
    CALL RANDOM_NUMBER(rand)
    x_pos(i) = rand*x_range+x_min
    CALL RANDOM_NUMBER(rand)
    y_pos(i) = rand*y_range+y_min
    CALL RANDOM_NUMBER(rand)
    z_pos(i) = rand*z_range+z_min
    overlap = 0
    IF (SQRT(x_pos(i)**2+y_pos(i)**2) > (x_max+part_dia/2.0d0)) overlap = 1
    IF (i > 1 .AND. overlap < 1) THEN
      checker: DO j = 1,i
        delta = SQRT((x_pos(i) - x_pos(j))**2 + (y_pos(i) - y_pos(j))**2 + (z_pos(i) - z_pos(j))**2)
        IF (delta <= part_dia .AND. i /= j) THEN
          overlap = 1
          EXIT checker
        END IF
      END DO checker
    END IF
    IF (overlap == 0) THEN
      i = i + 1
    END IF
  END DO main

! Printing points.dat file
! File Format:
!-----------------------------------------------------------------------------
!L1:       I16 [Number of Particles]
!L2-Npar: A16 [Part Material] 4F23.16 [x-coord] [y-coord] [z-coord] [Part Dia]
!-----------------------------------------------------------------------------

  fname    = 'C:\Users\avery\OneDrive\Desktop\35testpoints.dat'
  part_mat = 'St'
  OPEN(UNIT=11, FILE = fname, STATUS = 'NEW', ACTION = 'WRITE')
  WRITE(11,'(I16)') part_num
  DO i = 2, (part_num+1)
    WRITE(11,'(A16,4E23.16)') part_mat,x_pos(i-1),y_pos(i-1),z_pos(i-1),part_dia
  END DO

CALL CPU_TIME(finish)

WRITE(*,'("Time = ",f25.3," minutes.")') ((finish-start)/60)

END PROGRAM RandomPart
