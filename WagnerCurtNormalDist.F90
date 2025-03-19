PROGRAM WagnerPartCurtain
IMPLICIT NONE

INTEGER*4          :: part_num_max, part_num, i, j, k, nslices 
CHARACTER(len=200) :: fname        
CHARACTER(len=5)   :: part_mat 

REAL*8 :: pvf, part_dia, part_vol, x_min, x_max, y_min, y_max, x_center, &
          rand1, rand2, z_min, z_max, x_range, y_range, z_range, total_vol, pi, delta, &
          x_bmin, x_bmax, dx, pvf_max, pvf_min, pvf_dist, max_part_dia, delta_part_dia, &
          pVolTotal, y_bmin, y_bmax, z_bmin, z_bmax, x_minSlice, pVolTarget, &
          min_part_vol, min_part_dia, x_mean, pvf_mean, z_dist, z_dist2,x_sigma

REAL*8, DIMENSION(:), allocatable :: x_pos, y_pos, z_pos, partDia
  
! ***Program Start***  
  pi             = 3.14159265359D0
  ! particle Curtain parameters
  max_part_dia   = 130.0D-6
  delta_part_dia = 30.00D-6
  min_part_dia   = max_part_dia - delta_part_dia
  min_part_vol   = pi/6.0D0*min_part_dia**3
  part_mat       = 'St'
  
  ! Particle boundaries
  x_bmin  = -0.85D-3 
  x_bmax  =  0.85D-3
  x_range = x_bmax-x_bmin
  x_mean  = (x_bmin + x_bmax)/2.0D0
  x_sigma = 0.3D-3 ! From paper
  y_bmin  = 0.00D-3
  y_bmax  = 10.0D-3 
  y_range = y_bmax - y_bmin
  z_bmin  = 0.00D-3 
  z_bmax  = 10.0D-3 
  z_range = z_bmax - z_bmin
  total_vol = (x_range)*(y_range)*(z_range)
  
  !! Avg particle volume fraction definition and # Particles defined
  pvf_mean     = 17.0D-2
  part_num_max = pvf_mean*x_range*y_range*z_range/min_part_vol
  ALLOCATE (x_pos(1:part_num_max))
  ALLOCATE (y_pos(1:part_num_max))
  ALLOCATE (z_pos(1:part_num_max))
  ALLOCATE (partDia(1:part_num_max))

  i = 1
  part_num = 0
  pVolTotal = 0.0D0
  pVolTarget = pvf_mean*total_vol
  x_min = x_bmin 
  y_min = y_bmin + partDia(i)/2.0D0 
  z_min = z_bmin + partDia(i)/2.0D0
    
 
  main: DO 
    CALL RANDOM_NUMBER(rand1)
    partDia(i) = max_part_dia - delta_part_dia*rand1
    
    ! x has a Gaussian distribution
    CALL RANDOM_NUMBER(rand1)
    CALL RANDOM_NUMBER(rand2)
    z_dist  = SQRT(-2.0*LOG(rand1)) * SIN(2.0*pi*rand2)
    z_dist2 = SQRT(-2.0*LOG(rand1)) * COS(2.0*pi*rand2) 
    IF (MOD(i,2)==1) z_dist = z_dist2
    x_pos(i) = x_mean + x_sigma*z_dist
    
    ! y and z have pure random distribution
    CALL RANDOM_NUMBER(rand1)
    y_pos(i) = rand1*(y_range-partDia(i))+y_min
    CALL RANDOM_NUMBER(rand1)
    z_pos(i) = rand1*(z_range-partDia(i))+z_min
    !IF (x_pos(i) < (x_bmin+partDia(i)) .OR. x_pos(i) > (x_bmax-partDia(i))) CYCLE main
    checker: DO j = 1,i
      delta = SQRT((x_pos(i) - x_pos(j))**2 + (y_pos(i) - y_pos(j))**2 &
              + (z_pos(i) - z_pos(j))**2)
      IF (delta <= (1.01*(partDia(i)+partDia(j))/2) .AND. i /= j) CYCLE main
    END DO checker
    pVolTotal = pVolTotal + pi/6.0D0*partDia(i)**3
    part_num = part_num + 1
    i = i + 1
    IF (pVolTotal >= pVolTarget) THEN
      write(*,*) 'Curtain complete'
      EXIT main
    END IF
  END DO main

  ! Shift mean to x = 0.0 (or whatever value you'd like!)
  !DO i = 1,part_num
  !  x_pos(i) = x_pos(i) - 1.0D-3
  !END DO


! Printing points.dat file
! File Format:
!-----------------------------------------------------------------------------
!L1:       I16 [Number of Particles]
!L2-Npar: A16 [Part Material] 4F23.16 [x-coord] [y-coord] [z-coord] [Part Dia]
!-----------------------------------------------------------------------------

  fname    = 'points.dat'
  part_mat = 'St'
  OPEN(UNIT=11, FILE = fname, STATUS = 'NEW', ACTION = 'WRITE')
  WRITE(11,'(I16)') part_num
  DO i = 1,part_num
    WRITE(11,'(A5,11x,4(E23.16,1x))') part_mat,x_pos(i),y_pos(i),z_pos(i),partDia(i)
  END DO
  WRITE(*,*) 'points.dat file written!'
  DEALLOCATE (x_pos)
  DEALLOCATE (y_pos)
  DEALLOCATE (z_pos)
  DEALLOCATE (partDia) 
END PROGRAM WagnerPartCurtain
