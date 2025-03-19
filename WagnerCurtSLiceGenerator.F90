PROGRAM WagnerPartCurtain
IMPLICIT NONE

INTEGER*4          :: part_num_max, part_num, i, j, k, nslices 
CHARACTER(len=200) :: fname        
CHARACTER(len=5)   :: part_mat !If len changes, change A(len) in write statement format      
REAL*8 :: pvf, part_dia, part_vol, x_min, x_max, y_min, y_max, x_center, &
          rand, z_min, z_max, x_range, y_range, z_range, total_vol, pi, delta, &
          x_bmin, x_bmax, dx, pvf_max, pvf_min, pvf_dist, max_part_dia, delta_part_dia, &
          pVolTotal, y_bmin, y_bmax, z_bmin, z_bmax, x_minSlice, pVolTarget, &
          min_part_vol, min_part_dia
REAL*8, DIMENSION(:), allocatable :: x_pos, y_pos, z_pos, partDia
REAL*8 :: x_mean, x_std, u1, u2, z1 !for normal pvf distribution
  pi             = 3.14159265359D0
  ! particle Curtain parameters
  max_part_dia   = 130.0D-6
  delta_part_dia = 30.00D-6
  min_part_dia   = max_part_dia - delta_part_dia
  min_part_vol   = pi/6.0D0*min_part_dia**3
  part_mat       = 'St'
  
  ! varying vf dimension
  x_bmin  = 0.00D-3 
  x_bmax  = 1.70D-3
  x_range = x_bmax - x_bmin
  nslices = 7
  dx      = x_range/nslices !slice length in x for varying vf

  ! Wagner full test window is 68.3mm x 68.3mm
  ! Assuming we're only looking at a section of the curtain
  ! defining constant vf dimensions
  y_bmin  = 0.00D-3
  y_bmax  = 10.0D-3 
  y_range = y_bmax - y_bmin
  z_bmin  = 0.00D-3 
  z_bmax  = 10.0D-3 
  z_range = z_bmax - z_bmin
  
  !! linear particle volume fraction distribution
  pvf_min        = 0.10D0
  pvf_max        = 0.22D0  
  part_num_max   = pvf_max*x_range*y_range*z_range/min_part_vol
  ALLOCATE (x_pos(1:part_num_max))
  ALLOCATE (y_pos(1:part_num_max))
  ALLOCATE (z_pos(1:part_num_max))
  ALLOCATE (partDia(1:part_num_max))

  i = 1
  part_num = 0
  WRITE(*,*) 'Number of varying vf slices: ', nslices
  DO k = 1,nslices+1
    ! set x-coordinate for current slice iteration
    x_minSlice = x_bmin + (k-1)*dx
    IF (x_minSlice > x_bmax) exit
    x_max   = x_bmin + (k)*dx
    IF (x_max > x_bmax) x_max = x_bmax
    x_range = x_max - x_minSlice
    x_center = (x_max+x_minSlice)/2
    total_vol = (x_range)*(y_range)*(z_range)
 
    ! set pvf for current slice - linear
    IF (x_center <  (x_bmax-x_bmin)/2) THEN 
      pvf = pvf_min + (pvf_max-pvf_min)*(x_center - x_bmin)/((x_bmax-x_bmin)/2)
    ELSE
      pvf = pvf_min + (pvf_max-pvf_min)*(x_bmax - x_center)/((x_bmax-x_bmin)/2)
    END IF
    pVolTotal = 0.0D0
    pVolTarget = pvf*total_vol
    
    main: DO 
      CALL RANDOM_NUMBER(rand)
      partDia(i) = max_part_dia - delta_part_dia*rand
      x_min = x_minSlice + partDia(i)/2.0D0 
      y_min = y_bmin + partDia(i)/2.0D0 
      z_min = z_bmin + partDia(i)/2.0D0
      CALL RANDOM_NUMBER(rand)
      x_pos(i) = rand*(x_range-partDia(i))+x_min
      CALL RANDOM_NUMBER(rand)
      y_pos(i) = rand*(y_range-partDia(i))+y_min
      CALL RANDOM_NUMBER(rand)
      z_pos(i) = rand*(z_range-partDia(i))+z_min
      checker: DO j = 1,i
        delta = SQRT((x_pos(i) - x_pos(j))**2 + (y_pos(i) - y_pos(j))**2 &
                + (z_pos(i) - z_pos(j))**2)
        IF (delta <= (1.01*(partDia(i)+partDia(j))/2) .AND. i /= j) CYCLE main
      END DO checker
      pVolTotal = pVolTotal + pi/6.0D0*partDia(i)**3
      part_num = part_num + 1
      i = i + 1
      IF (pVolTotal >= pVolTarget) THEN
        write(*,*) 'slice complete'
        EXIT main
      END IF
    END DO main
    WRITE(*,*) 'Current # of Particles Inserted: ',part_num
    WRITE(*,*) 'x-coord: ', x_center, 'VF:', pvf
  END DO !k

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
