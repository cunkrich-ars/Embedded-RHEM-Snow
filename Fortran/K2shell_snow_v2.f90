! Usage:
!
! k2_snow_v2 parameter-file delt gage-id cligen-file soil-class slope aspect
!
!   k2_snow_v2 is the executable file

program K2shell_snow

  USE snow
  USE, INTRINSIC::ISO_C_BINDING

  use runpars
  use elpars
  use multip
  use kinsed_def

  IMPLICIT NONE

  integer             :: i, j, k, l, nd
  integer             :: arg_len, ierr, dir_len
  character(len=2)    :: flag
  character(len=10)   :: block_name, block_names(5000)
  character(len=150)  :: arg
  character(len=150)  :: parfile
  real                :: t(5000), z(5000), sat

  REAL (C_DOUBLE), ALLOCATABLE, TARGET :: t_snow(:)
  REAL (C_DOUBLE), ALLOCATABLE, TARGET :: z_snow(:)

  INTEGER (C_LONG), TARGET :: date(3)

  REAL (C_DOUBLE)  :: sat_snow, ice_snow
  INTEGER (C_LONG) :: ret, nd_snow

  INTEGER :: yr1, yr2, iel, nel
  REAL    :: q_sum, s_sum, q_tot, s_tot, q_snow, s_snow
  REAL    :: rmks_save, rmn_save, rmg_save

! TEMP
  CHARACTER(LEN=100) :: vgage, cligen_file, soil_class, out_dir
  REAL (C_DOUBLE) :: slope, aspect
! TEMP

  q_sum = 0.
  s_sum = 0.
  q_tot = 0.
  s_tot = 0.

  first_event = .TRUE.

10 format(a) ! for general character i/o

! Default values

  cour = .false.
  sed  = .TRUE.
  tabl = .false.

  CALL GET_COMMAND_ARGUMENT( 1, parfile )

  CALL GETARG( 2, out_dir, dir_len )
  IF( out_dir(dir_len:dir_len) /= '\' )THEN
    dir_len = dir_len + 1
    out_dir(dir_len:dir_len) = '\'
  END IF

  CALL GET_COMMAND_ARGUMENT( 3, arg )
  READ( arg, * ) delt

  CALL GET_COMMAND_ARGUMENT( 4, arg, arg_len )
  vgage = arg(1:arg_len)//C_NULL_CHAR

  CALL GET_COMMAND_ARGUMENT( 5, arg, arg_len )
  cligen_file = arg(1:arg_len)//C_NULL_CHAR

  CALL GET_COMMAND_ARGUMENT( 6, arg, arg_len )
  soil_class = arg(1:arg_len)//C_NULL_CHAR

  CALL GET_COMMAND_ARGUMENT( 7, arg )
  READ( arg, * ) slope

  CALL GET_COMMAND_ARGUMENT( 8, arg )
  READ( arg, * ) aspect

  open(files(1), file = parfile, status = 'old', iostat = ierr)
   
  j = INDEX( parfile, '.' )
  OPEN( 98, FILE = out_dir(1:dir_len)//parfile(1:j)//'csv', STATUS = 'UNKNOWN' )
  OPEN( 99, FILE = out_dir(1:dir_len)//parfile(1:j-1)//'_events.csv', STATUS = 'UNKNOWN' )
  WRITE( 98, '("Year, Runoff(mm), Sediment Yield(Mg/ha)")')
  WRITE( 99, '("year, month, day, sat, rainfall volume (mm), runoff volume (mm), sediment yield (kg/ha)")' )

  delt = delt * 60. ! minutes
     
  rmks  = 1.0
  rmn   = 1.0
  rmcv  = 1.0
  rmg   = 1.0
  rmin  = 1.0
  rmcoh = 1.0
  rmspl = 1.0
  rmwco = 1.0
  rmlen = 1.0
    
  rmn_chan  = 1.0
  rmsat     = 1.0

! Read global parameter block
    
  call reader(files(1), block_name, ierr)
  if(ierr .gt. 0) call errxit('K2shell_snow', "Invalid global block")
    
  call getr4('C', 0, clen, ierr)
  if(ierr .ne. 0) call errxit('K2shell_snow', "char. length not found")
    
  units  = 1
  dlab   = 'm.'
  dilb   = 'mm'
  conv   = 1000.
  wt     = 1000.  ! kg/m3 water density
  arconv = 10000. ! m^2 per ha.
  wconv  = 1000.  ! kg/tonne
  grav   = 9.81
  bdep   = 656.

  if(sed) call sed00()

! TEMP
!  vgage = 'wy485055'//C_NULL_CHAR
!  cligen_file = 'wy485055.stm'//C_NULL_CHAR
!  soil_class = 'Loam'//C_NULL_CHAR
!  slope = 5.14
!  aspect = 270.
! TEMP

  ret = snow_run(vgage, cligen_file, soil_class, slope, aspect)

  IF( ret /= 0 )THEN
    WRITE( *, '(I6)' ) ret
    STOP ' Error from snow_run'
  END IF

! Event loop

  DO

!   Get rain/melt, sat & ice from snow model

    ret = snow_next_event(C_LOC(date))

    IF(date(1) == 0) EXIT

    nd_snow = snow_npoints()

    nd = nd_snow

    ALLOCATE(t_snow(nd), z_snow(nd))

    ret = snow_times( C_LOC(t_snow) )
    ret = snow_depths( C_LOC(z_snow) )

    t(1:nd) = t_snow * 60. ! seconds
    z(1:nd) = z_snow / 1000. ! meters

    sat_snow = snow_sat()
    ice_snow = snow_ice()

    sat = sat_snow
    ice = ice_snow

!    OPEN( 1, FILE = '0001-05-03_precip.txt', STATUS = 'OLD' )
!    READ( 1, * ) sat
!    READ( 1, * ) ice
!    READ( 1, * ) nd
!    DO j = 1, nd
!      READ( 1, * ) t(j), z(j)
!      t(j) = t(j) * 60.
!      z(j) = z(j) / 1000.
!    END DO
!    date(1) = 1
!    date(2) = 5
!    date(3) = 3

    yr2 = date(1)

    IF( first_event ) yr1 = yr2

!   Compute number of time steps

    limit = INT( t(nd) * 2. / delt ) + 1 ! t, delt in seconds

    CALL clerk(limit)

    do j = 1, 4
      dtm(j) = delt
    end do

!   Report progress

    write (*, "('+ Event ', I2.2,'/',I2.2,'/',I3.3)") date(2), date(3), date(1)
    
!   Element processing loop

    iel = 0
    iplane = 0
    ichan = 0
    ifilt = 0
    ised = 0

    DO

      IF( first_event )THEN
        call reader(files(1), block_name, ierr)
        if(ierr .gt. 0) then
          if(ierr .eq. 1)then
            nel = iel
            exit ! end of file
          else
            call errxit('K2shell_snow', "error reading parameter file")
          end if
        else
          iel = iel + 1
          block_names(iel) = block_name
        end if
      ELSE
        iel = iel + 1
        IF( iel > nel) EXIT
      END IF
     
      nchan  = 1
    
!     Elements

      if(block_names(iel)(1:5) .eq. 'PLANE') then
        iplane = iplane + 1
        call plane(t, z, nd, sat)
      else if(block_names(iel)(1:7) .eq. 'CHANNEL') then
        call channl(t, z, nd, sat)
      else if(block_names(iel)(1:4) .eq. 'PIPE') then
        call pipe()                      
      else if(block_names(iel)(1:6) .eq. 'INJECT') then
        call inject()
      else if(block_names(iel)(1:4) .eq. 'POND') then
        call pond(t, z, nd)
      else if(block_names(iel)(1:5) .eq. 'URBAN') then                       
        call urban(t, z, nd, sat)
      else if(block_names(iel)(1:5) .eq. 'ADDER') then                       
        call adder()
      else if(block_names(iel)(1:8) .eq. 'DIVERTER') then                       
        call diverter()
      else
        call errxit(block_names(iel), 'invalid element')
      end if
      
! End of element loop

    END DO

!   Write year, month, day, sat, rainfall volume (mm), runoff volume (mm), sediment yield (kg/ha)

    q_snow = qbal(10)/qbal(1)*conv
    s_snow = twso_snow*wt/qbal(1)*arconv

    IF( q_snow > 0. )THEN
      WRITE(99,'(I4.4,",",I2.2,",",I2.2,",",F12.4,",",F12.4,",",F12.4,",",F12.4)') yr2, date(2), date(3), satin, qbal(2)/qbal(1)*conv, q_snow, s_snow
    END IF

    IF( yr2 == yr1 )THEN
      q_sum = q_sum + q_snow
      s_sum = s_sum + s_snow
    ELSE
      WRITE( 98, '(I4,",",F6.1,",",F8.3)') yr1, q_sum, s_sum/1000. ! mm, Mg/ha
      q_tot = q_tot + q_sum
      s_tot = s_tot + s_sum
      q_sum = q_snow
      s_sum = s_snow
      yr1 = yr2
    END IF

    DEALLOCATE(t_snow, z_snow)
    
    first_event = .FALSE.

! End of event loop

  END DO

  WRITE( 98, '(" Avg,",F6.1,",",F8.3)') q_tot/300., s_tot/1000./300.

  CLOSE( 98 )
  CLOSE( 99 )

  CLOSE( files(1) )

  ret = snow_finale()

  IF( ret /= 0 )THEN
    STOP ' Error from snow_finale'
  END IF

end program K2shell_snow


logical function char2log(char)

! Converts character values Y|N or y|n into true|false
! Returns char as uppercase

  character(LEN=1) :: char

  if(char .eq. 'Y' .or. char .eq. 'y') then
    char2log = .true.
    char = 'Y'
  else
    char2log = .false.
    char = 'N'
  end if
end

subroutine errxit (id, message)

! Write error message to stdout & output file

   use runpars

   character(len=*) id, message

   write(files(3), 10) id, trim(message)
   write(*, 10)        id, trim(message)
10 format(//, ' error - ', a, /, 1x, a)
   stop ' '
end
