! Run k2_snow from a command shell
!
! Usage:
!
! k2_snow -s file
!
!   k2_snow is the K2shell_snow executable file
!   -s suppresses screen output (optional)
!   file is the run file (default = kin.fil)
!
! The run file has the following comma-separated values:
!
! parameter-file, output-file, title, tfin, delt, cour?, sed?, multiplier-file?, table?
!
! The title must be enclosed in double quotes.
!
! Entries after delt are optional, will default to 'N'

! Citation:
! Goodrich, D.C., I.S. Burns, C.L. Unkrich, D.J. Semmens, D.P. Guertin, M. Hernandez, S. Yatheendradas, J.R. Kennedy,
! and L.R. Levick. 2012. K2/AGWA: model use, calibration, and validation. Transactions of the ASABE 55(4):1561-1574.
! 
! Disclaimer:
! No warranties, expressed or implied, are made that this program is free from errors or will meet the requirements
! for any particular application.  The U.S. Department of Agriculture disclaims all liability for direct or
! consequential damages resulting from the use of this program.

program K2shell_snow

  USE snow
  USE, INTRINSIC::ISO_C_BINDING

  use runpars
  use elpars
  use multip
  use kinsed_def

  IMPLICIT NONE

  integer             :: i, j, k, l, nd
  integer             :: arg_len, ierr
  integer             :: tvals(8), t0, te
  logical             :: char2log, extra_mults
  character(len=1)    :: bell, c
  character(len=2)    :: flag
  character(len=10)   :: block_name, block_names(5000)
  character(len=11)   :: version
  character(len=150)  :: runfile, fields(20)
  character(len=1000) :: line
  real                :: t(5000), z(5000), sat
  real                :: rmks_chan
  real                :: rmg_chan

  REAL (C_DOUBLE), ALLOCATABLE, TARGET :: t_snow(:)
  REAL (C_DOUBLE), ALLOCATABLE, TARGET :: z_snow(:)

  INTEGER (C_LONG), TARGET :: date(3)

  REAL (C_DOUBLE)  :: sat_snow, ice_snow
  INTEGER (C_LONG) :: ret, nd_snow

  INTEGER :: yr1, yr2, iel, nel
  REAL    :: q_sum, s_sum, q_tot, s_tot, q_snow, s_snow
  REAL    :: rmks_save, rmn_save, rmg_save

! TEMP
  CHARACTER(LEN=100) :: vgage, cligen_file, soil_class
  REAL (C_DOUBLE) :: slope, aspect
! TEMP

  q_sum = 0.
  s_sum = 0.
  q_tot = 0.
  s_tot = 0.

  first_event = .TRUE.

10 format(a) ! for general character i/o

  version = '11-Mar-2023'
  bell    = char(7)
  silent  = .FALSE.

! Default values

  cour = .false.
  sed  = .false.
  tabl = .false.

  if(COMMAND_ARGUMENT_COUNT() == 2) then
    call GET_COMMAND_ARGUMENT( 1, flag, arg_len )
    if(flag == '-s') silent = .TRUE.
    call GET_COMMAND_ARGUMENT( 2, runfile, arg_len )
  else if(COMMAND_ARGUMENT_COUNT() == 1) then
    call GET_COMMAND_ARGUMENT( 1, runfile, arg_len )
    if( arg_len == 2 )then
      silent = .TRUE.
      runfile = 'kin.fil'
      arg_len = 7
    end if
  else
    runfile = 'kin.fil'
    arg_len = 7
  end if

  open(4, file = runfile(1:arg_len), status = 'old', iostat = ierr)
  if(ierr .ne. 0) call errxit('Kineros2', "Can't open "//runfile(1:arg_len))

  CALL DATE_AND_TIME( values = tvals )
  t0 = tvals(5) * 3600 + tvals(6) * 60 + tvals(7)

  read(4, 10, iostat = ierr) line
   
! Replace double quotes around title field with spaces
   
  i = index(line, '"')
  line(i:i) = ' '
  j = index(line, '"')
  line(j:j) = ' '
   
! Temporarily replace any commas in title field with bell characters
     
  k = index(line(i:j), ',')
  do while( k > 0)
    k = i + k - 1
    line(k:k) = bell
    k = index(line(i:j), ',')
  end do
   
! Identify comma-separated fields
     
  fields = ' '
  i = 1
  j = 1
  k = index(line, ',')
  do while( k > 0)
    line(k:k) = ' '
    fields(i) = line(j:k-1)
    i = i + 1
    j = k + 1
    k = index(line, ',')
  end do
  fields(i) = line(j:)
   
! Restore commas to title field
     
  k = index(fields(3), bell)
  do while( k > 0)
    fields(3)(k:k) = ','
    k = index(fields(3), bell)
  end do

! Parameter file   

  open(files(1), file = trim(fields(1)), status = 'old', iostat = ierr)
  if(ierr .ne. 0) call errxit('K2shell', "Can't open "//trim(fields(1)))
   
! Output files

  OPEN( 98, FILE = fields(2), STATUS = 'UNKNOWN' )
  j = INDEX( fields(2), '.' )
  k = LEN_TRIM( fields(2) )
  OPEN( 99, FILE = fields(2)(1:j-1)//'_events'//fields(2)(j:k), STATUS = 'UNKNOWN' )
  WRITE( 98, '("Year, Runoff(mm), Sediment Yield(Mg/ha)")')
  WRITE( 99, '("year, month, day, sat, rainfall volume (mm), runoff volume (mm), sediment yield (kg/ha)")' )

! Decode character input into floating and logical values (module runpars)
   
  read(fields(4) , *) tfin
  read(fields(5) , *) delt
  
  delt = delt * 60. ! minutes
     
! If no value for a logical variable is specified, it remains at its default value
     
  if(len_trim(fields(6))  > 0) cour = char2log(fields(6))
  if(len_trim(fields(7))  > 0) sed  = char2log(fields(7))
  if(len_trim(fields(9)) > 0) tabl = char2log(fields(9))
   
! Two options for parameter multipliers:
   
! 1) Use file specified in run file
! 2) No multipliers 
   
  if(len_trim(fields(8)) == 0 .or. verify(fields(8),'N ') == 0 .or. verify(fields(8),'n ') == 0) then ! no multipliers
   
    rmks_save = 1.0
    rmn_save  = 1.0
    rmg_save  = 1.0

    rmcv  = 1.0
    rmin  = 1.0
    rmcoh = 1.0
    rmspl = 1.0
    
    rmks_chan = 1.0
    rmg_chan  = 1.0
    rmn_chan  = 1.0
    rmwco     = 1.0
    rmlen     = 1.0
    rmsat     = 1.0

  else
    
    open(7, file = fields(8), status = 'old', iostat = ierr)
    if(ierr .ne. 0) call errxit('K2shell_snow', "Can't open "//fields(9))
    
    read(7, *) rmks_save
    read(7, *) rmn_save
    read(7, *) rmcv
    read(7, *) rmg_save
    read(7, *) rmin
    read(7, *) rmcoh
    read(7, *) rmspl

    read(7, *, iostat = ierr) rmks_chan

    if(ierr .ne. 0)THEN
      rmks_chan = 1.0
      extra_mults = .FALSE.
    else
      extra_mults = .TRUE.
    end if

    read(7, 10, iostat = ierr) rmg_chan
    if(ierr .ne. 0) rmg_chan = 1.0

    read(7, 10, iostat = ierr) rmn_chan
    if(ierr .ne. 0) rmn_chan = 1.0

    read(7, 10, iostat = ierr) rmwco
    if(ierr .ne. 0) rmwco = 1.0

    read(7, 10, iostat = ierr) rmlen
    if(ierr .ne. 0) rmlen = 1.0

    read(7, 10, iostat = ierr) rmsat
    if(ierr .ne. 0) rmsat = 1.0

    close(7)
    
  end if
    
! Read global parameter block
    
  call reader(files(1), block_name, ierr)
  if(ierr .gt. 0) call errxit('K2shell_snow', "Invalid global block")
    
  call getr4('C', 0, clen, ierr)
  if(ierr .ne. 0) call errxit('K2shell_snow', "char. length not found")
    
  call getstr ('U', 0, c, l, ierr)
  if(ierr .ne. 0) call errxit('K2shell_snow', "units not specified")
    
  if(c .eq. 'm' .or. c .eq. 'M') then ! metric
    units  = 1
    dlab   = 'm.'
    dilb   = 'mm'
    conv   = 1000.
    wt     = 1000.  ! kg/m3 water density
    arconv = 10000. ! m^2 per ha.
    wconv  = 1000.  ! kg/tonne
    grav   = 9.81
    bdep   = 656.
  else if(c .eq. 'e' .or. c .eq. 'E') then ! english
    units  = 2
    dlab   = 'ft'
    dilb   = 'in'
    conv   = 12.
    wt     = 62.4    
    arconv = 43560.
    wconv  = 2000.
    grav   = 32.2
    bdep   = 200.
  else
    call errxit('K2shell_snow', "units not recognized")
  end if

  if(sed) call sed00()

! TEMP
  vgage = 'wy485055'//C_NULL_CHAR
  cligen_file = 'wy485055.stm'//C_NULL_CHAR
  soil_class = 'Loam'//C_NULL_CHAR
  slope = 5.14
  aspect = 270.
! TEMP

  ret = snow_run(vgage, cligen_file, soil_class, slope, aspect)

  IF( ret /= 0 )THEN
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

    yr2 = date(1)

    IF( first_event ) yr1 = yr2

!   Report progress

    if(.not. silent) write (*, "('+ Event ', I2.2,'/',I2.2,'/',I3.3)") date(2), date(3), date(1)
    
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
    
!     Compute number of time steps

      limit = INT( t(nd) * 2. / delt ) + 1 ! t, delt in seconds
      CALL clerk(limit)

      do j = 1, 4
        dtm(j) = delt
      end do

!     Multipliers

      IF( extra_mults .AND. block_names(iel)(1:7) == 'CHANNEL' )THEN
        rmks = rmks_chan
        rmg  = rmg_chan
        rmn  = rmn_chan
      ELSE
        rmks = rmks_save
        rmg  = rmg_save
        rmn  = rmn_save
      END IF

!     Elements

      if(block_names(iel)(1:5) .eq. 'PLANE') then
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

  CALL DATE_AND_TIME( values = tvals )
  te = tvals(5) * 3600 + tvals(6) * 60 + tvals(7) - t0
  write (*, "(' Elapsed time (sec): ', i6.6)") te ! / 60

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
