! Array Visualizer code added...search on 'fagl'

C   Code type: FORTRAN subroutine

C   Compiler: Fortran77

C   Programmed by: C. Unkrich

C   Date: 4/95

C   Revised: 7/95 mods for microtopography by R.E. Smith
C!! slight revision in calculating depth to send to infilt routine 2/12/02

C    Laminar flow added 8/2005 by C. Unkrich - laminar flow is used up to a
C    transition depth if a 'LAM' tag & value appear in the input block

C   Description:

C     Simulates overland flow using a four point implicit finite-difference
C     approximation to the kinematic wave equation written for one-dimensional
C     unsteady flow on a sloping plane. Rainfall is modeled as a spatially
C     uniform time, intensity step function. The infiltration function is
C     coupled to the routing equations to provide a realistic interaction
C     between rainfall, runoff, and infiltration. Modifed for micro-topograpy,
C     using hydraulic radius rather than mean depth. An algorithm has been
C     added which monitors the wave celerity and adjusts the computational
C     time increment to maintain numerical accuracy.

C     For output, volumetric rainfall rate from the upstream area for each
C     time step is added to the volumetric rate over the plane and stored.
C     Dividing by contributing area gives the integrated rainfall rate over
C     the area represented by the plane and its upstream contributors.
C------------------------------------------------------------------------------

C   Arguments:

C     dtm(4)        real        smallest dt to satisfy the Courant condition
C                               during each quarter of the simulation.
C
C     xnu           real        kinematic viscosity (m**2/s or ft**2/s)
C -----------------------------------------------------------------------------
C  for other variables see module definition file
C------------------------------------------------------------------------------

C   Input Block (input file labels in parenthesis):

C     id (ID)          int      identifier,

C     xleng (L)        real     length (m or ft),

C     width (W)        real     width (m or ft),

C     slope (SL)       real     slope,

C     rs (MA or CH)    real     Manning or Chezy discharge coefficient,

C     x (X)            real     x coordinate,

C     y (Y)            real     y coordinate,

C     dintr (IN)       real     interception depth (mm or in),

C     cover (CA)       real     canopy cover fraction,

C     amr (RE)         real     average microtopographic relief (mm or in),

C     ams (SPA)        real     average microtopographic spacing (m or ft),

C     fmin (KS)        real     saturated hydraulic conductivity (mm or in/hr),
C------------------------------------------------------------------------------

C   Subroutines/functions:

C     interp        (rain.for)
C     infilt        (infilt.for)
C     infil0             "
C     iter          (miscel.for)
C     sed0          (kinsed.for)
C     kinsed             "
C     sedfin             "
C     errh2         (plane.for)
C     erpaub             "
C     qfunc              "
C     wprf               "
C     zfu                "
C     dqdaf              "
C     progss        (K2RunW.for)
C     getr4         (reader.for)
C     getstr             "
C     new           (miscel.for)
C     old                "
C     get                "
C     store              "
C     stora              "
C     geta               "
C     qwrt          (writer.for)
C     errxit        (K2RunW.for)
C------------------------------------------------------------------------------


      subroutine plane (t, z, nd, sat)

      use runpars
      use elpars
      use multip
      use itpars
      use kinsed_def
      
C**      USE AVDEF
C                                                                     arguments
C------------------------------------------------------------------------------

!      logical sed, diag, cour

!      integer limit, units
!      dimension qbal(10)

C                                                               local variables
C------------------------------------------------------------------------------

      logical :: up, unif, laminar, skip

      integer :: nd, ierr, nk, j, jp1, ms, upq, outq(2), i, k, idup,
     &           inflow(12), nu, upr, outr, ngage, nzro=0

      character(LEN=20) :: msg
      character(LEN=6) :: source
C** new units character  22/4/03
      character(LEN=2) :: chl

      common /plane1/ dt, dx, j, jp1, qlav, alpha, p, h1(20), h2(20), 
     &                q1(20), q2(20), h2c, ams, wpc, c1, c2, bw, qb,
     &                omega
C
      real,dimension(10,2) :: qup
      real,dimension(20,2) :: ql
      real,dimension(20) ::  topw, slp, fav, qt, afac, h1avg, xc, qz
      real,dimension(5000) :: t, z, r, zn
      real, parameter :: dtmin = 0.
      real, parameter :: zro = 0.
      real :: qtrans
C                                               slopes and horizontal distances
      real,dimension(20)   :: slope, dist
      integer              :: npairs

      CHARACTER(LEN=20), SAVE :: msg_array(1000)

      INTEGER, SAVE :: id_array(1000)
      INTEGER, SAVE :: nk_array(1000)
      INTEGER, SAVE :: idup_array(1000)
      INTEGER, SAVE :: nu_array(1000)
      INTEGER, SAVE :: ms_array(1000)

      LOGICAL, SAVE :: up_array(1000)

      REAL, SAVE :: r_ks_array(1000)
      REAL, SAVE :: r_cv_array(1000)
      REAL, SAVE :: xleng_array(1000)
      REAL, SAVE :: dx_array(1000)
      REAL, SAVE :: width_array(1000)
      REAL, SAVE :: area_array(1000)
      REAL, SAVE :: slp_array(20,1000)
      REAL, SAVE :: p_array(1000)
      REAL, SAVE :: alpha_array(1000)
      REAL, SAVE :: dintr_array(1000)
      REAL, SAVE :: cover_array(1000)

C------------------------------------------------------------------------------

      external errh2

      external erpaub

      data topw /20*1./
      data qt /20*0./
      data qup /20*0./

      omega = 0.8        !  0.8 to match chan
      tol = 1.e-7
      ltyp = 0
      pof = conv*3600.

      if (units .eq. 1) then
        h_lim = 2. ! meters
      else
        h_lim = 6. ! ft
      end if

C TEMP SNOW
      IF( first_event )THEN

        call getr4 ('KE', 1, r_ks, ierr)
        if (ierr .gt. 0) r_ks = 0.

        r_ks = r_ks / (3600. * conv)

        call getr4 ('CV', 1, r_cv, ierr)
        if (ierr .gt. 0) r_cv = 0.

        r_ks_array(iplane) = r_ks
        r_cv_array(iplane) = r_cv

C TEMP SNOW

C                                                          get input parameters
C------------------------------------------------------------------------------

C                                                                    IDENTIFIER
        call geti4 ('ID', 0, id, ierr)

        if (ierr .eq. 3) then
          call errxit ('PLANE', 'invalid identifier (ID)')
        else if (ierr .ge. 1) then
          call errxit ('PLANE', 'missing identifier (ID)')
        end if

        msg(1:16) = typname(0)
        write (msg(15:20), '(i6)') id
C
        do ms = 1, 5
          if (msg(ms:ms) .ne. ' ') go to 1
        end do
1       continue

        id_array(iplane) = id
        msg_array(iplane) = msg
        ms_array(iplane) = ms

C                                                                        LENGTH
        call getr4 ('L', 0, xleng, ierr)

        if (ierr .gt. 0) call errxit(msg(ms:14), 'length (L) not found')

        nk = max1 ((15. * xleng / clen), 5.)
C                                                       number of spatial nodes
        if (nk .gt. 15) nk = 15
C                                                            distance increment
        dx = xleng / (float(nk) - 1.)
C                                                                         WIDTH
        call getr4 ('W', 0, width, ierr)

        if (ierr .gt. 0) then
C                                                          area (sq.m or sq.ft)
          call getr4 ('AR', 0, area, ierr)

          if (ierr .gt. 0) then

            if (units .eq. 1) then
C                                                               area (hectares)
              call getr4 ('HA', 0, area, ierr)

              if (ierr .gt. 0) call errxit
     &                       (msg(ms:14), 'width (W) or area not found')

              area = area * 10000.

            else
C                                                                  area (acres)
              call getr4 ('AC', 0, area, ierr)

              if (ierr .gt. 0) call errxit
     &                       (msg(ms:14), 'width (W) or area not found')

              area = area * 43560.

            end if
          end if

          width = area / xleng

        else

          area = xleng * width

        end if

C TEMP SNOW
        xleng_array(iplane) = xleng
        nk_array(iplane) = nk
        dx_array(iplane) = dx
        width_array(iplane) = width
        area_array(iplane) = area
C TEMP SNOW
C                                                                         SLOPE
        call getr4 ('SL', 1, slope(1), ierr)

        if (ierr .gt. 0) call errxit(msg(ms:14), 'slope (SL) not found')
        if (slope(1) .le. 0.) call errxit (msg(ms:14), 'zero slope')

        dist(1) = 0.
        npairs = 1

        do i = 2, 20
          np1 = npairs + 1
          call getr4 ('SL', np1, slope(np1), ierr)
          if (ierr .gt. 0) go to 2
          if (slope(np1) .le. 0.) call errxit (msg(ms:14), 'zero slope')
          call getr4 ('SX', np1, dist(np1), ierr)
          if (ierr.gt.0) call errxit(msg(ms:14),'distance SX not found')
          dist(np1) = dist(np1) * xleng
          npairs = np1
        end do

2       continue

        if (npairs .gt. 1) then

          slp(1) = slope(1)
          a = (slope(2) - slope(1)) / (dist(2) - dist(1))
          k = 2

          do j = 2, nk - 1
            x = float(j-1) * dx
3           if (x .gt. dist(k) .and. k .lt. npairs) then
              k = k + 1
              a = (slope(k) - slope(k-1)) / (dist(k) - dist(k-1))
              go to 3
            end if
            slp(j) = a * (x - dist(k-1)) + slope(k-1)
          end do

          slp(nk) = slope(npairs)

        else

          do j = 1, nk
            slp(j) = slope(1)
          end do

        end if

C TEMP SNOW
        do j = 1, nk
          slp_array(j,iplane) = slp(j)
        end do
C TEMP SNOW

C                                                                      MANNING?
        call getr4 ('MA', 0, rs, ierr)

        if (ierr .eq. 0) then

          p = 1.66667
          fac = 1.
          if (units .eq. 2) fac = 1.49

          do j = 1, nk
            alpha = fac * sqrt (slp(j)) / (rs * rmn)
          end do

        else
C                                                                        CHEZY?
          call getr4 ('CH', 0, rs, ierr)

          if (ierr .eq. 0) then

            p = 1.5

            do j = 1, nk
              alpha = rs * rmn * sqrt (slp(j))
            end do

          else

            call errxit
     &      (msg(ms:14), 'hydraulic roughness (MA or CH) not found')

          end if
        end if

C TEMP SNOW
        p_array(iplane) = p
        alpha_array(iplane) = alpha
C TEMP SNOW
C                                                                  INTERCEPTION
        call getr4 ('IN', 0, dintr, ierr)

        if (ierr .gt. 0) dintr = 0.

C TEMP SNOW
        dintr_array(iplane) = dintr
C TEMP SNOW

        if (dintr .gt. 0.) then

          dintr = rmin * dintr
C                                                                  CANOPY COVER
          call getr4 ('CA', 0, cover, ierr)
C!! change .eq. to .le.
          if (ierr .gt. 0 .or. cover .le. 0.)
     &    call errxit(msg(ms:14),'canopy cover fraction (CA) not found')
C!! check for out of bounds:
          if(cover .gt. 1.0) then
            call errxit (msg(ms:14),' canopy cover cannot be > 1')
          end if

C TEMP SNOW
          cover_array(iplane) = cover
C TEMP SNOW

        end if
C                                                                        RELIEF
        amr = 0.
        ams = 0.

C      call getr4 ('RE', 0, amr, ierr)

C      if (ierr .eq. 0) then
C                                                                       SPACING
C        call getr4 ('SPA', 0, ams, ierr)

C        if (ierr .gt. 0) call errxit
C     &  (msg(ms:14), 'average microtopographic spacing (SPA) not found')

C        amr = amr / conv

C      end if

        up_array(iplane) = .false.
        nu_array(iplane) = 0
        upr = 0
C                                                     upstream plane identifier
C TEMP SNOW
        call geti4 ('UP', 0, idup_array(iplane), ierr)
C TEMP SNOW

        if (ierr .eq. 0) then
          up_array(iplane) = .true.
          nu_array(iplane) = 1
        end if

C TEMP SNOW
      END IF

      r_ks = r_ks_array(iplane)
      r_cv = r_cv_array(iplane)
      id = id_array(iplane)
      msg = msg_array(iplane)
      ms = ms_array(iplane)
      xleng = xleng_array(iplane)
      nk = nk_array(iplane)
      dx = dx_array(iplane)
      width = width_array(iplane)
      area = area_array(iplane)
      slp = slp_array(:,iplane)
      p = p_array(iplane)
      alpha = alpha_array(iplane)
      dintr = dintr_array(iplane)
      cover = cover_array(iplane)
      idup = idup_array(iplane)
      up = up_array(iplane)
      nu = nu_array(iplane)
C TEMP SNOW
C                                              get indices to storage locations
C                                              for rainfall, inflow and outflow
C------------------------------------------------------------------------------
      call new (id, 'RA', outr)

      if (up) then

        call old (idup, 'RA', upr, nr, ierr)
        call old (idup, 'RO', upq, nq, ierr)
        if (ierr .gt. 0) call errxit
     &  (msg(ms:14), 'upstream inflow not found')
        inflow(1) = idup
        
      else
        upr = 0
        upq = 0
      end if

      call new (id, 'RO', outq(1))
C                                                      update contributing area
C------------------------------------------------------------------------------
      coarea = 0.
      if (up) call geta (idup, coarea)
      coarea = coarea + area
      call stora (id, coarea)
C                                                                      rainfall
C------------------------------------------------------------------------------

      d1 = 0.
      t1 = 0.
      di = 0.
      dinact = 0.
      dintr = dintr / conv   ! convert to m.
      tfinl = float (limit - 1) * delt   !this is in seconds

C TEMP SNOW
      if (t(nd) .lt. tfinl) then
        t(nd+1) = tfinl + 1.
        z(nd+1) = z(nd)
        nd = nd + 1
      end if
      skip = .true.
C TEMP SNOW

      itlim = limit
      i = 0
      ndnew = nd
C!!
      zn(1) = 0.
      depth = z(nd)
      z(nd+1) = z(nd)

      do j = 1, nd

        i = i + 1

        if (i .eq. ndnew) go to 10
C                                              convert depth z into intensity r
        t2 = t(i+1)
        dt = t2 - t1
        d2 = z(i+1)
        dd = d2 - d1

        if (dintr .gt. d1) then
C                                            reduce intensity by cover fraction
          if (d2 .gt. dintr) then
C                                               reduction over partial interval
            r(i) = dd * (1. - cover) / dt
C!!  intermediate time
            adt = dt*(dintr - d1)/dd
C                                                             insert breakpoint
            r(i+1) = dd / dt
C                                                reset rain array indices
            do k = nd, i+1, -1
              t(k+1) = t(k)
              z(k+1) = z(k)
            end do
C!!  new 10/03  find array zn which is the z array less interception
            z(i+1) = z(i) + dd*adt/dt
C            
            zn(i+1) = zn(i) + r(i)*adt
            t(i+1) = t(i) + adt 
            i = i + 1
            zn(i+1) = zn(i) + r(i)*(dt-adt)
            ndnew = nd + 1
            di = dintr

          else

            r(i) = dd * (1. - cover) / dt
            di = d2
C!!  new 10/03
            zn(i+1) = zn(i) + r(i)*dt

          end if

        else
C                                                               no interception
          r(i) = dd / dt
          zn(i+1) = zn(i) + r(i)*dt    !! also new 10/03

        end if

C TEMP SNOW
      r_ = r(i)/r_ks
      if (r_ .lt. 0.2) then
         vfm = 1.
      else if (r_cv .lt. 0.1) then
         vfm = 1.
      else
        p_ = 1.8 / r_cv**.85
        vfm = (1. + (1. / r_)**p_)**(-1. / p_)
      end if
      if (r(i) .gt. r_ks * vfm ) skip = .false.
C TEMP SNOW

        if (tfinl .gt. t1 .and. tfinl .le. t2) then
C                                                  compute final rainfall depth
          dtfin = (t2 - tfinl) / dt
          depth = d2 - dtfin * dd

          if (dintr .gt. 0.) dinact = di * cover

        end if

        t1 = t2
        d1 = d2

      end do

10    continue

      if(ndnew .lt. 5000) then
        t(ndnew+1) = tfinl + 1.
        r(ndnew+1) = 0.
      end if
      
C                                                       final depth intercepted
      dintr = dinact
C                                                 compute integrated rain rates
C------------------------------------------------------------------------------
      qrmax = 0.
      k = 1
      qrf = r(1) * area
C!!  new 9/03
      lmn = limit - 1
      zlr = zn(1)
C!!
      do i = 1, lmn

C!! replaced 9/03:
C        if (time .gt. t(k+1)) then
C          k = k + 1
C          if(k .gt. ndnew) then   
C            stop ' rainfall undefined at long time '
CC            r(k) = 0.
CC            ndnew = k
C          end if
C          qrf = r(k) * area
C        end if
C!! new code for qrf to account for rainsteps less than DELT:
        timen = delt * float (i)   ! seconds
        do while(t(k) .lt. timen)
          k = k+1
        end do
C interpolate cumulative depth at time: 
C!!  use zn rather than z to account for interception  (10/03)
        km = k-1
        z2 = zn(km) + (zn(k) - zn(km))*(timen - t(km))/(t(k) - t(km))
        qrf = (z2 - zlr)/DELT * area
        zlr = z2
!! end new code
        qrain = qrf

        if (upr .gt. 0) then
          call get (upr, i, qr)
          qrain = qrain + qr
        end if

        call store (outr, i, qrain)
        if (qrain .gt. qrmax) qrmax = qrain

      end do
C!!  add value for i=limit (maybe unnecessary)
      call store (outr, limit, zro)

      do j = 1, 15
        qbal(j) = 0.
      end do
C                                                 plane are in m2 or ft2
      qbal(1) = area                  
C                                                               rainfall volume
      qbal(2) = area * depth
C                                                            volume intercepted
      qbal(8) = area * dintr

C                                     read / initiallize infiltration variables
C------------------------------------------------------------------------------

      call infil0 (id, 1, nk, msg(ms:14), diag, sat)

      if (sed) then
C                                                              ... and sediment
C------------------------------------------------------------------------------

        call sed0 ( msg(ms:14), inflow, nu, nzro,
     &              nk, omega, dx, slp, width)

      end if

C                        compute geometric variables related to microtopography
C------------------------------------------------------------------------------

      bf = 0.1

      if (amr .gt. 0.) then         ! relief height = amr
C
        bm = .5 * (1. - bf) * ams
        ss = amr / bm
        cpz = sqrt(1. + ss * ss)
        bw = bf * ams
        ac = amr * (bm + bw)
        wpc = bw + cpz * 2. * bm ! critical wetted perim for one rill
        rmc = ac / wpc
        h2c = ac / ams
        uqc = alpha * rmc**(p - 1.) * ac
        wz = bf
        fnw = ams / width
        c1 = 2. / ss
        c2 = c1 * cpz

      else

        h2c = 0.
        uqc = 0.
        bw = 1.
        wpc = 1.
        ams = 1.
        wz = 1.

      end if

C                                                   initiallize local variables
C------------------------------------------------------------------------------
C                                                                    Q=0 at t=0
      do j = 1, nk
        h1(j) = 0.
        h2(j) = 0.
        q1(j) = 0.
        q2(j) = 0.
        ql(j,1) = 0.
        fav(j) = r(1)
        topw(j) = wz
        xc(j) = 0.
      end do

      time = 0.
      xc_max = 0.
      k = 1
      vin = 0.
      unif = .true.
      nc = 1
      sumrf = 0.

      p1 = limit / 4
      p2 = limit / 2
      p3 = 3 * limit / 4
      qub1 = 0.                           ! default
C!!  new  RES      
      if (up) call get (upq, 1, qub1)
      vin = 0.5 * qub1
      q1s = qub1 / width
C estimate first a2u
      a2u = (q1s * ams / alpha * wpc**(p - 1.))**(1. / p)
C      h1(1) = a2u/ams
C!! end new      
C!! new location so that upstream record is not written over
      call store (outq(1), 1, zro)
      qout1 = 0.
      qp = 0.
      tp = 0.
      qsp = 0.
      tsp = 0.
      wsi = 0.
      wsmax = 0.
      vcm = 0.
      afac = 1.0
      itot = 0
      smql = 0.
C!! added balance component for upstream info.  RES 03/12/02
      smup = 0.
C Minimum simulation time
      if (up) then
        t_min = max (t(ndnew), float (nq - 1) * delt)
      else
        t_min = t(ndnew)
      end if

C**      fagl_flag = .FALSE.

C*TEMP SNOW
       IF( skip ) GO TO 100
C*TEMP SNOW

C                                 start fixed time step loop (i, delt, step)...
C------------------------------------------------------------------------------

      do i = 2, limit

        itm = i
        qub2 = 0.

        dvf_delt = 0.
        rf_delt = 0.
C                                                                        inflow
        if (up .and. i .le. nq) call get (upq, i, qub2)

        dqub = qub2 - qub1
        qup(1,1) = qub2
        quba = 0.5*(qub2 + qub1)
        step = float (i - 1) * delt
        dtp = delt

        elapsed_time = step

C                                start variable time step loop (k, dt, time)...
C------------------------------------------------------------------------------
C!!  change limiting error band: (could be even larger)

20      if (time .eq. t(k)) then
C                                                   last time = last breakpoint
C                                              next time = next fixed time step
          rf = r(k)
          dtp = step - time
          k = k + 1

        end if
C                                            time elapsed since previous "step"
        dstep = time - (step - delt)

        if (time .lt. t(k) .and. time + dtp .gt. t(k)) then
C                                                                   next time =
C                                                               next breakpoint
          dtp = t(k) - time
          time = t(k)

        else

          time = step

        end if
C                                         skip this time step if it's too small
        if (dtp .lt. DTMIN) then
            if (time .ne. step) then
              go to 20
          else
            cycle
            end if
        end if

        if (cour) then
C                                                     compute next dt using max
C                                                         celerity from last dt
          nc = nint (vcm * dtp / dx)
          if (nc .eq. 0) nc = 1
          dt = dtp / float (nc)

          if (dt .lt. DTMIN) then
            nc = int (dtp / DTMIN)
            if (nc .eq. 0) nc = 1
            dt = dtp / float (nc)
          end if

        else

          dt = dtp

        end if

        dte = 0.
        sumrf = sumrf  + rf*dt
C
        do j = 2, nk
          jm = j-1
          afac(j) = 1.
          qui = q1(jm)
          if(j .eq. 2) qui = 0.5*(qub1+qub2)/width
          h1avg(j) = 0.5*(h1(j) + h1(jm))
          if(h1avg(j) .gt. 0.) then
            if( amr .gt. 0.) then
C!!new formulation                     reduce infiltration during recession
              afac(j) = 0.5*(topw(j) + topw(jm))
            end if
          end if
C               temporarily set ql equal to inflow rate only for passing to INFILT
          ql(j,1) = qui/dx
        end do
C                                                 get average infiltration rate
C!!                               afac returns to represent coeff. for qlat
C!!    replace dt with dtp/ 2/18/03
        call infilt(afac, rf, h1avg, ql, dtp, fav, 1)

        do j = 2, nk
C                                               compute net lateral inflow rate
            diff = rf - fav(j)
            reldif = abs(diff)
            if(reldif .lt. 1.e-10) then      ! redefine ql() array
              ql(j,1) = 0.
            else
              ql(j,1) = diff
            end if
C
        end do

        if (amr .gt. 0.) then
C                                          reduce infiltration during recession 
          do j = 2, nk
              if(ql(j,1) .lt. 0.) then
                ql(j,1) = ql(j,1) * afac(j)
              end if
          end do
        end if
C                                                 start Courant loop (dt, m)...
C------------------------------------------------------------------------------

50      do m = 1, nc

          vcm = 0.
C                                           time elapsed since beginning of dtp
          dte = dte + dt
C                                                time elapsed since last "step"

          deltr = dstep + dte
          qub = qub1 + dqub * deltr / delt
C                                                            interpolate inflow
          if (qub .gt. 0.) then

            q2(1) = qub / width        ! units L^2/T
            unif = .false.
C                                               compute upstream boundary depth
C------------------------------------------------------------------------------
            if (amr .lt. 1.E-8) then
C                                                            no microtopography
              h2(1) = (q2(1) / alpha)**(1. / p)
C                                                                      celerity
              vc = p * alpha * h2(1)**(p - 1.)

            else
C                                                               microtopography
              qb = qub * fnw

              if (qb .ge. uqc) then

                a2u = (qb/ alpha * wpc**(p - 1.))**(1. / p)
                vc  = p * alpha * (a2u /wpc)**(p - 1.)

              else

                if(qub .gt. 5. * qub1)
     &            a2u = ((qb / alpha) * (2. * c2 * c2 / c1)**
     &                  (.5 * (p - 1.)))**(1. / (.5 * (p + 1.)))

                source = 'PLanUB'

                call iter (erpaub, a2u, 0., h_lim, ierr, source)

                if (ierr .gt. 0) call errxit
     &                           (msg(ms:14), 'no convergence (erpaub)')

                vc = dqdaf (alpha, a2u, bw, p, c1, c2)

              end if

              h2(1) = a2u / ams

            end if

          else

            q2(1) = 0.
            h2(1) = 0.
            vc = 0.

          end if

          if (vc .gt. vcm) vcm = vc
C                                                               zone A solution
C------------------------------------------------------------------------------
          if (unif) then

            vc1 = 0.

            do j = 2, nk

              h2(j) = h1(j) + ql(j,1) * dt

              if (h2(j) .le. 0.) then

                h2(j) = 0.
                q2(j) = 0.
                vc2   = 0.

              else if (h2(j) .ge. h2c) then

                a2j   = h2(j) * ams
                q2(j) = alpha * a2j**p / wpc**(p - 1.) / ams
                vc2   = p * alpha * (a2j /wpc)**(p - 1.)

              else

                a2j   = h2(j) * ams
                q2(j) = qfunc (alpha, a2j, bw, p, c1, c2) / ams
                vc2   = dqdaf (alpha, a2j, bw, p, c1, c2)

              end if

              vc    = (vc1 + vc2)/ 2.
              xc(j) = xc(j) + vc * dt
              vc1   = vc2

              if (vc .gt. vcm) vcm = vc

              if (xc(j) > xc_max) xc_max = xc(j)

            end do

            if (xc_max .lt. dx) then
C                                        zone A - bypass iterative calculations
              go to 70

            else
C                                                     solution is out of zone A
              unif = .false.

              go to 55

            end if

          end if

C                                                         start spatial loop...
C------------------------------------------------------------------------------

55        do j = 1, nk - 1

            jp1 = j + 1
            htest = h2(j) + h1(jp1) + h1(j)
            qlav = ql(jp1,1) !0.5 * (ql(j) + ql(jp1))
C                                                               test for runoff
            if (htest .le. 1.e-8 .and. qlav .le. 1.e-11) then

              h2(jp1) = 0.
              q2(jp1) = 0.
C!!
              ql(j,1) = 0.
              vc = 0.

            else
C                                                iterative solution for h2(j+1)
              h2jp1 = (h2(j) + h1(jp1)) / 2.

              write(source,"('Plan',i2)") jp1

              call iter (errh2, h2jp1, 0., h_lim, ierr, source)

              if (ierr .gt. 0) then
                if(diag) write(99,*) jp1, hmin,h2jp1,ql(j,1),ql(jp1,1)
                call errxit(msg(ms:14), 'no convergence (errh2)')
              else if (j .eq. 2) then
                itot = itot - ierr
              end if

C              if(j .eq. 1) write(99,'(i3," itrs")') ierr

              if (h2jp1 .lt. 1.E-7) h2jp1 = 0.

              h2(jp1) = h2jp1
              a2jp1   = h2jp1 * ams
C                                                              compute celerity
              if (h2jp1 .ge. h2c) then
                vc = p*alpha*(a2jp1 /wpc)**(p-1.)
              else
                vc = dqdaf(alpha, a2jp1, bw, p, c1, c2)
              end if

            end if
C                                           vcm is max. celerity for current dt
            if (vc .gt. vcm) vcm = vc

          end do
C                                                       ... end of spatial loop
C------------------------------------------------------------------------------

C                                                       check Courant condition
C------------------------------------------------------------------------------

          if (vcm * dt / dx .gt. 1. .and. dt .gt. DTMIN) then
C                                                                     reduce dt
            dtprev = dt
            dtr = dtp - dte + dt
            nc = int (2. * vcm * dtr / dx)
            if (nc .eq. 0) nc = 1
            dt = dtr / float (nc)

            if (dt .lt. DTMIN) then
              nc = int (dtr / DTMIN)
              if (nc .eq. 0) nc = 1
              dt = dtr / float (nc)
            end if

            if (cour) then
C                                                recompute with new subinterval
              dte = dte - dtprev
              go to 50

            else
C                                             find minimum dt over each quarter
              if (i .le. p1) then
                if (dt .lt. dtm(1)) dtm(1) = dt
              else if (i .le. p2) then
                if (dt .lt. dtm(2)) dtm(2) = dt
              else if (i .le. p3) then
                if (dt .lt. dtm(3)) dtm(3) = dt
              else
                if (dt .lt. dtm(4)) dtm(4) = dt
              end if
C                                                               don't recompute
              nc = 1
              dt = dtp

            end if
          end if

70        continue
C
C!! added for balance component RES  03/12/02:
          smup = smup + 0.5*(q1(1)+q2(1))*dt
C
          hsum = 0.5*(h2(1) + h2(nk))
C
          do j = 1, nk
C!!           moved up from sed option below              calculate top widths
            if (h2(j) .lt. h2c) then
              a2u = h2(j) * ams
              z2 = zfu (bw, a2u, c1)
              topw(j) = z2 / ams
            else
              topw(j) = 1.
            end if

            if(j .ge. 2) then
C                                                            infiltrated volume
              if(j .lt. nk) hsum = hsum + h2(j)

              if (h2(j) .eq. 0.) then
                dvf = 0.5 * ((h2(j-1) - h1(j) - h1(j-1)) / dt -
     &                       (q2(j-1) + q1(j-1) - q1(j)) / dx) + rf
              else
                dvf = rf -  ql(j,1)    !)0.5 * (ql(j-1) +
              end if

              qbal(6) = qbal(6) + dt * dvf
              smql = smql +  ql(j,1) * dt * dx

              dvf_delt = dvf_delt + dt * dvf
              rf_delt = rf_delt + dt * rf
              
            end if

            h1(j) = h2(j)
            q1(j) = q2(j)
            
          end do
C
          storh = hsum*dx
          smin = qbal(6)*dx
          qbal(10) = qbal(10) + dt * 0.5 * (q1(nk) + q2(nk))  !outflow summation
C!! added smup to balance calcs  RES  03/12/02
          sterr =  smql + smup - qbal(10) - storh   !sumrf*xleng - smin

C                                                     ending time of current dt
          timec = time - dtp + dte

          if (q2(nk) .gt. qp) then
C                                                                peak discharge
            qp = q2(nk)
            tp = timec

          end if

          qpw = qp * width
C          
C                                                                      sediment
C------------------------------------------------------------------------------
          if (sed) then

C                                                                        kinsed
            call kinsed ( i, dt, deltr, ql, h2, q2,
     &                   topw, rf, qup, qt, timec)

          end if
C                                                                outflow volume
        end do
C                                                        ...end of Courant loop
C------------------------------------------------------------------------------

        if (time .ne. step) go to 20
C                                             ...end of variable time step loop
C------------------------------------------------------------------------------

C                                                               store discharge
        qout2 = width * q2(nk)

        call store (outq(1), i, qout2)
C!!  revised RES 3/12/02                                                                 inflow volume
        if(i .lt. limit) vin = vin + qub2

        qub1 = qub2

        dvf_delt = dvf_delt*1000.*3600./delt/FLOAT(nk-1)
        IF( dvf_delt .LT. 0. ) dvf_delt = 0.

C After rainfall and upstream inflow cease, stop simulating when no water
C remains on the surface

        if(time .gt. t_min .and. hsum .eq. 0.) go to 100

      end do
C                                               ...end of fixed time step loop
C------------------------------------------------------------------------------

100   continue

C                                                  finish sediment computations
      if (sed) call sedfin
C                                                                 inflow volume
      vin = (delt / 2.) * (2. * vin + qub2)
C                                                                outflow volume
      qbal(10) = qbal(10) * width

      do j = 2, nk - 1
C                                                                storage volume
        qbal(9) = qbal(9) + h2(j)

      end do

      qbal(9) = (dx / 2.) * (h2(1) + 2. * qbal(9) + h2(nk)) * width

C                                                            infiltrated volume
      qbal(6) = max (qbal(6) * width * dx, 0.)
C                                                         pass values to writer

C      call qwrt (id, 1, msg, qpw, tp, vin, coarea, outq, 0, outr, qrmax)

      return

      end

C------------------------------------------------------------------------------

      subroutine errh2 (h2jp1, ferr, dferr)

C   Computes the residual ferr and derivative dferr at h2jp1 of the error
C   function form of the finite difference approximation of the kinematic
C   equation for flow on a sloping plane.  Includes microtopog. geometry

      common /plane1/ dt, dx, j, jp1, qlav, alpha, p, h1(20), h2(20), 
     &                q1(20), q2(20), h2c, ams, wpc, c1, c2, bw, qb,
     &                omega

      dhdt = .5 * (h2jp1 - h1(jp1) + h2(j) - h1(j)) / dt

      a2jp1 = max( 0., h2jp1 * ams )

      if (h2jp1 .ge. h2c) then

        q2(jp1) = alpha * a2jp1**p/ wpc**(p - 1.) / ams
        dferr = .5 /dt + omega * alpha * p * (a2jp1/ wpc)**(p - 1.) /dx
      else

        q2(jp1) = qfunc (alpha, a2jp1, bw, p, c1, c2) / ams
        dferr   = .5 /dt + omega * 
     &            dqdaf(alpha, a2jp1, bw, p, c1, c2) /dx / ams

      end if
      
      dqdx = (omega*(q2(jp1)-q2(j)) + (1.- omega)*(q1(jp1)-q1(j))) / dx 

      ferr = dhdt + dqdx - qlav

      return

      end

C------------------------------------------------------------------------------

      subroutine erpaub (aub, ferr, dferr)

C   Computes the residual and derivative of the normal discharge function
C   to obtain the upstream area for a given upstream flow for trapezoidal
C   microtopography

      common /plane1/ dt, dx, j, jp1, qlav, alpha, p, h1(20), h2(20), 
     &                q1(20), q2(20), h2c, ams, wpc, c1, c2, bw, qb,
     &                omega

      ferr  = qfunc(alpha, aub, bw, p, c1, c2) - qb
      dferr = dqdaf (alpha, aub, bw, p, c1, c2)

      return

      end

C------------------------------------------------------------------------------

      real function qfunc (al, a, w, p, c1, c2)

C   Computes discharge as a function of area 'a' and bottom width 'w'

      if(a .gt. 0.) then
        qfunc = al * a * (a / wprf (w, a, c1, c2))**(p - 1.)
      else
        qfunc = 0.
      end if

      return

      end

C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

      real function dqdaf (al, a, w, p, c1, c2)

C   Computes derivative of 'q' with resp to 'a'

      if(a .le. 0.) then
        dqdaf = 0.
      else
        wp    = wprf (w, a, c1, c2)
        d     = zfu (w, a, c1)
        dqdaf = al*(a /wp)**(p-1.)*(p -((p-1.)*a*c2/ (d*wp)))
      end if

      return

      end

C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

      real function wprf (w, a, c1, c2)

C   Computes the wetted perimeter of a trazpezoid of bottom width 'w'

      wprf = w + c2 * (zfu(w, a, c1) - w) / c1

      return

      end

C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

      real function zfu (w, a, c1)

C   Computes geometric shape parameter 'z'

      if( a .ge. 0.) then
        zfu = sqrt (w * w + a * 2. * c1)
      else
        zfu = w
      end if

      return

      end
