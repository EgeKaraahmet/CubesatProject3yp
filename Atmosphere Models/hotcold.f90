!+
! PROGRAM TestNonStandard
! ------------------------------------------------------------------------------
! PURPOSE - Test the procedures in modules MILSTD210A and ModifiedStandard by
!  building tables of the properties of these atmospheres at selected altitudes.
! AUTHOR - Ralph L. Carmichael, Public Domain Aeronautical Software

! NOTE Needs the Modified standard for comparison
!
!     REVISION HISTORY
!   DATE  VERS PERSON  STATEMENT OF CHANGES
! 03Dec10 0.1    RLC   Original coding
! 12Dec10 0.2    RLC   Renamed files; MIL atmos keep constant temp. above 100kft
!-------------------------------------------------------------------------------


!+
MODULE ModifiedStandard
! ------------------------------------------------------------------------------
! PURPOSE - Hold the standard atmosphere procedures
IMPLICIT NONE
!-------------------------------------------------------------------------------

CONTAINS

!+
SUBROUTINE ModifiedStandardAtmosphere(alt, sigma, delta, theta, deltaT)
!   -------------------------------------------------------------------------
! PURPOSE - Compute the properties of the 1976 standard atmosphere to 86 km.
!  allowing for a user-specified increment to the temperature to account for
!  weather variations (hot day, cold day, etc.)

  IMPLICIT NONE
!============================================================================
!     A R G U M E N T S                                                     |
!============================================================================
  REAL,INTENT(IN)::  alt    ! geometric altitude, km.                        
  REAL,INTENT(OUT):: sigma  ! density/sea-level standard density              
  REAL,INTENT(OUT):: delta  ! pressure/sea-level standard pressure           
  REAL,INTENT(OUT):: theta  ! temperature/sea-level standard temperature
  REAL,INTENT(IN),OPTIONAL:: deltaT ! temp increment, kelvins
!============================================================================
!     L O C A L   C O N S T A N T S                                         |
!============================================================================
  REAL,PARAMETER:: REARTH = 6378.0                 ! radius of the Earth (km)
  REAL,PARAMETER:: GMR = 34.163195                             ! gas constant
  INTEGER,PARAMETER:: NTAB=8       ! number of entries in the defining tables
!============================================================================
!     L O C A L   V A R I A B L E S                                         |
!============================================================================
  INTEGER:: i,j,k                                                  ! counters
  REAL:: h                                       ! geopotential altitude (km)
  REAL:: tgrad, tbase      ! temperature gradient and base temp of this layer
  REAL:: tlocal                                           ! local temperature
  REAL:: deltah                             ! height above base of this layer
!============================================================================
!     L O C A L   A R R A Y S   ( 1 9 7 6   S T D.  A T M O S P H E R E )   |
!============================================================================
  REAL,DIMENSION(NTAB),PARAMETER:: htab= &
                          (/0.0, 11.0, 20.0, 32.0, 47.0, 51.0, 71.0, 84.852/)
  REAL,DIMENSION(NTAB),PARAMETER:: ttab= &
          (/288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65, 186.946/)
  REAL,DIMENSION(NTAB),PARAMETER:: ptab= &
               (/1.0, 2.233611E-1, 5.403295E-2, 8.5666784E-3, 1.0945601E-3, &
                                     6.6063531E-4, 3.9046834E-5, 3.68501E-6/)
  REAL,DIMENSION(NTAB),PARAMETER:: gtab= &
                                (/-6.5, 0.0, 1.0, 2.8, 0.0, -2.8, -2.0, 0.0/)
!----------------------------------------------------------------------------
  h=alt*REARTH/(alt+REARTH)      ! convert geometric to geopotential altitude

  i=1 
  j=NTAB                                  ! setting up for binary search
  DO
    k=(i+j)/2                                              ! integer division
    IF (h < htab(k)) THEN
      j=k
    ELSE
      i=k
    END IF   
    IF (j <= i+1) EXIT
  END DO

  tgrad=gtab(i)                                     ! i will be in 1...NTAB-1
  tbase=ttab(i)
  deltah=h-htab(i)
  tlocal=tbase+tgrad*deltah
  theta=tlocal/ttab(1)                                    ! temperature ratio

  IF (tgrad == 0.0) THEN                                     ! pressure ratio
    delta=ptab(i)*EXP(-GMR*deltah/tbase)
  ELSE
    delta=ptab(i)*(tbase/tlocal)**(GMR/tgrad)
  END IF
  
  IF (Present(deltaT)) THEN
    theta=(tlocal+deltaT)/ttab(1)                     ! temperature increment
  END IF

  sigma=delta/theta                                           ! density ratio
  RETURN
END Subroutine ModifiedStandardAtmosphere   ! ----------------------------------

!+
FUNCTION MetricViscosity(theta) RESULT(visc)
!   ----------------------------------------------------------------------------
! PURPOSE - Compute viscosity using Sutherland's formula.
!        Returns viscosity in kg/(meter-sec)

  REAL,INTENT(IN) :: theta                ! temperature/sea-level temperature  
  REAL:: visc
  REAL:: temp                              ! temperature in deg Kelvin

  REAL,PARAMETER:: TZERO = 288.15    ! sea level standard temp, kelvins
  REAL,PARAMETER:: BETAVISC = 1.458E-6 ! viscosity term, N sec/(sq.m sqrt(deg K)
  REAL,PARAMETER:: SUTH = 110.4              ! Sutherland's constant, kelvins

!-------------------------------------------------------------------------------
  temp=TZERO*theta
  visc=BETAVISC*Sqrt(temp*temp*temp)/(temp+SUTH)
  RETURN
END Function MetricViscosity   ! -----------------------------------------------



END Module ModifiedStandard   ! ================================================

!+
MODULE MILSTD210A
! ------------------------------------------------------------------------------
! PURPOSE - Define the four atmospheres defined in MIL-STD-210A.
!  The four atmospheres: hot, cold, polar, and tropical are defined in US units.
!  The temperatures are from tables at selected geometric altitudes from sea
!  level to 100,000 ft. The pressure is assumed to be the same as standard and
!  the density follows from the perfect gas law.
USE ModifiedStandard
IMPLICIT NONE

  REAL,PARAMETER:: FT2METERS = 0.3048     ! mult. ft. to get meters (exact)
  REAL,PARAMETER:: SQFT2SQMETERS = 0.3048**2   ! convert square feet to square meters
  REAL,PARAMETER:: CUFT2CUMETERS = 0.3048**3   ! convert cubic feet to cubic meters
  REAL,PARAMETER:: SLUGS2KG = 14.5939029
  REAL,PARAMETER:: GUS = 32.174  ! gravity in ft s^-2
  REAL,PARAMETER:: LBS2NEWTONS = 4.44822165
  REAL,PARAMETER:: PASCAL2PSF = 0.02089
  REAL,PARAMETER:: PSF2PASCALS = 47.880258
  REAL,PARAMETER:: SLUGPERCUFT2KGPERCUM = 515.379
  REAL,PARAMETER:: KGPERCUM2SLUGSPERCUFT = 1.0/SLUGPERCUFT2KGPERCUM

  REAL,PARAMETER:: PI = 3.14159265
  REAL,PARAMETER:: REARTH = 6356.766               ! radius of the Earth (km)
!  REAL,PARAMETER:: GMR = 34.163195                             ! gas constant
  REAL,PARAMETER:: GZERO = 9.80665 !  accel. of gravity, m/sec^2


  REAL,PARAMETER:: KELVIN2RANKINE = 1.8           ! mult kelvins to get deg R
  REAL,PARAMETER:: PSF2NSM = 47.880258          ! mult lb/sq.ft to get N/sq.m
  REAL,PARAMETER:: SCF2KCM = 515.379         ! mult slug/cu.ft to get kg/cu.m
  REAL,PARAMETER:: TZERO = 288.15            ! temperature at sealevel, deg K
  REAL,PARAMETER:: PZERO = 101325.0            ! pressure at sealevel, N/sq.m
  REAL,PARAMETER:: RHOZERO = 1.2250            ! density at sealevel, kg/cu.m
  REAL,PARAMETER:: RSTAR = 8314.32             ! perfect gas constant
  REAL,PARAMETER:: ASOUNDZERO = 340.294   ! speed of sound at sealevel, m/sec

  INTEGER,PARAMETER:: NDIM=80
  REAL,DIMENSION(NDIM),PARAMETER:: ALTKFT = (/ &
     0.0,  1.0,  2.0,  3.0,  4.0,  5.0,  6.0,  7.0,  8.0,  9.0, &
    10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, &
    20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, &
    30.0, 31.0, 32.0, 33.0, 34.0, 35.0, 36.0, 37.0, 38.0, 39.0, &
    40.0, 41.0, 42.0, 43.0, 44.0, 45.0, 46.0, 47.0, 48.0, 49.0, &
    50.0, 51.0, 52.0, 53.0, 54.0, 55.0, 56.0, 57.0, 58.0, 60.0, &
    62.0, 64.0, 66.0, 68.0, 70.0, 72.0, 74.0, 76.0, 78.0, 80.0, &
    82.0, 84.0, 86.0, 88.0, 90.0, 92.0, 94.0, 96.0, 98.0,100.0 /)

  REAL,DIMENSION(NDIM),PARAMETER:: HOT = (/ &
    562.7, 558.9, 555.1, 551.2, 547.3, 543.4, 539.5, 535.5, 531.5, 527.5, &
    523.6, 519.9, 516.1, 512.3, 508.5, 504.6, 500.7, 496.8, 492.8, 488.9, &
    485.2, 481.5, 477.7, 474.0, 470.2, 466.4, 462.6, 458.7, 454.8, 451.0, &
    447.4, 443.8, 440.2, 436.5, 432.9, 429.6, 426.3, 423.0, 419.6, 416.2, &
    414.9, 415.4, 415.8, 416.2, 416.6, 417.1, 417.6, 418.0, 418.5, 419.0, &
    419.5, 419.8, 420.0, 420.2, 420.4, 420.6, 420.7, 420.9, 421.1, 421.5, &
    421.9, 422.3, 422.6, 423.6, 425.0, 426.4, 427.8, 429.2, 430.6, 432.0, &
    433.6, 435.1, 436.7, 438.4, 439.9, 441.4, 442.9, 444.6, 446.3, 448.1 /)

  REAL,DIMENSION(NDIM),PARAMETER:: COLD = (/ &
    399.7, 413.2, 426.7, 440.4, 444.7, 444.7, 444.7, 444.7, 444.7, 444.7, &
    444.7, 443.9, 440.6, 437.3, 434.0, 430.6, 427.3, 423.9, 420.5, 417.0, &
    413.6, 410.1, 406.5, 403.0, 399.4, 395.8, 392.2, 388.6, 384.9, 381.1, &
    377.4, 374.7, 374.7, 374.7, 374.7, 374.7, 374.7, 374.7, 374.7, 374.7, &
    374.7, 374.7, 374.7, 371.5, 366.4, 361.1, 355.8, 350.4, 345.0, 340.5, &
    336.8, 334.7, 334.7, 334.7, 334.7, 334.7, 334.7, 334.7, 334.7, 334.7, &
    337.6, 343.7, 349.3, 354.4, 359.2, 363.6, 365.4, 364.9, 364.4, 363.8, &
    363.1, 362.3, 361.5, 360.8, 360.0, 359.2, 358.4, 357.6, 356.7, 355.8 /)

  REAL,DIMENSION(NDIM),PARAMETER:: TROPICAL = (/ &
    549.5, 545.6, 541.7, 537.8, 534.0, 530.1, 526.2, 522.3, 518.4, 514.6, &
    510.7, 506.8, 502.9, 499.1, 495.2, 491.3, 487.5, 483.6, 479.7, 475.8, &
    472.0, 468.1, 464.2, 460.4, 456.5, 452.7, 448.8, 444.9, 441.7, 437.2, &
    433.4, 429.5, 425.6, 421.8, 417.9, 414.1, 410.2, 406.4, 402.6, 398.8, &
    395.1, 391.4, 387.7, 384.1, 380.5, 376.9, 373.4, 369.9, 366.5, 363.0, &
    359.6, 356.3, 352.9, 349.6, 348.6, 350.7, 352.9, 355.1, 357.2, 361.7, &
    366.1, 370.7, 375.3, 379.9, 384.2, 386.8, 389.4, 392.1, 394.7, 397.4, &
    400.1, 402.8, 405.5, 408.2, 410.9, 413.6, 416.3, 418.9, 421.6, 424.3 /)

  REAL,DIMENSION(NDIM),PARAMETER:: POLAR = (/ &
    444.0, 447.0, 450.1, 453.1, 453.5, 453.0, 452.4, 451.9, 451.3, 450.8, &
    450.0, 447.2, 444.3, 441.5, 438.7, 435.9, 433.0, 430.2, 427.4, 424.5, &
    421.7, 418.8, 416.0, 413.1, 410.3, 407.4, 404.5, 401.7, 398.8, 395.9, &
    393.0, 392.5, 392.2, 392.0, 391.7, 391.4, 391.2, 390.9, 390.7, 390.4, &
    390.1, 389.9, 389.6, 389.4, 389.1, 388.8, 388.6, 388.3, 388.1, 387.8, &
    387.5, 387.3, 387.0, 386.8, 386.5, 386.2, 386.0, 385.7, 385.5, 385.0, &
    384.4, 383.9, 383.4, 382.9, 382.4, 381.9, 381.4, 380.9, 380.3, 379.8, &
    379.3, 378.8, 378.3, 378.3, 378.3, 378.3, 378.3, 378.3, 378.3, 378.3 /)

!-------------------------------------------------------------------------------
CONTAINS

!+
PURE FUNCTION MIL210Temperature(katm,h) RESULT(tr)
! ------------------------------------------------------------------------------
! PURPOSE - Compute the temperature at a given point in one of the non-standard
!  atmospheres from MIL-STD-210A.
! NOTE - These atmospheres are only defined to 100000 feet. If a higher altitude
! is entered, the temperature at 100000 feet is returned. 
  INTEGER,INTENT(IN):: katm   ! =1 HOT; =2 COLD; =3 POLAR; =4 TROPICAL
  REAL,INTENT(IN):: h   ! geopotential altitude, kft

  REAL:: fract
  INTEGER:: i,j,k
  REAL:: tr
!-------------------------------------------------------------------------------

  IF (h < ALTKFT(1) ) THEN
    i=1
  ELSE IF (h >= ALTKFT(NDIM)) THEN
    i=NDIM-1
  ELSE
    i=1
    j=NDIM
    DO
      k=(i+j)/2   ! integer division
      IF (h < ALTKFT(k)) THEN
        j=k
      ELSE
        i=k
      END IF
      IF (j <= i+1) EXIT
    END DO
  END IF
  
  fract=(h-ALTKFT(i))/(ALTKFT(i+1)-ALTKFT(i))
  fract=MIN(fract,1.0)  ! can't be less than 0 or more than 1

  SELECT CASE(katm)
    CASE(1)
      tr=HOT(i) + fract*(HOT(i+1)-HOT(i))
    CASE(2)
      tr=COLD(i) + fract*(COLD(i+1)-COLD(i))
    CASE(3)
      tr=POLAR(i) + fract*(POLAR(i+1)-POLAR(i))
    CASE(4)
      tr=TROPICAL(i) + fract*(TROPICAL(i+1)-TROPICAL(i))
  END SELECT

  RETURN
END Function MIL210Temperature   ! ------------------------------------------------------

!+
SUBROUTINE MIL210Atmosphere(katm, z,sigma,delta,theta)
! ------------------------------------------------------------------------------
! PURPOSE - Compute density,pressure, and temperature at a given geometric
!  altitude in the hot atmosphere defined by MIL-STD-210A
  INTEGER,INTENT(IN):: katm  ! =1
  REAL,INTENT(IN):: z        ! geometric altitude, in feet
  REAL,INTENT(OUT):: sigma   ! ratio of density to standard sea level density
  REAL,INTENT(OUT):: delta   ! ratio of pressure to standard sea level pressure
  REAL,INTENT(OUT):: theta   ! ratio of temperature to standard sea level temperature

  REAL:: hkft    ! geopotential altitude, kft
  REAL:: hkm     ! geopotential altitude, km
  INTEGER:: i,j,k
  REAL:: tr,tk
  REAL:: zkm, zkft
!-------------------------------------------------------------------------------
  zkm = z*FT2METERS/1000
  hkm=zkm*REARTH/(zkm+REARTH)

  CALL ModifiedStandardAtmosphere(hkm,sigma,delta,theta)

  hkft=hkm/FT2METERS
  tr=MIL210Temperature(katm,hkft)
  tk=tr/KELVIN2RANKINE
  theta=tk/TZERO

  sigma=delta/theta

  RETURN
END Subroutine MIL210Atmosphere   ! --------------------------------------------

END Module MILSTD210A   ! ======================================================


!+
MODULE TestNonStandardProcedures
! ------------------------------------------------------------------------------
! PURPOSE -
USE MILSTD210A
USE ModifiedStandard
IMPLICIT NONE

  CHARACTER(LEN=*),PARAMETER:: HEAD1 = &
    "     H     T     wt dens   press    mass dens     eta         c  "
  CHARACTER(LEN=*),PARAMETER:: HEAD2 = &
    "    kft    R     lb/cuft    psf     slug/cuft    sqft/s      ft/s"

  CHARACTER(LEN=*),PRIVATE,PARAMETER:: FMT1 = &
    '(F7.1,F7.1,F10.6,F9.1,F12.7,F12.7,F9.1,F9.2)'

  INTEGER,PRIVATE:: k
  REAL,PRIVATE,PARAMETER,DIMENSION(101):: HTAB = (/ (REAL(k),k=0,100) /)
  REAL,PRIVATE,PARAMETER,DIMENSION(101):: ZTAB = HTAB*REARTH/(REARTH-HTAB)
  REAL,PRIVATE:: sigma,delta,theta
!-------------------------------------------------------------------------------

CONTAINS
 
!+
SUBROUTINE BuildMIL210Atmosphere(efu,katm)
! ------------------------------------------------------------------------------
  INTEGER,INTENT(IN):: efu
  INTEGER,INTENT(IN):: katm   ! =1 HOT; =2 COLD; =3 POLAR; =4 TROPICAL

  CHARACTER(LEN=8),DIMENSION(4):: ATM_TITLE = &
     (/ "HOT     ", "COLD    ", "POLAR   ", "TROPICAL" /)
  REAL:: c
  REAL:: eta   ! kinematic viscosity  sq.ft/s
  REAL:: hkm
  INTEGER:: k
  REAL:: omega   ! weight density, lbs per cu. ft.
  REAL:: p
  REAL:: rho     ! mass density, kg per cu.m
  REAL:: rhous   ! mass density, slugs per cu. ft.
  REAL:: tr
  REAL:: visc  ! dynamic viscosity, US units
  REAL:: zkm,zft
!-------------------------------------------------------------------------------
  WRITE(efu,'(///1X,A)') TRIM(ATM_TITLE(katm))//" ATMOSPHERE"
  WRITE(efu,'(A)') HEAD1
  WRITE(efu,'(A)') HEAD2
  
  DO k=1,SIZE(HTAB)
    hkm=FT2METERS*HTAB(k)
    zkm=hkm*REARTH/(REARTH-hkm)
    zft=1000*zkm/FT2METERS
    CALL MIL210Atmosphere(katm, zft, sigma,delta,theta)
    write(3,*) 'htab,hkm,zkm,zft', HTAB(k),hkm,zkm,zft
    write(3,*) 'sigma,delta,theta', sigma,delta,theta
    tr=theta*TZERO*KELVIN2RANKINE
!    rho = sigma*RHOZERO               ! SI units
    rhous=sigma*RHOZERO*KGPERCUM2SLUGSPERCUFT   ! US units
    omega=GUS*rhous
    p=delta*PZERO*PASCAL2PSF
    write(3,*) 'omega,p', omega,p
    visc=(1.0/PSF2NSM)*MetricViscosity(theta)  ! US units
    eta=visc/rhous   
    write(3,*) 'visc,eta', visc,eta
    c=(ASOUNDZERO/FT2METERS)*SQRT(theta)
    WRITE(efu,FMT1) HTAB(k), tr, omega, p, rhous, eta, c, zft/1000
  END DO

  RETURN
END Subroutine BuildMIL210Atmosphere   ! ---------------------------------------

END Module TestNonStandardProcedures   ! =======================================

!+
PROGRAM TestNonStandard
! ------------------------------------------------------------------------------
USE TestNonStandardProcedures
IMPLICIT NONE

  INTEGER,PARAMETER:: OUT=1, FIG=2, DBG=3, GNU=4, PS=7, SVG=8, HTML=9
  INTEGER:: k

!-------------------------------------------------------------------------------
  OPEN(UNIT=OUT,FILE='hotcold.out',STATUS='REPLACE',ACTION='WRITE')
!  OPEN(UNIT=FIG,FILE='ns.fig',STATUS='REPLACE',ACTION='WRITE')
  OPEN(UNIT=DBG,FILE='hotcold.dbg',STATUS='REPLACE',ACTION='WRITE')
!  OPEN(UNIT=GNU,FILE='ns.gnu',STATUS='REPLACE',ACTION='WRITE')
!  OPEN(UNIT=HTML,FILE='ex210.html',STATUS='REPLACE',ACTION='WRITE')
  WRITE(DBG,*) 'Checking Non-standard atmospheres'

  DO k=1,4
    CALL BuildMIL210Atmosphere(OUT,k)
  END DO  
  WRITE(*,*) 'The files ns.dbg and ns.fig may be erased if the graphics are OK.'
  STOP
END Program TestNonStandard   ! ================================================
