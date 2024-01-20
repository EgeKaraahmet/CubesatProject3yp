!+
!PROGRAM BigTables                                      ! /atmos/bigtables.f90
! ------------------------------------------------------------------------------
! PURPOSE - Write a file in HTML 5 format that displays tables of atmospheric 
!   properties from 0 to 1000 km geometric altitude.
!   There are six different tables...
!   The tables are based on those from the reference document, Part 4.

!   1. Geopotential altitude, temperature, temperature ratio, 
!      pressure, pressure ratio, density, density ratio, 
!        sound speed, gravity for geometric altitudes in kilometers.
!
!   2. Geopotential altitude, dynamic viscosity, kinematic viscosity,  
!       pressure scale height, number density, mean particle speed,  
!       mean collision frequency, mean free path,  
!       thermal conductivity coefficient and molecular weight for  
!       geometric altitudes in kilometers.
!  
!   3. Non-dimensional ratios of dynamic viscosity,  
!       kinematic viscosity, number density, particle speed,  
!       mean collision frequency, mean free path, thermal conductivity,  
!       sound speed and molecular weight for altitudes in kilometers.  
!       All quantities are referenced to sea-level values.
  
!   4. Geopotential altitude, temperature, temperature ratio,  
!       pressure, pressure ratio, mass density, density ratio,  
!       weight density, sound speed and gravity for geometric altitudes  
!       in thousand feet. Table entries in US customary units.
  
!   5. Geopotential altitude, dynamic viscosity, kinematic viscosity,  
!       pressure scale height, number density, mean particle speed,  
!       mean collision frequency, mean free path,  
!       thermal conductivity coefficient and molecular weight for  
!       geometric altitudes in thousand feet.  
!       Table entries in US customary units.

!   6. Non-dimensional ratios of dynamic viscosity,  
!       kinematic viscosity, number density, particle speed,  
!       mean collision frequency, mean free path, thermal conductivity,  
!       sound speed and molecular weight for altitudes in thousand feet.  
!       All quantities are referenced to sea-level values.

!
! REFERENCE - "The U.S. Standard Atmosphere, 1976" published by the 
!    Government Printing Office. Also filed as NASA TM X-74335.

! AUTHOR - Ralph L. Carmichael, Public Domain Aeronautical Software


!     REVISION HISTORY
!   DATE  VERS PERSON  STATEMENT OF CHANGES
! 28Feb95  1.0   RLC   Assembled several old codes
!  7Jul95  1.1   RLC   Renamed routines, coded viscosity in US units
! 26Jul95  1.2   RLC   Converted f77 program to f90 style.
!                      Made Module PhysicalConstants
!  6Aug95  1.3   RLC   Replaced 1962 tables with 1976 tables
! 25Aug95  1.4   RLC   Added MODIFIER
! 19Sep95  1.5   RLC   Added a little precision to pressure tables
!  3Oct95  1.6   RLC   Numerous style changes to run under Elf90
! 29Nov96  1.7   RLC   Error message if unable to open output file
!  1Jan01  1.8   RLC   Converted to BigTable, going to 1000 km.
! 06Jul18  1.81  RLC   Revised contents of tables
! 31Dec10  1.9   RLC   Big revision of module atmosphere
! 06Jan11  1.95  RLC   Output in HTML. Created MakeTable1.
! 08Jan11  1.96  RLC   Added MakeTable3
! 11Jul18  1.97  RLC   Redesigned to put all 6 tables in one file
! 24Jul18  1.98  RLC   Added molecular weight calculations 
! 18Mar22  2.0   RLC   Made sure REARTH = 6356.766km

! NOTES - It is best to compile with a switch that makes the default real 
!  precision 64 bits. With gfortran, this is done with
!     -fdefault-real-8 

!-------------------------------------------------------------------------------


!+
MODULE Atmosphere1976
! ---------------------------------------------------------------------------
! PURPOSE - Compute properties of the U.S. Standard Atmosphere 1976
! AUTHORS - Steven S. Pietrobon. Small World 
!           Ralph L. Carmichael, Public Domain Aeronautical Software
! REFERENCE - "The U.S. Standard Atmosphere, 1976" published by the 
!    Government Printing Office. Also filed as NASA TM X-74335.
!    NASA SP-398 is another useful reference.
!
!     REVISION HISTORY
!   DATE  VERS PERSON  STATEMENT OF CHANGES
! 28Feb95  0.1   RLC   Assembled several old codes
!  1Aug00  0.2   RLC   Copied from old Tables76
! 23Aug00  0.3   RLC   Added NitrogenNumber using QUANC8 (later removed)
! 24Aug00  0.4   RLC   Added KineticTemperature
! 30Aug00  0.5   RLC   Corrected numerous errors
! 30Dec00  0.6   RLC   Adapted UpperAtmosphere from Pietrobon's Unofficial
!                        Australian Standard Atmosphere
! 30Jun18  0.7   RLC   Use tables of sigma and delta rather than p and rho 
! 19Jul18  0.8   RLC   Added subroutine TransportRatios
! 24Jul18  0.9   RLC   Added function MolecularWeight 
!-------------------------------------------------------------------------------
IMPLICIT NONE
  CHARACTER(LEN=*),PUBLIC,PARAMETER:: ATM76_VERSION = "0.9 (24 July 2018)"
  REAL,PRIVATE,PARAMETER:: PI = 3.14159265358979323846264338
!
! various choices for the radius of the ideal spherical earth...
! The value of POLAR_RADIUS is used, to agree with the reference document.  
  REAL,PARAMETER,PRIVATE:: POLAR_RADIUS = 6356.766
  REAL,PARAMETER,PRIVATE:: EQUATORIAL_RADIUS = 6378.1370
  REAL,PARAMETER,PRIVATE:: RADIUS_AT_45DEG_LAT = 6367.4895
  REAL,PARAMETER,PRIVATE:: AUTHALIC_RADIUS = 6371.0012               ! same area 
  REAL,PARAMETER,PRIVATE:: VOLUMETRIC_RADIUS = 6371.0008           ! same volume
  REAL,PARAMETER:: REARTH = POLAR_RADIUS              ! radius of the Earth (km)

! Some physical constants and properties of air
  REAL,PARAMETER:: GZERO = 9.80665         !  sea level accel. of gravity, m/s^2
  REAL,PARAMETER:: MOLWT_ZERO = 28.9644    ! molecular weight of air at sea level
  REAL,PARAMETER:: RSTAR = 8314.32          ! perfect gas constant, N-m/(kmol-K)
  REAL,PARAMETER:: GMR=1000*GZERO*MOLWT_ZERO/RSTAR ! hydrostatic constant, kelvins/km
! NOTE - Observe the factor of 1000 in computing GMR. Without it, GMR would have 
!   units of kelvins per meter. But the temperature gradient in LowerAtmosphere
!   is given in kelvins per kilometer.   GMR=34.163195 kelvins/km

  REAL,PARAMETER:: BETAVISC = 1.458E-6       ! viscosity term, N s/(sq.m sqrt(K)
  REAL,PARAMETER:: SUTH = 110.4                       ! Sutherland's constant, K
  REAL,PARAMETER:: AVOGADRO =  6.022169E26           ! 1/kmol, Avogadro constant
  REAL,PARAMETER:: BOLTZMANN = 1.380622E-23           ! Nm/K, Boltzmann constant
  
! Some sea level values of various quantities  
  REAL,PARAMETER:: TZERO = 288.15                   ! temperature at sea level, K
  REAL,PARAMETER:: PZERO = 101325.0               ! pressure at sea level, N/sq.m
  REAL,PARAMETER:: RHOZERO = 1.2250               ! density at sea level, kg/cu.m
  REAL,PARAMETER:: ASOUNDZERO = 340.294      ! speed of sound at sea level, m/sec
  REAL,PARAMETER:: MUZERO = 1.7894E-5   ! viscosity at sea-level, 
  REAL,PARAMETER:: ETAZERO = 1.4607E-5 ! kinetic viscosity at sea-level, 
  REAL,PARAMETER:: KAPPA_ZERO = 0.025326            ! thermal coeff. at sea-level
                             ! units of KAPPA_ZERO are watts per meter per kelvin
  REAL,PARAMETER:: CROSS = 3.65E-10 ! collision cross-section, diatomic air m^2
  REAL,PARAMETER:: NUM_DENSITY_ZERO = RHOZERO*AVOGADRO/MOLWT_ZERO !  1/m^3
  REAL,PARAMETER:: PART_SPEED_ZERO_SQ = (8.0/PI)*RSTAR*TZERO/MOLWT_ZERO! m^2/s^2
  REAL,PARAMETER:: PART_SPEED_ZERO = SQRT((8.0/PI)*RSTAR*TZERO/MOLWT_ZERO) ! m/s

  REAL,PARAMETER:: FREE_PATH_ZERO = &
                 RSTAR*TZERO/(SQRT(2.0)*PI*AVOGADRO*CROSS*CROSS*PZERO)       ! m
  REAL,PARAMETER:: PRESSURE_SCALE_HEIGHT_ZERO = RSTAR*TZERO/(MOLWT_ZERO*GZERO)! m

! Some conversion constants for other than SI values
!!! ??  REAL,PARAMETER:: KILOWATTHOURS2KCAL = 860.0            ! exact. the definition
  REAL,PARAMETER:: BTU2KCAL = 0.2520                      ! mult BTU to get kcal
  REAL,PARAMETER:: FT2METERS = 0.3048          ! mult. ft. to get meters (exact)
  REAL,PARAMETER:: KELVIN2RANKINE = 1.8              ! mult kelvins to get deg R
  REAL,PARAMETER:: PSF2NSM = 47.880258             ! mult lb/sq.ft to get N/sq.m
  REAL,PARAMETER:: SCF2KCM = 515.379            ! mult slug/cu.ft to get kg/cu.m
  REAL,PARAMETER:: INHG2NSM = 3376.85                  ! mult inHg to get N/sq.m
  REAL,PARAMETER:: SLUG2KG = 14.594
  REAL,PARAMETER:: BTU2JOULE = 1055.0
  REAL,PARAMETER:: KCAL2JOULE = 4184.0
  REAL,PARAMETER:: INWATER2PASCAL = 248.84
  REAL,PARAMETER:: MMHG2PASCAL = 133.3224
  REAL,PARAMETER:: BARREL2LITER = 158.9873
  REAL,PARAMETER:: CALORIE2JOULE = 4.1868
  REAL,PARAMETER:: TORR2PASCAL = 133.3224   
  REAL,PARAMETER:: BAR2PASCAL = 1.0E5      

  PUBLIC:: Atmosphere
  PRIVATE:: EvaluateCubic
  PRIVATE:: KineticTemperature
  PRIVATE:: LowerAtmosphere
  PRIVATE:: UpperAtmosphere
  PUBLIC:: SimpleAtmosphere
  PUBLIC:: ThermalConductivity
  PUBLIC:: Viscosity


CONTAINS

!+
PURE FUNCTION EvaluateCubic(a,fa,fpa, b,fb,fpb, u) RESULT(fu)
! ---------------------------------------------------------------------------
! PURPOSE - Evaluate a cubic polynomial defined by the function and the
!   1st derivative at two points.
!       Usually  a < u < b, but not necessary.  
  REAL,INTENT(IN):: u   ! point where function is to be evaluated
  REAL,INTENT(IN):: a,fa,fpa   ! a, f(a), f'(a)  at first point
  REAL,INTENT(IN):: b,fb,fpb   ! b, f(b), f'(b)  at second point


  REAL:: fu                    ! computed value of f(u)

  REAL:: d,t,p
!----------------------------------------------------------------------------
  d=(fb-fa)/(b-a)   ! b==a will be a fatal error 
  t=(u-a)/(b-a)
  p=1.0-t

  fu = p*fa + t*fb - p*t*(b-a)*(p*(d-fpa)-t*(d-fpb))
  RETURN
END Function EvaluateCubic   ! -------------------------------------------------

!+
PURE FUNCTION KineticTemperature(z) RESULT(t)
!   ----------------------------------------------------------------------------
! PURPOSE - Compute kinetic temperature above 86 km.

  REAL,INTENT(IN)::  z     ! geometric altitude, km.                        
  REAL:: t     ! kinetic temperature, K

! TABLE 5 - DEFINITION OF KINETIC TEMPERATURE FROM 86 km to 1000 km
  REAL,PARAMETER:: C1 = -76.3232  ! uppercase A in document
  REAL,PARAMETER:: C2 = 19.9429   ! lowercase a in document
  REAL,PARAMETER:: C3 = 12.0
  REAL,PARAMETER:: C4 = 0.01875   ! lambda in document
  REAL,PARAMETER:: TC = 263.1905
  REAL,PARAMETER:: Z7 =  86.0   ! not used  
  REAL,PARAMETER:: T7=186.8673
  REAL,PARAMETER:: z8 =  91.0
  REAL,PARAMETER:: Z9 = 110.0,  T9=240.0
  REAL,PARAMETER:: Z10= 120.0, T10=360.0
!  REAL,PARAMETER:: Z11= 500.0, T11=999.2356   ! not used
  REAL,PARAMETER:: T12=1000.0
  
! NOTE - REARTH is a module variable  
   ! same as T sub infinity

  REAL:: xx,yy
!-------------------------------------------------------------------------------
  IF (z <= Z8) THEN
    t=T7                    ! Eq (25), p.11
  ELSE IF (z < Z9) THEN  
    xx=(z-Z8)/C2                        ! from Appendix B, p.223
    yy=SQRT(1.0-xx*xx)
    t=TC+C1*yy            ! Eq (27), p.11
  ELSE IF (z <= Z10) THEN
    t=T9+C3*(z-Z9)      ! Eq (29)
  ELSE
    xx=(REARTH+Z10)/(REARTH+z)
    yy=(T12-T10)*EXP(-C4*(z-Z10)*xx)
    t=T12-yy       ! Eq. (31)
  END IF

  RETURN
END Function KineticTemperature   ! --------------------------------------------

!+  
PURE FUNCTION ThermalConductivity(t) RESULT(k)   
!   ----------------------------------------------------------------------------
! PURPOSE - Compute the coefficient of thermal conductivity at temperature t.
! This is a coding of section 1.3.13 on p. 19 of the reference document.

  REAL,INTENT(IN) :: t                ! temperature, K

  REAL:: k  ! coefficient of thermal conductivity, watts per meter per kelvin
    ! some refer to the units as joules per second per meter per kelvin
  REAL,PARAMETER:: C1=2.64638E-3, C2=245.4
! NOTE - This empirical equation is given in the reference document.
!-------------------------------------------------------------------------------
  k=C1*SQRT(t*t*t)/(t+C2*10.0**(-12.0/t))
  RETURN
END Function ThermalConductivity   ! -------------------------------------------

!+
SUBROUTINE UpperAtmosphere(alt, sigma, delta, theta)
!   ----------------------------------------------------------------------------
! PURPOSE - Compute the properties of the 1976 standard atmosphere from
!   86 km. to 1000 km.

  REAL,INTENT(IN)::  alt    ! geometric altitude, km.                        
  REAL,INTENT(OUT):: sigma  ! density/sea-level standard density              
  REAL,INTENT(OUT):: delta  ! pressure/sea-level standard pressure           
  REAL,INTENT(OUT):: theta  ! temperature/sea-level standard temperature

  INTEGER,PARAMETER:: NTABLE = 25

! altitude table (km)
  REAL,PARAMETER,DIMENSION(NTABLE):: Z = (/      &
     86.0,  93.0, 100.0, 107.0, 114.0, &
    121.0, 128.0, 135.0, 142.0, 150.0, &
    160.0, 170.0, 180.0, 190.0, 200.0, &
    220.0, 260.0, 300.0, 400.0, 500.0, &
    600.0, 700.0, 800.0, 900.0,1000.0 /)

! pressure table  (DELTA = P/P0)
  REAL,PARAMETER,DIMENSION(SIZE(Z)):: DELTA_TABLE = (/                  &
    3.6850E-6, 1.0660E-6, 3.1593E-7, 1.0611E-7, 4.3892E-8,  &
    2.3095E-8, 1.3997E-8, 9.2345E-9, 6.4440E-9, 4.4828E-9,  &
    2.9997E-9, 2.0933E-9, 1.5072E-9, 1.1118E-9, 8.3628E-10, &
    4.9494E-10, 1.9634E-10, 8.6557E-11, 1.4328E-11, 2.9840E-12, &
    8.1056E-13, 3.1491E-13, 1.6813E-13, 1.0731E-13, 7.4155E-14 /)
    

! density table  (SIGMA = RHO/RHO0)
  REAL,PARAMETER,DIMENSION(SIZE(Z)):: SIGMA_TABLE = (/                              &
    5.680E-6, 1.632E-6, 4.575E-7, 1.341E-7, 4.061E-8, &
    1.614E-8, 7.932E-9, 4.461E-9, 2.741E-9, 1.694E-9, &
    1.007E-9, 6.380E-10, 4.240E-10, 2.923E-10, 2.074E-10, &
    1.116E-10, 3.871E-11, 1.564E-11, 2.288E-12, 4.257E-13, &
    9.279E-14, 2.506E-14, 9.272E-15, 4.701E-15, 2.907E-15 /)
    
  REAL,PARAMETER,DIMENSION(SIZE(Z)):: LOGDELTA = LOG(DELTA_TABLE)
  REAL,PARAMETER,DIMENSION(SIZE(Z)):: LOGSIGMA = LOG(SIGMA_TABLE)

  REAL,PARAMETER,DIMENSION(SIZE(Z)):: DLOGDELTA = (/   & 
    -0.174061, -0.177924, -0.167029, -0.142755, -0.107859,   &
    -0.079322, -0.064664, -0.054879, -0.048260, -0.042767,   &
    -0.037854, -0.034270, -0.031543, -0.029384, -0.027632,   &
    -0.024980, -0.021559, -0.019557, -0.016735, -0.014530,   &
    -0.011314, -0.007677, -0.005169, -0.003944, -0.003612 /)
  REAL,PARAMETER,DIMENSION(SIZE(Z)):: DLOGSIGMA = (/   & 
    -0.172421, -0.182258, -0.178090, -0.176372, -0.154322,   &
    -0.113750, -0.090582, -0.075033, -0.064679, -0.056067,   &
    -0.048461, -0.043042, -0.038869, -0.035648, -0.033063,   &
    -0.029164, -0.024220, -0.021336, -0.017686, -0.016035,   &
    -0.014327, -0.011631, -0.008248, -0.005580, -0.004227 /)

  INTEGER:: i,j,k                                                     ! counters
!-------------------------------------------------------------------------------

  IF (alt > Z(SIZE(Z))) THEN          ! trap altitudes greater than 1000 km.
    delta=DELTA_TABLE(SIZE(Z))        ! return value for Z=1000.0
    sigma=SIGMA_TABLE(SIZE(Z))        !   ditto 
    theta=1000.0/TZERO
    RETURN
  END IF

  i=1 
  j=SIZE(Z)                                    ! setting up for binary search
  DO
    k=(i+j)/2                                              ! integer division
    IF (alt < Z(k)) THEN
      j=k
    ELSE
      i=k
    END IF   
    IF (j <= i+1) EXIT
  END DO

  delta=EXP(EvaluateCubic(Z(i), LOGDELTA(i),   DLOGDELTA(i),          &
                        Z(i+1), LOGDELTA(i+1), DLOGDELTA(i+1), alt))

  sigma=EXP(EvaluateCubic(Z(i), LOGSIGMA(i),   DLOGSIGMA(i),          &
                        Z(i+1), LOGSIGMA(i+1), DLOGSIGMA(i+1), alt))

  theta=KineticTemperature(alt)/TZERO

  RETURN
END Subroutine UpperAtmosphere   ! ------------------------------------------

!+
SUBROUTINE LowerAtmosphere(alt, sigma, delta, theta)
!   -------------------------------------------------------------------------
! PURPOSE - Compute the properties of the 1976 standard atmosphere to 86 km.

  REAL,INTENT(IN)::  alt    ! geometric altitude, km.                        
  REAL,INTENT(OUT):: sigma  ! density/sea-level standard density              
  REAL,INTENT(OUT):: delta  ! pressure/sea-level standard pressure           
  REAL,INTENT(OUT):: theta  ! temperature/sea-level standard temperature

  INTEGER,PARAMETER:: NTAB=8       ! number of entries in the defining tables
! NOTE - REARTH and GMR are module parameters (constants).  
! If you copy this subroutine for use elsewhere, then uncomment following:
!   REAL,PARAMETER:: REARTH = 6356.766, GMR = 34.163195

  INTEGER:: i,j,k                                                  ! counters
  REAL:: h                                       ! geopotential altitude (km)
  REAL:: tgrad, tbase      ! temperature gradient and base temp of this layer
  REAL:: tlocal                                           ! local temperature
  REAL:: deltah                             ! height above base of this layer

!===============================================================================
!      C O N S T A N T   A R R A Y S  O F  1 9 7 6   S T D.  A T M O S P H E R E
!===============================================================================
  REAL,DIMENSION(NTAB),PARAMETER:: HTAB = &
                          (/0.0, 11.0, 20.0, 32.0, 47.0, 51.0, 71.0, 84.852/)
  REAL,DIMENSION(NTAB),PARAMETER:: TTAB = &
          (/288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65, 186.946/)
  REAL,DIMENSION(NTAB),PARAMETER:: PTAB = &
               (/1.0, 2.233611E-1, 5.403295E-2, 8.5666784E-3, 1.0945601E-3, &
                                     6.6063531E-4, 3.9046834E-5, 3.68501E-6/)
  REAL,DIMENSION(NTAB),PARAMETER:: GTAB = &
                                (/-6.5, 0.0, 1.0, 2.8, 0.0, -2.8, -2.0, 0.0/)
!----------------------------------------------------------------------------
  h=alt*REARTH/(alt+REARTH)      ! convert geometric to geopotential altitude

  i=1 
  j=NTAB                                          ! setting up for binary search
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

  sigma=delta/theta                                           ! density ratio
  RETURN
END Subroutine LowerAtmosphere   ! ------------------------------------------

!+
SUBROUTINE SimpleAtmosphere(alt,sigma,delta,theta)
!   -------------------------------------------------------------------------
! PURPOSE - Compute the characteristics of the atmosphere below 20 km.

! NOTES-Correct to 20 km. Only approximate above there

  IMPLICIT NONE

  REAL,INTENT(IN)::  alt    ! geometric altitude, km.
  REAL,INTENT(OUT):: sigma  ! density/sea-level standard density             
  REAL,INTENT(OUT):: delta  ! pressure/sea-level standard pressure            
  REAL,INTENT(OUT):: theta  ! temperature/sea-level standard temperature   

! NOTE - REARTH and GMR are module parameters (constants).  
! If you copy this subroutine for use elsewhere, then uncomment following:
!   REAL,PARAMETER:: REARTH = 6356.766, GMR = 34.163195

  REAL:: h   ! geopotential altitude
!----------------------------------------------------------------------------
  h=alt*REARTH/(alt+REARTH)      ! convert geometric to geopotential altitude

  IF (h < 11.0) THEN
    theta=1.0+(-6.5/288.15)*h                                   ! Troposphere
    delta=theta**(GMR/6.5)
  ELSE
    theta=216.65/288.15                                        ! Stratosphere
    delta=0.2233611*EXP(-GMR*(h-11.0)/216.65)
  END IF

  sigma=delta/theta
  RETURN
END Subroutine SimpleAtmosphere   ! --------------------------------------------

!+
PURE FUNCTION Viscosity(theta) RESULT(visc)
!   ----------------------------------------------------------------------------
! PURPOSE - Compute viscosity using Sutherland's formula.
!        Returns viscosity in kg/(meter-sec)

  REAL,INTENT(IN) :: theta                ! temperature/sea-level temperature  
  REAL:: visc
  REAL:: temp                              ! temperature in kelvins
! NOTE - BETAVISC and TZERO are module parameters (constants).  
! If you copy this subroutine for use elsewhere, uncomment the following:
!  REAL,PARAMETER:: BETAVISC = 1.458E-6, SUTH = 110.4    
!-------------------------------------------------------------------------------
  temp=TZERO*theta
  visc=BETAVISC*Sqrt(temp*temp*temp)/(temp+SUTH)
  RETURN
END Function Viscosity   ! -----------------------------------------------------

!+
SUBROUTINE Atmosphere(alt,sigma,delta,theta)
!   ----------------------------------------------------------------------------
! PURPOSE - Compute the characteristics of the U.S. Standard Atmosphere 1976

  IMPLICIT NONE
  REAL,INTENT(IN)::  alt    ! geometric altitude, km.
  REAL,INTENT(OUT):: sigma  ! density/sea-level standard density             
  REAL,INTENT(OUT):: delta  ! pressure/sea-level standard pressure            
  REAL,INTENT(OUT):: theta  ! temperature/sea-level standard temperature   
!-------------------------------------------------------------------------------
  IF (alt > 86.0) THEN
    CALL UpperAtmosphere(alt,sigma,delta,theta)
  ELSE
    CALL LowerAtmosphere(alt,sigma,delta,theta)
  END IF
  RETURN
END Subroutine Atmosphere   ! --------------------------------------------------

!+
SUBROUTINE TransportRatios(alt, gRatio,dynamicViscosityRatio, &
  kinematicViscosityRatio, numberDensityRatio, collisionFrequencyRatio, &
  particleSpeedRatio, &
  meanFreePathRatio,thermalConductivityRatio, soundSpeedRatio, &
  molecularWeightRatio )
! ------------------------------------------------------------------------------
! PURPOSE - Compute various non-dimensional ratios.
  REAL,INTENT(IN):: alt   ! geometric altitude in kilometers
  REAL,INTENT(OUT):: gRatio 
  REAL,INTENT(OUT):: dynamicViscosityRatio
  REAL,INTENT(OUT):: kinematicViscosityRatio
  REAL,INTENT(OUT):: numberDensityRatio
  REAL,INTENT(OUT):: collisionFrequencyRatio
  REAL,INTENT(OUT):: particleSpeedRatio 
  REAL,INTENT(OUT):: meanFreePathRatio
  REAL,INTENT(OUT):: thermalConductivityRatio 
  REAL,INTENT(OUT):: soundSpeedRatio 
  REAL,INTENT(OUT):: molecularWeightRatio
  
  REAL:: sigma,theta,delta 
  REAL:: kappa
  REAL:: molwt
!-------------------------------------------------------------------------------
  CALL Atmosphere(alt,sigma,delta,theta)
  molwt=MolecularWeight(alt)
  molecularWeightRatio=molwt/MOLWT_ZERO
  
  gRatio=(REARTH/(alt+REARTH))**2
  dynamicViscosityRatio=Viscosity(theta)/MUZERO 
  kinematicViscosityRatio=dynamicViscosityRatio/sigma
  numberDensityRatio=sigma*MOLWT_ZERO/molwt
  particleSpeedRatio=SQRT(theta/molecularWeightRatio)
  meanFreePathRatio=molecularWeightRatio/sigma
  collisionFrequencyRatio=particleSpeedRatio/meanFreePathRatio
  kappa=ThermalConductivity(theta*TZERO)
  thermalConductivityRatio=kappa/KAPPA_ZERO  
  soundSpeedRatio=SQRT(theta)

   
  RETURN 
END Subroutine TransportRatios   ! ---------------------------------------------

!+
PURE FUNCTION MolecularWeight(altKm) RESULT(mw)
! ------------------------------------------------------------------------------
! PURPOSE - Compute molecular weight of air at altitude by interpolation
!    in tables.
  REAL,INTENT(IN):: altKm
  INTEGER,PARAMETER:: NTABLE = 25

! tables of altitude (km), mol.wt., and d(mol.wt.)/d(altKm)
  REAL,PARAMETER,DIMENSION(NTABLE):: Z = (/      &
     86.0,  93.0, 100.0, 107.0, 114.0, &
    121.0, 128.0, 135.0, 142.0, 150.0, &
    160.0, 170.0, 180.0, 190.0, 200.0, &
    220.0, 260.0, 300.0, 400.0, 500.0, &
    600.0, 700.0, 800.0, 900.0,1000.0 /)
    
  REAL,PARAMETER,DIMENSION(NTABLE):: M = (/    &
    28.95, 28.82, 28.40, 27.64, 26.79,   &
    26.12, 25.58, 25.09, 24.62, 24.10,   &
    23.49, 22.90, 22.34, 21.81, 21.30,   &
    20.37, 18.85, 17.73, 15.98, 14.33,   &
    11.51,  8.00,  5.54,  4.40,  3.94 /)  

  REAL,PARAMETER,DIMENSION(NTABLE):: MP = (/                 &
    -0.001340, -0.036993, -0.086401, -0.123115, -0.111136,   &
    -0.083767, -0.072368, -0.068190, -0.066300, -0.063131,   &
    -0.059786, -0.057724, -0.054318, -0.052004, -0.049665,   &
    -0.043499, -0.032674, -0.023804, -0.014188, -0.021444,   &
    -0.034136, -0.031911, -0.017321, -0.006804, -0.003463 /)


  
  REAL:: mw
  INTEGER:: i,j,k
!-------------------------------------------------------------------------------
  IF (altKm <= 86.0) THEN 
    mw = MOLWT_ZERO
  ELSE 
    i=1 
    j=SIZE(Z)                                    ! setting up for binary search
    DO
      k=(i+j)/2                                              ! integer division
      IF (altKm < Z(k)) THEN
        j=k
      ELSE
        i=k
      END IF   
      IF (j <= i+1) EXIT
    END DO
    mw = EvaluateCubic(Z(i), M(i), MP(i), Z(i+1), M(i+1), MP(i+1), altKm)
  END IF
  
  RETURN 
END Function MolecularWeight   ! -----------------------------------------------      
   
END Module Atmosphere1976   ! ==================================================



!+
MODULE BigTablesProcedures
! ------------------------------------------------------------------------------
! PURPOSE - Write the files tables1.html, tables2.html, etc. corresponding to 
!  the tables in "U.S. Standard Atmosphere, 1976"
USE Atmosphere1976
IMPLICIT NONE

  INTEGER,PARAMETER:: SP = SELECTED_REAL_KIND(6,30)
  INTEGER,PARAMETER:: DP = SELECTED_REAL_KIND(14,80)
  INTEGER,PARAMETER:: WP = DP   ! working precision

  INTEGER,PRIVATE,PARAMETER:: HTML=1

  CHARACTER(LEN=*),PARAMETER:: TD1="<td>", TD2="</td>"
  CHARACTER(LEN=*),PARAMETER:: TDC1='<td class="center">'
  CHARACTER(LEN=*),PARAMETER:: BLANKCELL = '<td></td>'
  
  CHARACTER(LEN=*),PARAMETER:: FMT1 = '(3A)'
  CHARACTER(LEN=*),PARAMETER:: PDAS = &
     "Public Domain Aeronautical Software (PDAS) &nbsp;"
  
  INTEGER,PARAMETER:: USMIN=0, USMAX=3000, USDEL=10
  

!-------------------------------------------------------------------------------

CONTAINS


!+
FUNCTION GetDateStr() RESULT (s)
!   ---------------------------------------------------------------------------
! PURPOSE - Return a string with the current date
  IMPLICIT NONE
  CHARACTER(LEN=*),PARAMETER:: MONTH='JanFebMarAprMayJunJulAugSepOctNovDec'
  CHARACTER(LEN=*),PARAMETER:: FMT='(I2,1X,A3,I5)'
  CHARACTER(LEN=11):: s
  INTEGER,DIMENSION(8):: v
!------------------------------------------------------------------------------
  CALL DATE_AND_TIME(VALUES=v)
  WRITE(s,FMT) v(3), MONTH(3*v(2)-2:3*v(2)), v(1)
  RETURN
END FUNCTION GetDateStr   ! ----------------------------------------------------





!+
SUBROUTINE MakeMainPage()
! ------------------------------------------------------------------------------
! PURPOSE - Open the output file bt3.htm and write the HTML header, the CSS
!  style sheets, and the javascript functions.
!  The <head> element contains the <meta charset=utf-8 />, <title>, <style> 
!  and <script> elements.

!-------------------------------------------------------------------------------
  OPEN(UNIT=HTML, FILE='bigtablesf.html', STATUS='REPLACE', ACTION='WRITE')
  WRITE(HTML,FMT1) '<!DOCTYPE html>'
  WRITE(HTML,FMT1) '<html lang="en">'
  WRITE(HTML,FMT1) '<head>'
  WRITE(HTML,FMT1) '<meta charset="utf-8" />'
  WRITE(HTML,FMT1) '<title>Table of the  U.S. Standard Atmosphere 1976</title>'
  
  WRITE(HTML,FMT1) '<style>'
  WRITE(HTML,FMT1) 'body {font-size: 99%;}'
  WRITE(HTML,FMT1) 'table {border: none;}'
  WRITE(HTML,FMT1) 'td {border: none; padding: 0 0 0.2em 0.2em; text-align: right;}'
  WRITE(HTML,FMT1) 'td.center {text-align:center}'
  WRITE(HTML,FMT1) 'td.right {text-align:right}'
  WRITE(HTML,FMT1) 'tr.gray {background-color: #f0f0f0;}'
  
  WRITE(HTML,FMT1) 'header {width: 100%; border-bottom: 3px solid blue; padding: 3px; }'
  WRITE(HTML,FMT1) 'footer {width: 100%; font-size: 80%; '
  WRITE(HTML,FMT1) '   padding: 3px; border-top: 3px solid red; }'
  WRITE(HTML,FMT1) 'h1 { font-family: Helvetica,Arial,sans-serif;'
  WRITE(HTML,FMT1) '   text-align: center; font-weight: bold; font-size: 2.0em;'
  WRITE(HTML,FMT1) '   color: blue; margin-top: 0.5em; }'

!...Define the properties of the seven buttons
  WRITE(HTML,FMT1) '#MainButton   {margin: 10px; color: black;'
  WRITE(HTML,FMT1) 'background-color: white; font: 1em Arial,sans-serif;}'
  WRITE(HTML,FMT1) '#Table1Button {margin: 10px; color: black;'
  WRITE(HTML,FMT1) 'background-color: red; font: 1em Arial,sans-serif;}'
  WRITE(HTML,FMT1) '#Table2Button {margin: 10px; color: black;'
  WRITE(HTML,FMT1) 'background-color: orange; font: 1em Arial,sans-serif;}'
  WRITE(HTML,FMT1) '#Table3Button {margin: 10px; color: black;'
  WRITE(HTML,FMT1) 'background-color: yellow; font: 1em Arial,sans-serif;}'
  WRITE(HTML,FMT1) '#Table4Button {margin: 10px; color: white;'
  WRITE(HTML,FMT1) 'background-color: green; font: 1em Arial,sans-serif;}'
  WRITE(HTML,FMT1) '#Table5Button {margin: 10px; color: white;'
  WRITE(HTML,FMT1) 'background-color: blue; font: 1em Arial,sans-serif;}'
  WRITE(HTML,FMT1) '#Table6Button {margin: 10px; color: white;'
  WRITE(HTML,FMT1) 'background-color: violet; font: 1em Arial,sans-serif;}'
  
!... Begin with all 6 tables turned off  
  WRITE(HTML,FMT1) '#table1 {display:none;}'
  WRITE(HTML,FMT1) '#table2 {display:none;}'
  WRITE(HTML,FMT1) '#table3 {display:none;}'
  WRITE(HTML,FMT1) '#table4 {display:none;}'
  WRITE(HTML,FMT1) '#table5 {display:none;}'
  WRITE(HTML,FMT1) '#table6 {display:none;}'
    
  WRITE(HTML,FMT1) 'div.banner {width: 100%; background: blue; color: white;'
  WRITE(HTML,FMT1) '  font-size: 120%; font-family: Verdana,Arial,Helvetica,'
  WRITE(HTML,FMT1) '  sans-serif; text-align: right;}'

  WRITE(HTML,FMT1) '</style>'

  WRITE(HTML,FMT1) '<script>' 

  WRITE(HTML,FMT1) 'function HideAll() {' 
  WRITE(HTML,FMT1) '  document.getElementById("main").style.display="none";'
  WRITE(HTML,FMT1) '  document.getElementById("table1").style.display="none";'
  WRITE(HTML,FMT1) '  document.getElementById("table2").style.display="none";' 
  WRITE(HTML,FMT1) '  document.getElementById("table3").style.display="none";' 
  WRITE(HTML,FMT1) '  document.getElementById("table4").style.display="none";' 
  WRITE(HTML,FMT1) '  document.getElementById("table5").style.display="none";' 
  WRITE(HTML,FMT1) '  document.getElementById("table6").style.display="none";'   
  WRITE(HTML,FMT1) '}   // --- End of function HideAll' 

  WRITE(HTML,FMT1) 'function ShowMain() {' 
  WRITE(HTML,FMT1) '  HideAll();'
  WRITE(HTML,FMT1) '  document.getElementById("main").style.display="block";'
  WRITE(HTML,FMT1) '}   // --- End of function ShowMain' 

  WRITE(HTML,FMT1) 'function ShowTable1() {' 
  WRITE(HTML,FMT1) '  HideAll();'
  WRITE(HTML,FMT1) '  document.getElementById("table1").style.display="block";' 
  WRITE(HTML,FMT1) '}   // --- End of function ShowTable1' 

  WRITE(HTML,FMT1) 'function ShowTable2() {' 
  WRITE(HTML,FMT1) '  HideAll();' 
  WRITE(HTML,FMT1) '  document.getElementById("table2").style.display="block";' 
  WRITE(HTML,FMT1) '}   // --- End of function ShowTable2' 

  WRITE(HTML,FMT1) 'function ShowTable3() {' 
  WRITE(HTML,FMT1) '  HideAll();' 
  WRITE(HTML,FMT1) '  document.getElementById("table3").style.display="block";' 
  WRITE(HTML,FMT1) '}   // --- End of function ShowTable3' 

  WRITE(HTML,FMT1) 'function ShowTable4() {' 
  WRITE(HTML,FMT1) '  HideAll();'  
  WRITE(HTML,FMT1) '  document.getElementById("table4").style.display="block";' 
  WRITE(HTML,FMT1) '}   // --- End of function ShowTable4' 

  WRITE(HTML,FMT1) 'function ShowTable5() {' 
  WRITE(HTML,FMT1) '  HideAll();'
  WRITE(HTML,FMT1) '  document.getElementById("table5").style.display="block";' 
  WRITE(HTML,FMT1) '}   // --- End of function ShowTable5' 

  WRITE(HTML,FMT1) 'function ShowTable6() {' 
  WRITE(HTML,FMT1) '  HideAll();'
  WRITE(HTML,FMT1) '  document.getElementById("table6").style.display="block";' 
  WRITE(HTML,FMT1) '}   // --- End of function ShowTable6' 

  WRITE(HTML,FMT1) '</script>' 

  WRITE(HTML,FMT1) '</head>' 
  WRITE(HTML,FMT1) '<body>'
 
  WRITE(HTML,FMT1) '<a href="index.html">PDAS home</a> &gt; ', &
                   '<a href="atmos.html">Standard Atmosphere</a> &gt; ', &
                    'BigTables'  
 
  WRITE(HTML,FMT1) '<div class="banner">', PDAS, '</div>'
  WRITE(HTML,FMT1) '<header><h1>Tables of the U.S. Standard Atmosphere, 1976</h1></header>'

  WRITE(HTML,FMT1) '<!-- This is an example of a form that is never submitted.'
  WRITE(HTML,FMT1) '  All of the computation is done by the included'
  WRITE(HTML,FMT1) '  JavaScript functions. No server calculations needed. -->'
  WRITE(HTML,FMT1) '<form>'
  WRITE(HTML,FMT1) '<button id="MainButton" type="button"', &
     ' onclick="ShowMain();">Main</button>'
  WRITE(HTML,FMT1) '<button id="Table1Button" type="button"', & 
     ' onclick="ShowTable1();">Table 1</button>'
  WRITE(HTML,FMT1) '<button id="Table2Button" type="button"', & 
     ' onclick="ShowTable2();">Table 2</button>'
  WRITE(HTML,FMT1) '<button id="Table3Button" type="button"', & 
     ' onclick="ShowTable3();">Table 3</button>'
  WRITE(HTML,FMT1) '<button id="Table4Button" type="button"', & 
     ' onclick="ShowTable4();">Table 4</button>'
  WRITE(HTML,FMT1) '<button id="Table5Button" type="button"', & 
     ' onclick="ShowTable5();">Table 5</button>'
  WRITE(HTML,FMT1) '<button id="Table6Button" type="button"', & 
     ' onclick="ShowTable6();">Table 6</button>'
  WRITE(HTML,FMT1) '</form>'

  WRITE(HTML,FMT1) '<div id="main">'
  WRITE(HTML,FMT1) '<p>Make tables of atmospheric properties'
  WRITE(HTML,FMT1) 'from 0 to 1000 km.</p>'
  WRITE(HTML,FMT1) '<p>The tables are based on those from the reference'
  WRITE(HTML,FMT1) 'document, Part 4.</p>'
  WRITE(HTML,FMT1) '<ol>'
  WRITE(HTML,FMT1) '<li>Geopotential altitude, temperature, temperature ratio,'
  WRITE(HTML,FMT1) 'pressure, pressure ratio, density, density ratio,'
  WRITE(HTML,FMT1) 'sound speed, garvity for geometric altitudes in kilometers.</li>'

  WRITE(HTML,FMT1) '<li>Geopotential altitude, dynamic viscosity, kinematic viscosity,'
  WRITE(HTML,FMT1) 'pressure scale height, number density, mean particle speed,'
  WRITE(HTML,FMT1) 'mean collision frequency, mean free path,'
  WRITE(HTML,FMT1) 'thermal conductivity coefficient and molecular weight for'
  WRITE(HTML,FMT1) 'geometric altitudes in kilometers.</li>'
  
  WRITE(HTML,FMT1) '<li>Non-dimensional ratios of dynamic viscosity,'
  WRITE(HTML,FMT1) 'kinematic viscosity, number density, particle speed,'
  WRITE(HTML,FMT1) 'mean collision frequency, mean free path, thermal conductivity,'
  WRITE(HTML,FMT1) 'sound speed and molecular weight for altitudes in kilometers.'
  WRITE(HTML,FMT1) 'All quantities are referenced to sea-level values.</li>'
  
  WRITE(HTML,FMT1) '<li>Geopotential altitude, temperature, temperature ratio,'
  WRITE(HTML,FMT1) 'pressure, pressure ratio, mass density, density ratio,'
  WRITE(HTML,FMT1) 'weight density, sound speed and gravity for geometric altitudes'
  WRITE(HTML,FMT1) 'in thousand feet. Table entries in US customary units.</li>'
  
  WRITE(HTML,FMT1) '<li>Geopotential altitude, dynamic viscosity, kinematic viscosity,'
  WRITE(HTML,FMT1) 'pressure scale height, number density, mean particle speed,'
  WRITE(HTML,FMT1) 'mean collision frequency, mean free path,'
  WRITE(HTML,FMT1) 'thermal conductivity coefficient and molecular weight for'
  WRITE(HTML,FMT1) 'geometric altitudes in thousand feet.'
  WRITE(HTML,FMT1) 'Table entries in US customary units.</li>'

  WRITE(HTML,FMT1) '<li>Non-dimensional ratios of dynamic viscosity,'
  WRITE(HTML,FMT1) 'kinematic viscosity, number density, particle speed,'
  WRITE(HTML,FMT1) 'mean collision frequency, mean free path, thermal conductivity,'
  WRITE(HTML,FMT1) 'sound speed and molecular weight for altitudes in thousand feet.'
  WRITE(HTML,FMT1) 'All quantities are referenced to sea-level values.</li>'

  WRITE(HTML,FMT1) '</ol>'
  WRITE(HTML,FMT1) '<p>To return to this page, press Main</p>'
  WRITE(HTML,FMT1) '</div>'
  RETURN
END Subroutine MakeMainPage   ! ------------------------------------------------  
!+
SUBROUTINE MakeTable1()
! ------------------------------------------------------------------------------
! PURPOSE - Create Table 1
  
  INTEGER:: i                       ! geometric altitude in kilometers (INTEGER)
  REAL:: altKm                         ! geometric altitude in kilometers (REAL)
  REAL:: hkm                               ! geopotential altitude in kilometers
  REAL:: sigma,delta,theta, density
  REAL:: asound, g
!-------------------------------------------------------------------------------
  
 
  WRITE(HTML,FMT1) '<div id="table1">'
  WRITE(HTML,FMT1) "<h2>Table 1 Temperature, Pressure, & Density (SI units)</h2>"
  WRITE(HTML,FMT1) "<h3>0 to 1000 Km in steps of 5 Km</h3>"
  WRITE(HTML,FMT1) "<table>"
 
  DO i=0,1000,5
    IF (MODULO(i,200) == 0) THEN  ! add title rows every 40 rows
      WRITE(HTML,FMT1) "<tr>"     
      WRITE(HTML,FMT1) TDC1,'Z',TD2, TDC1,'H',TD2
      WRITE(HTML,FMT1) TDC1,'T',TD2, TDC1,'T/T<sub>0</sub>',TD2
      WRITE(HTML,FMT1) TDC1,'p',TD2, TDC1,'p/p<sub>0</sub>',TD2
      WRITE(HTML,FMT1) TDC1,'&rho;',TD2, TDC1,'&rho;/&rho;<sub>0</sub>',TD2
      WRITE(HTML,FMT1) TDC1,'c',TD2, TDC1,'g',TD2
      WRITE(HTML,FMT1) "</tr>"

      IF (i==0) THEN 
        WRITE(HTML,FMT1) "<tr>"   ! only printed at top of table
        WRITE(HTML,FMT1) TDC1,'km',TD2, TDC1,'km',TD2, TDC1,'K',TD2, BLANKCELL 
        WRITE(HTML,FMT1) TDC1,'Pa',TD2, BLANKCELL
        WRITE(HTML,FMT1) TDC1,'kg m<sup>-3</sup>',TD2, BLANKCELL
        WRITE(HTML,FMT1) TDC1,'m/s',TD2, TDC1, 'm s<sup>-2</sup>',TD2
        WRITE(HTML,FMT1) "</tr>"
      END IF 
    END IF
    IF (MODULO(i,10)==0) THEN
      WRITE(HTML,FMT1) '<tr class="gray">'   ! make every other row gray
    ELSE
      WRITE(HTML,FMT1) '<tr>'
    END IF

    altKm=REAL(i)
    hkm=altKm*REARTH/(altKm+REARTH)      ! convert geometric to geopotential
    g=GZERO*(REARTH/(altKm+REARTH))**2
    
    CALL Atmosphere(altKm, sigma,delta,theta)
    density=sigma*RHOZERO
    asound=ASOUNDZERO*SQRT(theta)

    WRITE(HTML,'(A,I0,A)')     TD1, INT(altKm), TD2
    WRITE(HTML,'(A,F0.1,A)')   TD1, hkm, TD2
    WRITE(HTML,'(A,F0.3,A)')   TD1,theta*TZERO,TD2
    WRITE(HTML,'(A,F6.4,A)')   TD1,theta,TD2
    WRITE(HTML,'(A,ES12.4,A)') TD1,delta*PZERO,TD2, TD1,delta,TD2
    WRITE(HTML,'(A,ES12.4,A)') TD1,sigma*RHOZERO,TD2, TD1,sigma,TD2
    WRITE(HTML,'(A,F0.2,A)')   TD1,asound,TD2
    WRITE(HTML,'(A,F0.4,A)')   TD1,g,TD2

    WRITE(HTML,FMT1) "</tr>"

  END DO
  WRITE(HTML,'(A)') '</table>'
  WRITE(HTML,'(A)') '</div>'
  
  RETURN
END Subroutine MakeTable1   ! ---------------------------------------------------


!+
SUBROUTINE MakeTable2()
! ------------------------------------------------------------------------------
! PURPOSE - Create Table 2
  
  INTEGER:: i
  REAL:: altKm                                          ! altitude in kilometers
  REAL:: hkm                               ! geopotential altitude in kilometers
  REAL:: sigma,delta,theta
!  REAL:: temp,pressure,density, asound
  REAL:: dynamicViscosity,kinematicViscosity
  REAL:: kappa
  REAL:: collisionFrequency
  REAL:: g,gRatio
  REAL:: meanFreePath
!  REAL:: mw              ! molecular weight
  REAL:: numDensity
  REAL:: partSpeed
  REAL:: pressureScaleHeight
  
  REAL:: dynamicViscosityRatio
  REAL:: kinematicViscosityRatio, numberDensityRatio, collisionFrequencyRatio
  REAL:: particleSpeedRatio
  REAL:: meanFreePathRatio,thermalConductivityRatio, soundSpeedRatio
  REAL:: molecularWeightRatio


!-------------------------------------------------------------------------------
  
  WRITE(HTML,FMT1) '<div id="table2">'
  WRITE(HTML,FMT1) "<h2>Table 2 - Transport Properties (SI units)</h2>"
  WRITE(HTML,FMT1) "<h3>0 to 1000 Km in steps of 5 km</h3>"
  WRITE(HTML,FMT1) "<table>"

  WRITE(HTML,FMT1) "<tr>"   ! header row, only printed at top of table
  WRITE(HTML,FMT1) TDC1,'alt',TD2, TDC1,'H',TD2
  WRITE(HTML,FMT1) TDC1,'Dynamic<br />viscosity',TD2
  WRITE(HTML,FMT1) TDC1,'Kinematic<br />viscosity',TD2
  WRITE(HTML,FMT1) TDC1,'Pressure<br />scale<br />height',TD2
  WRITE(HTML,FMT1) TDC1,'Number<br />density',TD2
  WRITE(HTML,FMT1) TDC1,'Particle<br />speed',TD2
  WRITE(HTML,FMT1) TDC1,'Collision<br />frequency',TD2
  WRITE(HTML,FMT1) TDC1,'Mean<br />free<br />path',TD2
  WRITE(HTML,FMT1) TDC1,'Thermal<br />conductivity',TD2
  WRITE(HTML,FMT1) TDC1,'Mol.<br />weight',TD2 
  WRITE(HTML,FMT1) "</tr>"

  DO i=0,1000,5
    IF (MODULO(i,200) == 0) THEN  ! add symbol row every 40 rows
      WRITE(HTML,FMT1) "<tr>"   ! symbol row
      WRITE(HTML,FMT1) TDC1,'Z',TD2, TDC1,'H',TD2
      WRITE(HTML,FMT1) TDC1,'&mu;',TD2
      WRITE(HTML,FMT1) TDC1,'&eta;',TD2
      WRITE(HTML,FMT1) TDC1,'H<sub>p</sub>',TD2
      WRITE(HTML,FMT1) TDC1,'n',TD2
      WRITE(HTML,FMT1) TDC1,'V',TD2
      WRITE(HTML,FMT1) TDC1,'&nu;',TD2
      WRITE(HTML,FMT1) TDC1,'L',TD2
      WRITE(HTML,FMT1) TDC1,'&kappa;',TD2
      WRITE(HTML,FMT1) TDC1,'M',TD2
      WRITE(HTML,FMT1) "</tr>"

      IF (i == 0) THEN 
        WRITE(HTML,FMT1) "<tr>"   ! units row, only printed at top
        WRITE(HTML,FMT1) TDC1,'km',TD2, TDC1,'km',TD2
        WRITE(HTML,FMT1) TDC1,'kg m<sup>-1</sup> s<sup>-1</sup>',TD2
        WRITE(HTML,FMT1) TDC1,'m<sup>2</sup>/s',TD2
        WRITE(HTML,FMT1) TDC1,'m',TD2, TDC1,'m<sup>-3</sup>',TD2
        WRITE(HTML,FMT1) TDC1,'m/s',TD2
        WRITE(HTML,FMT1) TDC1,'s<sup>-1</sup>',TD2
        WRITE(HTML,FMT1) TDC1,'m',TD2
        WRITE(HTML,FMT1) TDC1,'W/(m K)',TD2, BLANKCELL
        WRITE(HTML,FMT1) "</tr>"
      END IF   
    END IF   
    
    IF (MODULO(i,10)==0) THEN
      WRITE(HTML,FMT1) '<tr class="gray">'
    ELSE
      WRITE(HTML,FMT1) '<tr>'
    END IF

    altKm=REAL(i)
    hkm=altkm*REARTH/(altkm+REARTH)      ! convert geometric to geopotential
    
    gRatio=(REARTH/(altKm+REARTH))**2
    g=gRatio*GZERO
    CALL Atmosphere(altKm, sigma,delta,theta)
    CALL TransportRatios(altkm, gRatio,dynamicViscosityRatio, &
      kinematicViscosityRatio, numberDensityRatio, collisionFrequencyRatio, &
      particleSpeedRatio, meanFreePathRatio,thermalConductivityRatio, &
      soundSpeedRatio, molecularWeightRatio )

    dynamicViscosity=MUZERO*dynamicViscosityRatio
    kinematicViscosity=ETAZERO*kinematicViscosityRatio
    pressureScaleHeight=PRESSURE_SCALE_HEIGHT_ZERO*theta/(gRatio*molecularWeightRatio)
    numDensity=NUM_DENSITY_ZERO*numberDensityRatio
    partSpeed=PART_SPEED_ZERO*particleSpeedRatio 
    meanFreePath=FREE_PATH_ZERO*meanFreePathRatio 
    collisionFrequency=partSpeed/meanFreePath
    kappa=KAPPA_ZERO*thermalConductivityRatio 

    WRITE(HTML,'(A,I0,A)')     TD1, INT(altkm), TD2
    WRITE(HTML,'(A,F0.1,A)')   TD1, hkm, TD2 
    WRITE(HTML,'(A,ES12.4,A)') TD1,dynamicViscosity,TD2
    WRITE(HTML,'(A,ES12.4,A)') TD1,kinematicViscosity,TD2
    WRITE(HTML,'(A,I0,A)')     TD1,NINT(pressureScaleHeight),TD2
    WRITE(HTML,'(A,ES12.4,A)') TD1,numDensity,TD2
    WRITE(HTML,'(A,F0.2,A)')   TD1,partSpeed,TD2  
    WRITE(HTML,'(A,ES12.4,A)') TD1,collisionFrequency,TD2
    WRITE(HTML,'(A,ES12.4,A)') TD1,meanFreePath,TD2
    WRITE(HTML,'(A,ES12.4,A)') TD1,kappa, TD2
    WRITE(HTML,'(A,F0.3,A)')   TD1,molecularWeightRatio*MOLWT_ZERO, TD2 
    WRITE(HTML,FMT1) "</tr>"

  END DO
  WRITE(HTML,'(A)') '</table>','</div>'

  RETURN
END Subroutine MakeTable2   ! ---------------------------------------------------

!+
SUBROUTINE MakeTable3()
! ------------------------------------------------------------------------------  
  INTEGER:: i                        ! geometric altitude in kilometers(INTEGER)
  REAL:: altKm                                ! geometric altitude in kilometers
  REAL:: hkm                               ! geopotential altitude in kilometers
  REAL:: sigma,delta,theta
!  REAL:: kappa
  REAL:: gRatio,dynamicViscosityRatio
  REAL:: kinematicViscosityRatio, numberDensityRatio, collisionFrequencyRatio
  REAL:: particleSpeedRatio
  REAL:: meanFreePathRatio,thermalConductivityRatio, soundSpeedRatio
  REAL:: molecularWeightRatio


!-------------------------------------------------------------------------------
  
  WRITE(HTML,FMT1) '<div id="table3">'
  WRITE(HTML,FMT1) "<h2>Table 3 - Atmosphere Properties (non-dimensional)</h2>"
  WRITE(HTML,FMT1) "<h3>0 to 1000 Km in steps of 5 km</h3>"
  WRITE(HTML,FMT1) "<p>The units of alt and H are kilometers.</p>"
  WRITE(HTML,FMT1) "<table>"

  DO i=0,1000,5
    IF (MODULO(i,200)==0) THEN 
      WRITE(HTML,FMT1) "<tr>"     
      WRITE(HTML,FMT1) TDC1,'Z',TD2, TDC1,'H',TD2
      WRITE(HTML,FMT1) TDC1,'&mu;/&mu;<sub>0</sub>',TD2
      WRITE(HTML,FMT1) TDC1,'&eta;/&eta;<sub>0</sub>',TD2
      WRITE(HTML,FMT1) TDC1,'n/n<sub>0</sub>',TD2
      WRITE(HTML,FMT1) TDC1,'V/V<sub>0</sub>',TD2                      
      WRITE(HTML,FMT1) TDC1,'&nu;/&nu;<sub>0</sub>',TD2
      WRITE(HTML,FMT1) TDC1,'L/L<sub>0</sub>',TD2
      WRITE(HTML,FMT1) TDC1,'&kappa;/&kappa;<sub>0</sub>',TD2
      WRITE(HTML,FMT1) TDC1,'c/c<sub>0</sub>',TD2
      WRITE(HTML,FMT1) TDC1,'M/M<sub>0</sub>',TD2
      WRITE(HTML,FMT1) "</tr>"
    END IF 
    
    IF (MODULO(i,10)==0) THEN  ! Make every other row with light gray background.
      WRITE(HTML,FMT1) '<tr class="gray">'
    ELSE
      WRITE(HTML,FMT1) '<tr>'
    END IF

    altKm=REAL(i)
    hkm=altkm*REARTH/(altkm+REARTH)      ! convert geometric to geopotential
    
    CALL Atmosphere(altKm, sigma,delta,theta)  
    CALL TransportRatios(altkm, gRatio,dynamicViscosityRatio, &
      kinematicViscosityRatio, numberDensityRatio, collisionFrequencyRatio, &
      particleSpeedRatio, meanFreePathRatio,thermalConductivityRatio, &
      soundSpeedRatio, molecularWeightRatio )

    WRITE(HTML,'(A,I0,A)')     TD1, INT(altkm), TD2
    WRITE(HTML,'(A,F0.1,A)')   TD1, hkm, TD2    
    WRITE(HTML,'(A,F6.4,A)')   TD1, dynamicViscosityRatio, TD2
    WRITE(HTML,'(A,ES12.4,A)') TD1, kinematicViscosityRatio, TD2
    WRITE(HTML,'(A,ES12.4,A)') TD1, numberDensityRatio, TD2
    WRITE(HTML,'(A,F6.4,A)')   TD1, particleSpeedRatio, TD2
    WRITE(HTML,'(A,ES12.4,A)') TD1, collisionFrequencyRatio, TD2
    WRITE(HTML,'(A,ES12.4,A)') TD1, meanFreePathRatio, TD2
    WRITE(HTML,'(A,F6.4,A)')   TD1, thermalConductivityRatio, TD2
    WRITE(HTML,'(A,F6.4,A)')   TD1, soundSpeedRatio, TD2
    WRITE(HTML,'(A,F6.4,A)')   TD1, molecularWeightRatio, TD2
    WRITE(HTML,FMT1) "</tr>"

  END DO
  
  WRITE(HTML,'(A)') '</table>','</div>'

  RETURN
END Subroutine MakeTable3   ! ---------------------------------------------------

!+
SUBROUTINE MakeTable4()
! ------------------------------------------------------------------------------
! PURPOSE - Create Table4 with altitude, temperature, pressure, density in
!  US Customary units.
  
  INTEGER:: i                    ! geometric altitude in thousand feet (INTEGER)
  REAL:: altKm                                          ! altitude in kilometers
  REAL:: hkm                               ! geopotential altitude in kilometers
  REAL:: altKft                            ! geometric altitude in thousand feet
  REAL:: hkft                           ! geopotential altitude in thousand feet
  REAL:: gRatio,gus
  REAL:: soundspeed

  REAL:: sigma,delta,theta
!-------------------------------------------------------------------------------
  WRITE(HTML,FMT1) '<div id="table4">'

  WRITE(HTML,FMT1) &
    "<h2>Table 4 - Temperature, Pressure, & Density (US units)</h2>"
  WRITE(HTML,FMT1) "<h3>0 to 3000 Kft in steps of 10 Kft</h3>"
  WRITE(HTML,FMT1) "<table>"
  
  DO i=USMIN, USMAX, USDEL
    IF (MODULO(i,40*USDEL)==0) THEN 
      WRITE(HTML,FMT1) "<tr>"     
      WRITE(HTML,FMT1) TDC1,'Z',TD2, TDC1,'H',TD2
      WRITE(HTML,FMT1) TDC1,'T',TD2, TDC1,'T/T<sub>0</sub>',TD2
      WRITE(HTML,FMT1) TDC1,'Pressure',TD2, TDC1,'p/p<sub>0</sub>',TD2
      WRITE(HTML,FMT1) TDC1,'Mass<br />density',TD2
      WRITE(HTML,FMT1) TDC1,'&rho;/&rho;<sub>0</sub>',TD2
      WRITE(HTML,FMT1) TDC1,'Weight<br />density',TD2
      WRITE(HTML,FMT1) TDC1,'Sound<br />speed',TD2
      WRITE(HTML,FMT1) TDC1,'Gravity',TD2
      WRITE(HTML,FMT1) "</tr>"
    END IF   

    IF (i==0) THEN              ! only print at top of page
      WRITE(HTML,FMT1) "<tr>"
      WRITE(HTML,FMT1) TDC1,'kft',TD2, TDC1,'kft',TD2, TDC1,'R',TD2,BLANKCELL
      WRITE(HTML,FMT1) TDC1,'psf',TD2, BLANKCELL
      WRITE(HTML,FMT1) TDC1,'slug-ft<sup>-3</sup>',TD2, BLANKCELL
      WRITE(HTML,FMT1) TDC1,'lb-ft<sup>-3</sup>',TD2
      WRITE(HTML,FMT1) TDC1,'ft/s',TD2
      WRITE(HTML,FMT1) TDC1,'ft s<sup>-2</sup>',TD2
      WRITE(HTML,FMT1) "</tr>"
    END IF 
    
    IF (MODULO(i,2*USDEL)==0) THEN
      WRITE(HTML,FMT1) '<tr class="gray">'
    ELSE
      WRITE(HTML,FMT1) '<tr>'
    END IF

    altKft=REAL(i)
    altkm=FT2METERS*altKft
    hkm=altkm*REARTH/(altkm+REARTH)      ! convert geometric to geopotential
    hkft=hkm/FT2METERS
    CALL Atmosphere(altKm, sigma,delta,theta)
    gRatio=(REARTH/(altKm+REARTH))**2
    gus=gRatio*GZERO/FT2METERS       ! g in ft/sec**2   (about 32.2)
    soundspeed = (ASOUNDZERO/FT2METERS)*SQRT(theta)

    WRITE(HTML,'(A,I0,A)')     TD1, INT(altkft), TD2
    WRITE(HTML,'(A,F0.1,A)')   TD1, hkft, TD2
    WRITE(HTML,'(A,F0.3,A)')   TD1, theta*TZERO*KELVIN2RANKINE, TD2
    WRITE(HTML,'(A,F6.4,A)')   TD1,theta,TD2
    WRITE(HTML,'(A,ES12.4,A)') TD1,delta*PZERO/PSF2NSM,TD2    
    WRITE(HTML,'(A,ES12.4,A)') TD1,delta,TD2    
    WRITE(HTML,'(A,ES12.4,A)') TD1, sigma*RHOZERO/SCF2KCM, TD2
    WRITE(HTML,'(A,ES12.4,A)') TD1,sigma,TD2
    WRITE(HTML,'(A,ES12.4,A)') TD1, gus*sigma*RHOZERO/SCF2KCM, TD2    
    WRITE(HTML,'(A,F0.2,A)')   TD1,soundspeed,TD2  
    WRITE(HTML,'(A,F0.4,A)')   TD1,gus,TD2
    WRITE(HTML,FMT1) "</tr>"

  END DO
  
  WRITE(HTML,'(A)') '</table>', '</div>'
  
  RETURN
END Subroutine MakeTable4   ! ---------------------------------------------------


!+
SUBROUTINE MakeTable5()
! ------------------------------------------------------------------------------
! PURPOSE - Create Table 5 with gravity, presure scale height, number density, 
!  particle speed, collision frequency, mean free path
  
  INTEGER:: i
  REAL:: altKft                            ! geometric altitude in thousand feet
  REAL:: altKm                                ! geometric altitude in kilometers
  REAL:: hkft                           ! geopotential altitude in thousand feet
  REAL:: hkm                               ! geopotential altitude in kilometers
  REAL:: sigma,delta,theta
!  REAL:: temp,pressure,density, asound
  REAL:: dynamicViscosity,kinematicViscosity
  REAL:: kappa
  REAL:: collisionFrequency
  REAL:: g,gRatio
  REAL:: meanFreePath
  REAL:: mw
  REAL:: numDensity
  REAL:: partSpeed
  REAL:: pressureScaleHeight
  REAL:: dynamicViscosityRatio
  REAL:: kinematicViscosityRatio, numberDensityRatio, collisionFrequencyRatio
  REAL:: particleSpeedRatio
  REAL:: meanFreePathRatio,thermalConductivityRatio, soundSpeedRatio
  REAL:: molecularWeightRatio
  REAL:: xxx   ! temporary calculation value

!-------------------------------------------------------------------------------
  

  WRITE(HTML,FMT1) '<div id="table5">'
  WRITE(HTML,FMT1) "<h2>Table 5 - Transport Properties in US Customary units</h2>"
  WRITE(HTML,FMT1) "<h3>0 to 3000 Kft in steps of 10 Kft</h3>"
  WRITE(HTML,FMT1) "<table>"
 
  WRITE(HTML,FMT1) "<tr>"   ! only printed at top of table    
  WRITE(HTML,FMT1) TDC1,'alt',TD2, TDC1,'H',TD2
  WRITE(HTML,FMT1) TDC1,'Dynamic<br />viscosity',TD2
  WRITE(HTML,FMT1) TDC1,'Kinematic<br />viscosity',TD2
  WRITE(HTML,FMT1) TDC1,'Pressure<br />scale<br />height',TD2
  WRITE(HTML,FMT1) TDC1,'Number<br />density',TD2
  WRITE(HTML,FMT1) TDC1,'Particle<br />speed',TD2
  WRITE(HTML,FMT1) TDC1,'Collision<br />frequency',TD2
  WRITE(HTML,FMT1) TDC1,'Mean<br />free<br />path',TD2
  WRITE(HTML,FMT1) TDC1,'Thermal<br />conductivity',TD2
  WRITE(HTML,FMT1) TDC1,'Mol.<br />weight',TD2 
  WRITE(HTML,FMT1) "</tr>"


  DO i=USMIN, USMAX, USDEL
    IF (MODULO(i,40*USDEL)==0) THEN 
      WRITE(HTML,FMT1) "<tr>"     ! printed every 40 rows
      WRITE(HTML,FMT1) TDC1,'Z',TD2, TDC1,'H',TD2
      WRITE(HTML,FMT1) TDC1,'&mu;',TD2 
      WRITE(HTML,FMT1) TDC1,'&eta;',TD2
      WRITE(HTML,FMT1) TDC1,'H<sub>p</sub>',TD2
      WRITE(HTML,FMT1) TDC1,'n',TD2
      WRITE(HTML,FMT1) TDC1,'V',TD2
      WRITE(HTML,FMT1) TDC1,'&nu;',TD2
      WRITE(HTML,FMT1) TDC1,'L',TD2
      WRITE(HTML,FMT1) TDC1,'&kappa;',TD2 
      WRITE(HTML,FMT1) TDC1,'M',TD2
      WRITE(HTML,FMT1) "</tr>"
    END IF  

    IF (i == USMIN) THEN 
      WRITE(HTML,FMT1) "<tr>"   ! units row, only printed at top
      WRITE(HTML,FMT1) TDC1,'kft',TD2, TDC1,'kft',TD2
      WRITE(HTML,FMT1) TDC1,'slug ft<sup>-1</sup> s<sup>-1</sup>',TD2
      WRITE(HTML,FMT1) TDC1,'ft<sup>2</sup>/s',TD2
      WRITE(HTML,FMT1) TDC1,'ft',TD2, TDC1,'ft<sup>-3</sup>',TD2
      WRITE(HTML,FMT1) TDC1,'ft/s',TD2, TDC1,'s<sup>-1</sup>',TD2
      WRITE(HTML,FMT1) TDC1,'ft',TD2
      WRITE(HTML,FMT1) TDC1,'BTU/(s ft R)',TD2 
      WRITE(HTML,FMT1) BLANKCELL
      WRITE(HTML,FMT1) "</tr>"
    END IF 

    IF (MODULO(i, 2*USDEL) == 0) THEN
      WRITE(HTML,FMT1) '<tr class="gray">'
    ELSE
      WRITE(HTML,FMT1) '<tr>'
    END IF

    altKft=REAL(i)
    altKm=FT2METERS*altKft
    hkm=altKm*REARTH/(altKm+REARTH)      ! convert geometric to geopotential
    hkft=hkm/FT2METERS
    
    gRatio=(REARTH/(altKm+REARTH))**2
    g=gRatio*GZERO/FT2METERS
    CALL Atmosphere(altKm, sigma,delta,theta)
    mw=MolecularWeight(altKm) 
    CALL TransportRatios(altkm, gRatio,dynamicViscosityRatio, &
      kinematicViscosityRatio, numberDensityRatio, collisionFrequencyRatio, &
      particleSpeedRatio, meanFreePathRatio,thermalConductivityRatio, &
      soundSpeedRatio, molecularWeightRatio )

    dynamicViscosity=(MUZERO*FT2METERS/SLUG2KG)*dynamicViscosityRatio
    kinematicViscosity=(ETAZERO/FT2METERS**2)*kinematicViscosityRatio

    xxx=theta/(gRatio*molecularWeightRatio)
    pressureScaleHeight=(PRESSURE_SCALE_HEIGHT_ZERO/FT2METERS)*xxx
    numDensity=(NUM_DENSITY_ZERO*(FT2METERS**3))*numberDensityRatio
    partSpeed=(PART_SPEED_ZERO/FT2METERS)*particleSpeedRatio
    meanFreePath=(FREE_PATH_ZERO/FT2METERS)*meanFreePathRatio
    collisionFrequency=partSpeed/meanFreePath
    xxx=(KAPPA_ZERO*FT2METERS/(KELVIN2RANKINE*BTU2JOULE))
    kappa=xxx*thermalConductivityRatio 

    
    WRITE(HTML,'(A,I0,A)')     TD1, INT(altKft), TD2
    WRITE(HTML,'(A,F0.1,A)')   TD1, hkft, TD2
    WRITE(HTML,'(A,ES12.4,A)') TD1, dynamicViscosity, TD2
    WRITE(HTML,'(A,ES12.4,A)') TD1, kinematicViscosity, TD2
    WRITE(HTML,'(A,I0,A)')     TD1, NINT(pressureScaleHeight), TD2
    WRITE(HTML,'(A,ES12.4,A)') TD1, numDensity, TD2
    WRITE(HTML,'(A,F0.2,A)')   TD1, partSpeed, TD2
    WRITE(HTML,'(A,ES12.4,A)') TD1, collisionFrequency, TD2
    WRITE(HTML,'(A,ES12.4,A)') TD1, meanFreePath, TD2
    WRITE(HTML,'(A,ES12.4,A)') TD1, kappa, TD2
    WRITE(HTML,'(A,F0.3,A)')   TD1, mw, TD2
    
    WRITE(HTML,FMT1) "</tr>"

  END DO
  
  WRITE(HTML,'(A)') '</table>', '</div>'
  
  RETURN
END Subroutine MakeTable5   ! ---------------------------------------------------


!+
SUBROUTINE MakeTable6()
! ------------------------------------------------------------------------------
  INTEGER:: i
  REAL:: altKm                                          ! altitude in kilometers
  REAL:: hkm                               ! geopotential altitude in kilometers
  REAL:: altKft                                          ! altitude in kilo-feet
  REAL:: hkft                               ! geopotential altitude in kilo=feet
  
  REAL:: sigma,delta,theta
!  REAL:: temp,pressure,density, asound
!  REAL:: dynamicViscosity,kinematicViscosity
!  REAL:: kappa
  
  REAL:: gRatio,dynamicViscosityRatio
  REAL:: kinematicViscosityRatio, numberDensityRatio, collisionFrequencyRatio
  REAL:: particleSpeedRatio
  REAL:: meanFreePathRatio,thermalConductivityRatio, soundSpeedRatio
  REAL:: molecularWeightRatio
!-------------------------------------------------------------------------------
  
  WRITE(HTML,FMT1) '<div id="table6">'
  WRITE(HTML,FMT1) "<h2>Table 6 - Transport Properties (Non-dimensional)</h2>"
  WRITE(HTML,FMT1) "<h3>0 to 3000 Kft in steps of 10 Kft</h3>"
  WRITE(HTML,FMT1) "<p>The units of alt and H are thousands of feet.</p>"
  WRITE(HTML,FMT1) "<table>"

  DO i=USMIN,USMAX,USDEL
    IF (MODULO(i,40*USDEL)==0) THEN    
      WRITE(HTML,FMT1) "<tr>"     
      WRITE(HTML,FMT1) TDC1,'Z',TD2, TDC1,'H',TD2
      WRITE(HTML,FMT1) TDC1,'&mu;/&mu;<sub>0</sub>',TD2
      WRITE(HTML,FMT1) TDC1,'&eta;/&eta;<sub>0</sub>',TD2
      WRITE(HTML,FMT1) TDC1,'n/n<sub>0</sub>',TD2
      WRITE(HTML,FMT1) TDC1,'V/V<sub>0</sub>',TD2
      WRITE(HTML,FMT1) TDC1,'&nu;/&nu;<sub>0</sub>',TD2
      WRITE(HTML,FMT1) TDC1,'L/L<sub>0</sub>',TD2
      WRITE(HTML,FMT1) TDC1,'&kappa;/&kappa;<sub>0</sub>',TD2
      WRITE(HTML,FMT1) TDC1,'c/c<sub>0</sub>',TD2
      WRITE(HTML,FMT1) TDC1,'M/M<sub>0</sub>',TD2
      WRITE(HTML,FMT1) "</tr>"
    END IF 
    
    IF (MODULO(i,2*USDEL)==0) THEN  ! Make every other row with light gray background.
      WRITE(HTML,FMT1) '<tr class="gray">'
    ELSE
      WRITE(HTML,FMT1) '<tr>'
    END IF

    altKft=REAL(i)
    altKm=altKft*FT2METERS
    hkm=altkm*REARTH/(altkm+REARTH)      ! convert geometric to geopotential
    hkft=hkm/FT2METERS
    
    CALL Atmosphere(altKm, sigma,delta,theta)
    CALL TransportRatios(altkm, gRatio,dynamicViscosityRatio, &
      kinematicViscosityRatio, numberDensityRatio, collisionFrequencyRatio, &
      particleSpeedRatio, meanFreePathRatio,thermalConductivityRatio, &
      soundSpeedRatio, molecularWeightRatio )

    WRITE(HTML,'(A,I0,A)') TD1, INT(altKft), TD2
    WRITE(HTML,'(A,F0.1,A)') TD1, hkft, TD2   
    WRITE(HTML,'(A,F6.4,A)') TD1, dynamicViscosityRatio, TD2
    WRITE(HTML,'(A,ES12.4,A)') TD1, kinematicViscosityRatio, TD2
    WRITE(HTML,'(A,ES12.4,A)') TD1, numberDensityRatio, TD2
    WRITE(HTML,'(A,F6.4,A)') TD1, particleSpeedRatio, TD2
    WRITE(HTML,'(A,ES12.4,A)') TD1, collisionFrequencyRatio, TD2
    WRITE(HTML,'(A,ES12.4,A)') TD1, meanFreePathRatio, TD2
    WRITE(HTML,'(A,F6.4,A)') TD1, thermalConductivityRatio, TD2
    WRITE(HTML,'(A,F6.4,A)') TD1, soundSpeedRatio, TD2
    WRITE(HTML,'(A,F6.4,A)') TD1, molecularWeightRatio, TD2 
    WRITE(HTML,FMT1) "</tr>"

  END DO
  
  WRITE(HTML,'(A)') '</table>','</div>'

  RETURN
END Subroutine MakeTable6   ! ---------------------------------------------------


!+
SUBROUTINE MakeFooter(updater)
! ------------------------------------------------------------------------------
! PURPOSE - Write the footer and close the HTML file properly
  CHARACTER(LEN=*),INTENT(IN):: updater
  CHARACTER(LEN=11):: dateStr
!-------------------------------------------------------------------------------
  dateStr=GetDateStr()
  WRITE(HTML,FMT1) '<footer>'
  WRITE(HTML,FMT1) 'Last updated: ', dateStr, ' by '//updater
  WRITE(HTML,FMT1) '<a href="mailto:pdaerowebmaster@gmail.com">', &
                   'pdaerowebmaster AT gmail DOT com</a>'
  WRITE(HTML,FMT1) '</footer>'
  WRITE(HTML,FMT1) '<div class="banner">', PDAS, '</div>'
  WRITE(HTML,FMT1) '<a href="index.html">PDAS home</a> &gt; ', &
                   '<a href="atmos.html">Standard Atmosphere</a> &gt; ', &
                    'BigTables'  
  WRITE(HTML,FMT1) '</body>'
  WRITE(HTML,FMT1) '</html>'
  
  CLOSE(UNIT=HTML)

  RETURN 
END Subroutine MakeFooter   ! --------------------------------------------------  

!+
SUBROUTINE TextTable()
! ------------------------------------------------------------------------------

  CHARACTER(LEN=*),PARAMETER::HEAD1=' alt    sigma     delta   theta  temp'// &
    '   press    dens     a    visc  k.visc'
  CHARACTER(LEN=*),PARAMETER::HEAD2='  Km                            kelvins'// &
    ' N/sq.m   kg/cu.m m/s   kg/m-s sq.m/s'
  CHARACTER(LEN=*),PARAMETER:: FMT = &
    '(I4,2ES11.4,F7.4,F6.1,2ES10.3,F6.1,F6.2,2ES9.2)'
  INTEGER,PARAMETER::ITXT=1
  INTEGER:: errCode
  INTEGER:: i
  REAL:: altKm                                ! altitude in kilometers
  REAL:: sigma,delta,theta
  REAL:: temp,pressure,density, asound
  REAL:: visc,kinematicViscosity
!-------------------------------------------------------------------------------
  OPEN(UNIT=ITXT, FILE='bt3text.txt', STATUS='REPLACE', &
    IOSTAT=errCode, ACTION='WRITE', POSITION='REWIND')
  IF (errCode /= 0) THEN
    WRITE(*,*) "Unable to open output file bt3text.txt"
    STOP
  END IF

  DO i=0,1000,5
    IF (MODULO(i,5*40)==0) THEN   ! titles every 40 lines
      WRITE(ITXT,*) HEAD1
      WRITE(ITXT,*) HEAD2
    END IF 
    
    altKm=REAL(i)
    CALL Atmosphere(altKm, sigma,delta,theta)
    temp=TZERO*theta
    pressure=PZERO*delta
    density=RHOZERO*sigma
    asound=ASOUNDZERO*Sqrt(theta)
    visc=Viscosity(theta)
    kinematicViscosity=visc/density
    WRITE(ITXT,FMT) i, sigma, delta, theta, temp, pressure,  &
      density, asound, 1.0E6*visc, kinematicViscosity
  END DO
  Close(Itxt)
  RETURN
END Subroutine TextTable   ! ---------------------------------------------------

END Module BigTablesProcedures   ! =============================================

!+
PROGRAM BigTables
! ------------------------------------------------------------------------------
! PURPOSE - Make tables of atmospheric properties from 0 to 1000 km.
!
! REFERENCE - "The U.S. Standard Atmosphere, 1976" published by the 
!    Government Printing Office. Also filed as NASA TM X-74335.

! AUTHOR - Ralph L. Carmichael, Public Domain Aeronautical Software
USE ISO_FORTRAN_ENV

USE BigTablesProcedures
IMPLICIT NONE
! NOTES-

  CHARACTER(LEN=*),PARAMETER::GREETING =  &
    'bt3 - A Fortran program to make atmosphere tables to 1000 km'
  CHARACTER(LEN=*),PARAMETER::AUTHOR=  &
    ' Ralph L. Carmichael, Public Domain Aeronautical Software'
  CHARACTER(LEN=*),PARAMETER::MODIFIER=' '   ! put your name here if you mod it
  CHARACTER(LEN=*),PARAMETER::VERSION=' 2.0 (2022 March 18)'
  CHARACTER(LEN=*),PARAMETER::FAREWELL= &
    'File bigtablesf.html added to your directory.'
  REAL:: time1,time2   
!-------------------------------------------------------------------------------
  CALL CPU_TIME(time1)
  WRITE(*,*) GREETING
  WRITE(*,*) AUTHOR
  IF (MODIFIER /= ' ') WRITE(*,*) 'Modified by '//MODIFIER
  WRITE(*,*) VERSION
  WRITE(*,'(3A/A)') 'This file was compiled by ', compiler_version(), &
                    ' using the options ',        compiler_options()

! Check some of the program constants...
!  WRITE(*,'(A,ES12.4)') 'num density', NUM_DENSITY_ZERO, &
!    'particle speed sq', PART_SPEED_ZERO_SQ, &
!    'particle speed', SQRT(PART_SPEED_ZERO_SQ), &
!    'mean free path', FREE_PATH_ZERO, &
!    'scale height', PRESSURE_SCALE_HEIGHT_ZERO

  CALL MakeMainPage()
  CALL MakeTable1()
  CALL MakeTable2()
  CALL MakeTable3()
  CALL MakeTable4()
  CALL MakeTable5()
  CALL MakeTable6()
  IF (MODIFIER == ' ') THEN 
    CALL MakeFooter(AUTHOR)
  ELSE    
    CALL MakeFooter(MODIFIER)
  END IF   

  WRITE(*,*) FAREWELL
  WRITE(*,*) 'bigtables has terminated successfully.'
  CALL CPU_TIME(time2)
  WRITE(*,'(A,I0,A)') 'CPU time=', INT(1000*(time2-time1)),' ms'
  STOP


END Program BigTables   ! =======================================================

