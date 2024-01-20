#!/usr/bin/env python3
#!/usr/bin/python3

"""
bt9.py - Make tables  of atmospheric properties   (Python 3)
   Write a file in HTML 5 format that displays tables of atmospheric 
   properties from 0 to 1000 km geometric altitude.
   There are six different tables...





Adapted by
    Richard J. Kwan, Lightsaber Computing
from original programs by
    Ralph L. Carmichael, Public Domain Aeronautical Software

Revision History
Date         Vers Person Statement of Changes
2004 Oct 04  1.0  RJK    Initial program
2017 Jun 04  1.1  RLC    All indents are with spaces; all prints in ( )
                         New version for Python 3 that does integer div with //
2018 Aug 13  1.2  RLC    Each table in a <div> with id. Buttons to control view.
2019 Jun 25  1.3  RLC    LOGDELTA and LOGSIGMA with 5 decimals.
2019 Jun 29  1.31 RLC    Made text closer to F90 version
2019 Sep 28  1.4  RLC    Last cleanup for V17 PDAS
2022 Mar 18  1.5  RLC    REARTH = 6356.766
"""

import sys, math, time

version = "1.5 (2022 Mar 18)"
greeting = "bigtables - A Python program to display atmosphere tables"
author = "Ralph L. Carmichael, Public Domain Aeronautical Software"
modifier = ""   # put your name here if you change anything
farewell = "File bigtables.html added to your directory."
finalmess = "Normal termination of bigtables."
PDAS = "Public Domain Aeronautical Software (PDAS) &nbsp;"


PI = 3.14159265
#   P H Y S I C A L   C O N S T A N T S

FT2METERS = 0.3048      # mult. ft. to get meters (exact)
KELVIN2RANKINE = 1.8    # mult deg K to get deg R
PSF2NSM = 47.880258     # mult lb/sq.ft to get sq.m
SCF2KCM = 515.379       # mult slugs/cu.ft to get kg/cu.m
SLUG2KG = 14.594        # mult slugs to get kilograms
BTU2JOULE = 1055.0      # mult BTU to get joules
BETAVISC = 1.458E-6     # viscosity constant
SUTH    = 110.4         # Sutherland's constant, kelvins
AVOGADRO =  6.022169E26   # 1/kmol, Avogadro constant
BOLTZMANN = 1.380622E-23  # Nm/K, Boltzmann constant
REARTH  = 6356.766     # polar radius of the earth, kilometers
RSTAR = 8314.32          # perfect gas constant, N-m/(kmol-K)

# NOTE - Observe the factor of 1000 in computing GMR. Without it, GMR would have
#   units of kelvins per meter. But the temperature gradient in LowerAtmosphere
#   is given in kelvins per kilometer.   GMR=34.163195 kelvins/km

#   S E A   L E V E L   C O N D I T I O N S 
GZERO = 9.80665         #  sea level accel. of gravity, m/s^2
MOLWT_ZERO = 28.9644    # molecular weight of air at sea level
GMR = 1000*GZERO*MOLWT_ZERO/RSTAR # hydrostatic constant, kelvins/km
TZERO   = 288.15        # sea-level temperature, kelvins
PZERO   = 101325.0      # sea-level pressure, N/sq.m
RHOZERO = 1.225         # sea-level density, kg/cu.m
ASOUNDZERO = 340.294    # speed of sound at S.L.  m/sec
MUZERO = 1.7894E-5      # viscosity at sea-level,
ETAZERO = 1.4607E-5     # kinetic viscosity at sea-level,
KAPPA_ZERO = 0.025326   # thermal coeff. at sea-level
                        # units of KAPPA_ZERO are watts per meter per kelvin
CROSS = 3.65E-10 # collision cross-section, diatomic air m^2
NUM_DENSITY_ZERO = RHOZERO*AVOGADRO/MOLWT_ZERO #  1/m^3
PART_SPEED_ZERO = math.sqrt((8.0/PI)*RSTAR*TZERO/MOLWT_ZERO) # m/s
FREE_PATH_ZERO = RSTAR*TZERO/(math.sqrt(2.0)*PI*AVOGADRO*CROSS*CROSS*PZERO) # m
PRESSURE_SCALE_HEIGHT_ZERO = RSTAR*TZERO/(MOLWT_ZERO*GZERO)# m

# Initial and final values and increments of SI and US tables
SIMIN=0; SIMAX=1000; SIDEL=5;   USMIN=0; USMAX=3000; USDEL=10

# ------------------------------------------------------------------------------
# PURPOSE - Evaluate a cubic polynomial defined by the function and the
#   1st derivative at two points.
#       Usually  a < u < b, but not necessary.
#  REAL,INTENT(IN):: u   # point where function is to be evaluated
#  REAL,INTENT(IN):: a,fa,fpa   # a, f(a), f'(a)  at first point
#  REAL,INTENT(IN):: b,fb,fpb   # b, f(b), f'(b)  at second point
def EvaluateCubic(a,fa,fpa, b,fb,fpb, u):
    d = (fb-fa)/(b-a)   # b==a will be a fatal error
    t = (u-a)/(b-a)
    p = 1.0-t

    fu = p*fa + t*fb - p*t*(b-a)*(p*(d-fpa)-t*(d-fpb))
    return fu                                 # -- End of Function EvaluateCubic

#   ----------------------------------------------------------------------------
# PURPOSE - Compute kinetic temperature above 86 km.
#  REAL,INTENT(IN)::  z     # geometric altitude, km.

def KineticTemperature(z):
    C1 = -76.3232  # uppercase A in document
    C2 = 19.9429   # lowercase a in document
    C3 = 12.0
    C4 = 0.01875   # lambda in document
    TC = 263.1905
#    Z7 =  86.0
    T7 = 186.8673
    Z8 =  91.0
#    T8=T7
    Z9 = 110.0
    T9=240.0
    Z10= 120.0
    T10=360.0
#  REAL,PARAMETER:: Z11= 500.0, T11=999.2356   # not used
#    Z12=1000.0
    T12=1000.0   # T12 is also T-sub-infinity

# NOTE - REARTH is a global variable


    if z <= Z8:
        t = T7    # Eq. (25), p.11
    elif z < Z9:
        xx = (z-Z8)/C2                      
        yy = math.sqrt(1.0-xx*xx)
        t = TC+C1*yy   # Eq. (27)
    elif z <= Z10:
        t = T9+C3*(z-Z9)   # Eq. (29)
    else:
        xx = (REARTH+Z10)/(REARTH+z)
        yy = (T12-T10)*math.exp(-C4*(z-Z10)*xx)
        t = T12-yy   # Eq. (31)

    return t    # --------------------------- End of Function KineticTemperature


#   ----------------------------------------------------------------------------
# PURPOSE - Compute k, the coefficient of thermal conductivity at temperature t.
#  k has units of watts per meter per kelvin
# This is a coding of section 1.3.13 on p. 19 of the reference document.

def ThermalConductivity(t):      # t is temperature, K
    C1=2.64638E-3; C2=245.4
    k=C1*math.sqrt(t*t*t)/(t+C2*10.0**(-12.0/t))    # Eq. (53)

    return(k)   # -------------------------- End of Function ThermalConductivity




# ------------------------------------------------------------------------------
# PURPOSE -    Compute temperature, density, and pressure in the lower 86 km of
#  the standard atmosphere.
#    Input:
#    alt geometric altitude, km.
#    Return: (sigma, delta, theta)
#    sigma   density/sea-level standard density
#    delta   pressure/sea-level standard pressure
#    theta   temperature/sea-level std. temperature

def LowerAtmosphere(alt):
    htab = [ 0.0,  11.0, 20.0, 32.0, 47.0,
             51.0, 71.0, 84.852 ]
    ttab = [ 288.15, 216.65, 216.65, 228.65, 270.65,
             270.65, 214.65, 186.946 ]
    ptab = [ 1.0, 2.2336110E-1, 5.4032950E-2, 8.5666784E-3, 1.0945601E-3,
             6.6063531E-4, 3.9046834E-5, 3.68501E-6 ]
    gtab = [ -6.5, 0.0, 1.0, 2.8, 0, -2.8, -2.0, 0.0 ]

# NOTE - REARTH and GMR are global variables
    h = alt*REARTH/(alt+REARTH) # geometric to geopotential altitude

    i = 0; j = len(htab)
    while (j > i+1):
        k = (i+j)//2      # this is floor division in Python 3
        if h < htab[k]:
            j = k
        else:
            i = k
    tgrad = gtab[i]     # temp. gradient of local layer
    tbase = ttab[i]     # base  temp. of local layer
    deltah=h-htab[i]        # height above local base
    tlocal=tbase+tgrad*deltah   # local temperature
    theta = tlocal/ttab[0]  # temperature ratio

    if 0.0 == tgrad:
        delta=ptab[i]*math.exp(-GMR*deltah/tbase)
    else:
        delta=ptab[i]*math.pow(tbase/tlocal, GMR/tgrad)
    sigma = delta/theta
    return ( sigma, delta, theta )   # --------- end of function LowerAtmosphere


#   ----------------------------------------------------------------------------
# PURPOSE - Compute the properties of the 1976 standard atmosphere from
#   86 km. to 1000 km.
#    Input:
#    alt geometric altitude, km.
#    Return: (sigma, delta, theta)
#    sigma   density/sea-level standard density
#    delta   pressure/sea-level standard pressure
#    theta   temperature/sea-level std. temperature

def UpperAtmosphere(alt):
    # altitude table (km)  len(Z)=25
    Z = [86.0,  93.0, 100.0, 107.0, 114.0,
    121.0, 128.0, 135.0, 142.0, 150.0,
    160.0, 170.0, 180.0, 190.0, 200.0,
    220.0, 260.0, 300.0, 400.0, 500.0,
    600.0, 700.0, 800.0, 900.0,1000.0]

# pressure table  (DELTA = P/P0)
    DELTA_TABLE = [3.6850E-6, 1.0660E-6, 3.1593E-7, 1.0611E-7, 4.3892E-8,
    2.3095E-8, 1.3997E-8, 9.2345E-9, 6.4440E-9, 4.4828E-9,
    2.9997E-9, 2.0933E-9, 1.5072E-9, 1.1118E-9, 8.3628E-10,
    4.9494E-10, 1.9634E-10, 8.6557E-11, 1.4328E-11, 2.9840E-12,
    8.1056E-13, 3.1491E-13, 1.6813E-13, 1.0731E-13, 7.4155E-14]

# density table  (SIGMA = RHO/RHO0)
    SIGMA_TABLE = [5.680E-6, 1.632E-6, 4.575E-7, 1.341E-7, 4.061E-8,
    1.614E-8, 7.932E-9, 4.461E-9, 2.741E-9, 1.694E-9,
    1.007E-9, 6.380E-10, 4.240E-10, 2.923E-10, 2.074E-10,
    1.116E-10, 3.871E-11, 1.564E-11, 2.288E-12, 4.257E-13,
    9.279E-14, 2.506E-14, 9.272E-15, 4.701E-15, 2.907E-15 ]

#    LOGDELTA = [
#    -1.2511E+01, -1.3752E+01, -1.4968E+01, -1.6059E+01, -1.6942E+01,
#    -1.7584E+01, -1.8084E+01, -1.8500E+01, -1.8860E+01, -1.9223E+01,
#    -1.9625E+01, -1.9985E+01, -2.0313E+01, -2.0617E+01, -2.0902E+01,
#    -2.1427E+01, -2.2351E+01, -2.3170E+01, -2.4969E+01, -2.6538E+01,
#    -2.7841E+01, -2.8786E+01, -2.9414E+01, -2.9863E+01, -3.0233E+01]

#    LOGSIGMA = [
#    -1.2079E+01, -1.3326E+01, -1.4597E+01, -1.5825E+01, -1.7019E+01,
#    -1.7942E+01, -1.8652E+01, -1.9228E+01, -1.9715E+01, -2.0196E+01,
#    -2.0716E+01, -2.1173E+01, -2.1581E+01, -2.1953E+01, -2.2296E+01,
#    -2.2916E+01, -2.3975E+01, -2.4881E+01, -2.6803E+01, -2.8485E+01,
#    -3.0008E+01, -3.1318E+01, -3.2312E+01, -3.2991E+01, -3.3472E+01]

#    LOGDELTA = [
#    -1.25112E+01,-1.37516E+01,-1.49677E+01,-1.60588E+01,-1.69415E+01,
#    -1.75837E+01,-1.80844E+01,-1.85003E+01,-1.88601E+01,-1.92230E+01,
#    -1.96248E+01,-1.99845E+01,-2.03130E+01,-2.06173E+01,-2.09021E+01,
#    -2.14266E+01,-2.23512E+01,-2.31702E+01,-2.49688E+01,-2.65378E+01,
#    -2.78411E+01,-2.87865E+01,-2.94140E+01,-2.98631E+01,-3.02326E+01]
#    LOGSIGMA = [
#    -1.20786E+01,-1.33257E+01,-1.45975E+01,-1.58247E+01,-1.70193E+01,
#    -1.79420E+01,-1.86524E+01,-1.92279E+01,-1.97149E+01,-2.01962E+01,
#    -2.07163E+01,-2.11727E+01,-2.15813E+01,-2.19532E+01,-2.22964E+01,
#    -2.29161E+01,-2.39749E+01,-2.48812E+01,-2.68033E+01,-2.84850E+01,
#    -3.00084E+01,-3.13175E+01,-3.23118E+01,-3.29910E+01,-3.34717E+01]
#    LOGDELTA = [
#    -1.251124E+01, -1.375160E+01, -1.496774E+01, -1.605879E+01, -1.694153E+01,
#    -1.758365E+01, -1.808442E+01, -1.850032E+01, -1.886012E+01, -1.922302E+01,
#    -1.962475E+01, -1.998452E+01, -2.031301E+01, -2.061728E+01, -2.090206E+01,
#    -2.142658E+01, -2.235117E+01, -2.317022E+01, -2.496881E+01, -2.653776E+01,
#    -2.784105E+01, -2.878649E+01, -2.941404E+01, -2.986305E+01, -3.023262E+01]
#    LOGSIGMA = [
#    -1.207856E+01, -1.332570E+01, -1.459749E+01, -1.582468E+01, -1.701925E+01,
#    -1.794197E+01, -1.865236E+01, -1.922789E+01, -1.971494E+01, -2.019617E+01,
#    -2.071629E+01, -2.117268E+01, -2.158129E+01, -2.195324E+01, -2.229637E+01,
#    -2.291610E+01, -2.397492E+01, -2.488119E+01, -2.680334E+01, -2.848504E+01,
#    -3.000844E+01, -3.131750E+01, -3.231178E+01, -3.299100E+01, -3.347165E+01]
 
    LOGDELTA = DELTA_TABLE[:]
    LOGSIGMA = SIGMA_TABLE[:]



    DLOGDELTA = [-0.174061, -0.177924, -0.167029, -0.142755, -0.107859,
    -0.079322, -0.064664, -0.054879, -0.048260, -0.042767,
    -0.037854, -0.034270, -0.031543, -0.029384, -0.027632,
    -0.024980, -0.021559, -0.019557, -0.016735, -0.014530,
    -0.011314, -0.007677, -0.005169, -0.003944, -0.003612]

    DLOGSIGMA = [-0.172421, -0.182258, -0.178090, -0.176372, -0.154322,
    -0.113750, -0.090582, -0.075033, -0.064679, -0.056067,
    -0.048461, -0.043042, -0.038869, -0.035648, -0.033063,
    -0.029164, -0.024220, -0.021336, -0.017686, -0.016035,
    -0.014327, -0.011631, -0.008248, -0.005580, -0.004227]

    for k in range(0,25):
        LOGDELTA[k]=math.log(LOGDELTA[k])
        LOGSIGMA[k]=math.log(LOGSIGMA[k])

    if alt >= Z[24] :          # trap altitudes greater than 1000 km.
        delta=DELTA_TABLE[24]        # return value for Z=1000.0
        sigma=SIGMA_TABLE[24]        #   ditto
        theta=1000.0/TZERO
    else:
        i=0
        j=len(Z)-1                               # setting up for binary search
        while j > i+1 :
            k=(i+j)//2                           # integer (floor) division
            if alt < Z[k] :
                j=k
            else:
                i=k


#        if alt>900.0:
#            print(alt,i,j,k,Z[i])
        
        delta=math.exp(EvaluateCubic(Z[i], LOGDELTA[i], DLOGDELTA[i],
        Z[i+1], LOGDELTA[i+1], DLOGDELTA[i+1], alt))

        sigma=math.exp(EvaluateCubic(Z[i], LOGSIGMA[i],   DLOGSIGMA[i],
                        Z[i+1], LOGSIGMA[i+1], DLOGSIGMA[i+1], alt))

        theta=KineticTemperature(alt)/TZERO

    return (sigma, delta, theta)   # ----------- End of Function UpperAtmosphere

#   ----------------------------------------------------------------------------
# PURPOSE - Compute the properties of the 1976 standard atmosphere from
#   0 km. to 1000 km.
#    Input:
#    alt geometric altitude, km.
#    Return: (sigma, delta, theta)
#    sigma   density/sea-level standard density
#    delta   pressure/sea-level standard pressure
#    theta   temperature/sea-level std. temperature

def Atmosphere(alt):
    if alt > 86.0 :
        (sigma,delta,theta) = UpperAtmosphere(alt)
    else:
        (sigma,delta,theta) = LowerAtmosphere(alt)

    return(sigma,delta,theta)   # ------------------- End of Function Atmosphere


def Viscosity(theta):
    t=theta*TZERO
    return BETAVISC*math.sqrt(t*t*t)/(t+SUTH)   # ---- End of Function Viscosity

# ------------------------------------------------------------------------------
# PURPOSE - Compute molecular weight of air at altitude by interpolation
#    in tables.
def MolecularWeight(altKm):

# tables of altitude (km), mol.wt., and d(mol.wt.)/d(altKm)
    Z = [86.0,  93.0, 100.0, 107.0, 114.0,
    121.0, 128.0, 135.0, 142.0, 150.0,
    160.0, 170.0, 180.0, 190.0, 200.0,
    220.0, 260.0, 300.0, 400.0, 500.0,
    600.0, 700.0, 800.0, 900.0,1000.0]

    M = [28.95, 28.82, 28.40, 27.64, 26.79,
    26.12, 25.58, 25.09, 24.62, 24.10,
    23.49, 22.90, 22.34, 21.81, 21.30,
    20.37, 18.85, 17.73, 15.98, 14.33,
    11.51,  8.00,  5.54,  4.40,  3.94]

    MP = [ -0.001340, -0.036993, -0.086401, -0.123115, -0.111136,
    -0.083767, -0.072368, -0.068190, -0.066300, -0.063131,
    -0.059786, -0.057724, -0.054318, -0.052004, -0.049665,
    -0.043499, -0.032674, -0.023804, -0.014188, -0.021444,
    -0.034136, -0.031911, -0.017321, -0.006804, -0.003463]


    if altKm <= 86.0:
        mw = MOLWT_ZERO
    elif altKm >= 1000.0:
        mw = M[24]
    else:
        i=0
        j=len(Z)-1                                # setting up for binary search
#        print(j)
        while (j>i+1):
            k=(i+j)//2                                        # integer division
            if altKm < Z[k] :
                j=k
            else:
                i=k
#            print(i)
            mw = EvaluateCubic(Z[i], M[i], MP[i], Z[i+1], M[i+1], MP[i+1], altKm)

    return mw   # ------------------------------ End of Function MolecularWeight

def TransportRatios(alt):      # alt is geometric altitude in kilometers
    (densityRatio,pressureRatio,temperatureRatio) = Atmosphere(alt)
    molwt=MolecularWeight(alt)
    molecularWeightRatio=molwt/MOLWT_ZERO

    gRatio=(REARTH/(alt+REARTH))**2
    dynamicViscosityRatio=Viscosity(temperatureRatio)/MUZERO
    kinematicViscosityRatio=dynamicViscosityRatio/densityRatio
    numberDensityRatio=densityRatio*MOLWT_ZERO/molwt
    particleSpeedRatio=math.sqrt(temperatureRatio/molecularWeightRatio)
    meanFreePathRatio=molecularWeightRatio/densityRatio
    collisionFrequencyRatio=particleSpeedRatio/meanFreePathRatio
    kappa=ThermalConductivity(temperatureRatio*TZERO)
    thermalConductivityRatio=kappa/KAPPA_ZERO
    soundSpeedRatio=math.sqrt(temperatureRatio)


    return(gRatio,dynamicViscosityRatio,kinematicViscosityRatio,
    numberDensityRatio, collisionFrequencyRatio,particleSpeedRatio,
    meanFreePathRatio,thermalConductivityRatio, soundSpeedRatio,
    molecularWeightRatio)
                            # ------------------ End of Function TransportRatios


# ------------------------------------------------------------------------------
# PURPOSE - Open the output file bigtables.html and write the HTML header, the
#  CSS style sheets, and the javascript functions.
#  The <head> element contains the <meta charset=utf-8 />, <title>, <style>
#  and <script> elements.
# NOTE - When padding property has four values: top-right-bottom-left
#        When padding property has two values: top & bottom-right & left
#        When padding property has one value: all
#-------------------------------------------------------------------------------

def MakeMainPage(ihtml):
    ihtml.write('<!DOCTYPE html>\n<html lang="en">\n<head>\n')
    ihtml.write('<meta charset="utf-8" />\n')
    ihtml.write('<title>Table of the  U.S. Standard Atmosphere 1976</title>\n')

    ihtml.write('<style>\n')
    ihtml.write('body {font-size: 99%;}\n')
    ihtml.write('table {border: none;}\n')
    ihtml.write('td {border: none; padding: 0 0 0.2em 0.2em; text-align: right;}\n')
    ihtml.write('td.center {text-align:center}\n')
    ihtml.write('td.right {text-align:right}\n')
    ihtml.write('tr.gray {background-color: #f0f0f0;}\n')

    ihtml.write('header {width: 100%; ')
    ihtml.write('border-bottom: 3px solid blue; padding: 3px; }\n')
    ihtml.write('footer {width: 100%; font-size: 80%; \n')
    ihtml.write('   padding: 3px; border-top: 3px solid red; }\n')
    ihtml.write('h1 { font-family: Helvetica,Arial,sans-serif;\n')
    ihtml.write('   text-align: center; font-weight: bold; font-size: 2.0em;\n')
    ihtml.write('   color: blue; margin-top: 0.5em; }\n')

    ihtml.write('#MainButton   {margin: 10px; color: black;\n')
    ihtml.write('background-color: white; font: 1em Arial,sans-serif;}\n')
    ihtml.write('#Table1Button {margin: 10px; color: black;\n')
    ihtml.write('background-color: red; font: 1em Arial,sans-serif;}\n')
    ihtml.write('#Table2Button {margin: 10px; color: black;\n')
    ihtml.write('background-color: orange; font: 1em Arial,sans-serif;}\n')
    ihtml.write('#Table3Button {margin: 10px; color: black;\n')
    ihtml.write('background-color: yellow; font: 1em Arial,sans-serif;}\n')
    ihtml.write('#Table4Button {margin: 10px; color: white;\n')
    ihtml.write('background-color: green; font: 1em Arial,sans-serif;}\n')
    ihtml.write('#Table5Button {margin: 10px; color: white;\n')
    ihtml.write('background-color: blue; font: 1em Arial,sans-serif;}\n')
    ihtml.write('#Table6Button {margin: 10px; color: white;\n')
    ihtml.write('background-color: violet; font: 1em Arial,sans-serif;}\n')

    ihtml.write('#table1 {display:none;}\n')
    ihtml.write('#table2 {display:none;}\n')
    ihtml.write('#table3 {display:none;}\n')
    ihtml.write('#table4 {display:none;}\n')
    ihtml.write('#table5 {display:none;}\n')
    ihtml.write('#table6 {display:none;}\n')

    ihtml.write('div.banner {width: 100%; background: blue; color: white;\n')
    ihtml.write('  font-size: 120%; font-family: Verdana,Arial,Helvetica,\n')
    ihtml.write('  sans-serif; text-align: right;}\n')

    ihtml.write('</style>\n')   # ------------------------------- End of <style>

    ihtml.write('<script>\n')

    ihtml.write('function HideAll() {\n' )
    ihtml.write('  document.getElementById("main").style.display="none";\n')
    ihtml.write('  document.getElementById("table1").style.display="none";\n')
    ihtml.write('  document.getElementById("table2").style.display="none";\n')
    ihtml.write('  document.getElementById("table3").style.display="none";\n')
    ihtml.write('  document.getElementById("table4").style.display="none";\n')
    ihtml.write('  document.getElementById("table5").style.display="none";\n')
    ihtml.write('  document.getElementById("table6").style.display="none";\n')
    ihtml.write('}   // --- End of function HideAll\n')

    ihtml.write('function ShowMain() {\n')
    ihtml.write('  HideAll();\n')
    ihtml.write('  document.getElementById("main").style.display="block";\n')
    ihtml.write('}   // --- End of function ShowMain\n')

    ihtml.write('function ShowTable1() {\n')
    ihtml.write('  HideAll();\n')
    ihtml.write('  document.getElementById("table1").style.display="block";\n')
    ihtml.write('}   // --- End of function ShowTable1\n')

    ihtml.write('function ShowTable2() {\n')
    ihtml.write('  HideAll();\n')
    ihtml.write('  document.getElementById("table2").style.display="block";\n')
    ihtml.write('}   // --- End of function ShowTable2\n')

    ihtml.write('function ShowTable3() {\n')
    ihtml.write('  HideAll();\n')
    ihtml.write('  document.getElementById("table3").style.display="block";\n')
    ihtml.write('}   // --- End of function ShowTable3\n')

    ihtml.write('function ShowTable4() {\n')
    ihtml.write('  HideAll();\n')
    ihtml.write('  document.getElementById("table4").style.display="block";\n')
    ihtml.write('}   // --- End of function ShowTable4\n')

    ihtml.write('function ShowTable5() {\n')
    ihtml.write('  HideAll();\n')
    ihtml.write('  document.getElementById("table5").style.display="block";\n')
    ihtml.write('}   // --- End of function ShowTable5\n')

    ihtml.write('function ShowTable6() {\n')
    ihtml.write('  HideAll();\n')
    ihtml.write('  document.getElementById("table6").style.display="block";\n')
    ihtml.write('}   // --- End of function ShowTable6\n')

    ihtml.write('</script>\n')   # ----------------------------- End of <script>

    ihtml.write('</head>\n')
    ihtml.write('<body>\n')

    ihtml.write('<a href="index.html">PDAS home</a>'+ 
                    ' &gt; <a href="atmos.html">Standard Atmosphere</a>'+
                    ' &gt; BigTables\n' )
    ihtml.write('<div class="banner">' + PDAS +'</div>\n')
    ihtml.write('<header><h1>')
    ihtml.write('Tables of the U.S. Standard Atmosphere, 1976</h1></header>\n')

    ihtml.write('<!-- This is an example of a form that is never submitted.\n')
    ihtml.write('  All of the computation is done by the included\n')
    ihtml.write('  JavaScript functions. No server calculations needed. -->\n')
    ihtml.write('<form>\n')
    ihtml.write('<button id="MainButton" type="button"')
    ihtml.write(' onclick="ShowMain();">Main</button>\n')
    ihtml.write('<button id="Table1Button" type="button"')
    ihtml.write(' onclick="ShowTable1();">Table 1</button>\n')
    ihtml.write('<button id="Table2Button" type="button"')
    ihtml.write(' onclick="ShowTable2();">Table 2</button>\n')
    ihtml.write('<button id="Table3Button" type="button"')
    ihtml.write(' onclick="ShowTable3();">Table 3</button>\n')
    ihtml.write('<button id="Table4Button" type="button"')
    ihtml.write(' onclick="ShowTable4();">Table 4</button>\n')
    ihtml.write('<button id="Table5Button" type="button"')
    ihtml.write(' onclick="ShowTable5();">Table 5</button>\n')
    ihtml.write('<button id="Table6Button" type="button"')
    ihtml.write(' onclick="ShowTable6();">Table 6</button>\n')
    ihtml.write('</form>\n')

    ihtml.write('<div id="main">\n')
    ihtml.write('<p>Make tables of atmospheric properties\n')
    ihtml.write('from 0 to 1000 km.</p>\n')
    ihtml.write('<p>The tables are based on those from the reference\n')
    ihtml.write('document, Part 4.</p>\n')
 
    ihtml.write('<ol>\n')

    ihtml.write('<li>Geopotential altitude, temperature, temperature ratio,\n')
    ihtml.write('pressure, pressure ratio, density, density ratio,\n')
    ihtml.write('sound speed, garvity for geometric altitudes in kilometers.</li>\n')



    ihtml.write('<li>Geopotential altitude, dynamic viscosity, kinematic viscosity,\n')
    ihtml.write('pressure scale height, number density, mean particle speed,\n')
    ihtml.write('mean collision frequency, mean free path,\n')
    ihtml.write('thermal conductivity coefficient and molecular weight for\n')
    ihtml.write('geometric altitudes in kilometers.</li>\n')
 
    ihtml.write('<li>Non-dimensional ratios of dynamic viscosity,\n')
    ihtml.write('kinematic viscosity, number density, particle speed,\n')
    ihtml.write('mean collision frequency, mean free path, thermal conductivity,\n')
    ihtml.write('sound speed and molecular weight for altitudes in kilometers.\n')
    ihtml.write('All quantities are referenced to sea-level values.</li>\n')

    ihtml.write('<li>Geopotential altitude, temperature, temperature ratio,\n')
    ihtml.write('pressure, pressure ratio, mass density, density ratio,\n')
    ihtml.write('weight density, sound speed and gravity for geometric altitudes\n')
    ihtml.write('in thousand feet. Table entries in US customary units.</li>\n')

    ihtml.write('<li>Geopotential altitude, dynamic viscosity, kinematic viscosity,\n')
    ihtml.write('pressure scale height, number density, mean particle speed,\n')
    ihtml.write('mean collision frequency, mean free path,\n')
    ihtml.write('thermal conductivity coefficient and molecular weight for\n')
    ihtml.write('geometric altitudes in thousand feet.\n')
    ihtml.write('Table entries in US customary units.</li>\n')

    ihtml.write('<li>Non-dimensional ratios of dynamic viscosity,\n')
    ihtml.write('kinematic viscosity, number density, particle speed,\n')
    ihtml.write('mean collision frequency, mean free path, thermal conductivity,\n')
    ihtml.write('sound speed and molecular weight for altitudes in thousand feet.\n')
    ihtml.write('All quantities are referenced to sea-level values.</li>\n')

 
    ihtml.write('</ol>\n')

    ihtml.write('<p>To return to this page, press Main</p>\n')
    ihtml.write('</div>\n')

def WriteTextCell(ihtml,x):
    ihtml.write('<td class="center">' + x + '</td>\n')

def WriteIntegerCell(ihtml,k):
    ihtml.write('<td>' + str(k) + '</td>\n')

def WriteF0Cell(ihtml,x):
    ihtml.write('<td>' + '{:0.0f}'.format(x) + '</td>\n')

def WriteF1Cell(ihtml,x):
    ihtml.write('<td>' + '{:0.1f}'.format(x) + '</td>\n')

def WriteF2Cell(ihtml,x):
    ihtml.write('<td>' + '{:0.2f}'.format(x) + '</td>\n')

def WriteF3Cell(ihtml,x):
    ihtml.write('<td>' + '{:0.3f}'.format(x) + '</td>\n')

def WriteF4Cell(ihtml,x):
    ihtml.write('<td>' + '{:0.4f}'.format(x) + '</td>\n')


def WriteES4Cell(ihtml,x):
    ihtml.write('<td>' + '{:12.4E}'.format(x) + '</td>\n')

def WriteBlankCell(ihtml):
    ihtml.write('<td></td>\n')


def MakeTable1(ihtml):
    ihtml.write('<div id="table1">\n')
    ihtml.write("<h2>Table 1 Temperature, Pressure, & Density (SI units)</h2>\n")
    ihtml.write("<h3>0 to 1000 Km in steps of 5 Km</h3>\n")
    ihtml.write("<table>\n")

    for i in range(0,1001,5):
        if (i%200 == 0):
            ihtml.write("<tr>\n")
            WriteTextCell(ihtml,'Z')
            WriteTextCell(ihtml,'H')
            WriteTextCell(ihtml,'T')
            WriteTextCell(ihtml,'T/T<sub>0</sub>')
            WriteTextCell(ihtml,'p')
            WriteTextCell(ihtml,'p/p<sub>0</sub>')
            WriteTextCell(ihtml,'&rho;')
            WriteTextCell(ihtml,'&rho;/&rho;<sub>0</sub>')
            WriteTextCell(ihtml,'c')
            WriteTextCell(ihtml,'g')
            ihtml.write('</tr>\n')

        if (i==0):
            ihtml.write('<tr>\n') # only printed at top of table
            WriteTextCell(ihtml,'km')
            WriteTextCell(ihtml,'km')
            WriteTextCell(ihtml,'K')
            WriteBlankCell(ihtml)
            WriteTextCell(ihtml,'Pa')
            WriteBlankCell(ihtml)
            WriteTextCell(ihtml,'kg m<sup>-3</sup>')
            WriteBlankCell(ihtml)
            WriteTextCell(ihtml,'m/s')
            WriteTextCell(ihtml,'m s<sup>-2</sup>')
            ihtml.write('</tr>\n')

        if (i%10 == 0):
            ihtml.write('<tr class="gray">\n')   # make every other row gray
        else:
            ihtml.write('<tr>\n')

        altKm=float(i)
        hkm=altKm*REARTH/(altKm+REARTH)      # convert geometric to geopotential
        g=GZERO*(REARTH/(altKm+REARTH))**2

        (sigma, delta, theta) = Atmosphere(altKm)
 #       density=sigma*RHOZERO
        asound=ASOUNDZERO*math.sqrt(theta)

        WriteIntegerCell(ihtml,i)
        WriteF1Cell(ihtml,hkm)
        WriteF3Cell(ihtml, theta*TZERO)
        WriteF4Cell(ihtml, theta)
        WriteES4Cell(ihtml,delta*PZERO)
        WriteES4Cell(ihtml,delta)
        WriteES4Cell(ihtml,sigma*RHOZERO)
        WriteES4Cell(ihtml,sigma)
        WriteF2Cell(ihtml,asound)
        WriteF4Cell(ihtml,g)
        ihtml.write('</tr>\n')

    ihtml.write('</table>\n</div>\n')   # ----------- End of Function MakeTable1

def MakeTable2(ihtml):
    ihtml.write('<div id="table2">\n')
    ihtml.write("<h2>Table 2 - Transport Properties (SI units)</h2>\n")
    ihtml.write("<h3>0 to 1000 Km in steps of 5 km</h3>\n")
    ihtml.write("<table>\n")

    ihtml.write("<tr>\n")
    WriteTextCell(ihtml,'alt')
    WriteTextCell(ihtml,'H')
    WriteTextCell(ihtml,'Dynamic<br />viscosity')
    WriteTextCell(ihtml,'Kinematic<br />viscosity')
    WriteTextCell(ihtml,'Pressure<br />scale<br />height')
    WriteTextCell(ihtml,'Number<br />density')
    WriteTextCell(ihtml,'Particle<br />speed')
    WriteTextCell(ihtml,'Collision<br />frequency')
    WriteTextCell(ihtml,'Mean<br />free<br />path')
    WriteTextCell(ihtml,'Thermal<br />conductivity')
    WriteTextCell(ihtml,'Mol.<br />weight')
    ihtml.write('</tr>\n')

    for i in range(0,1001,5):
        if (i%(5*40) == 0):
            ihtml.write("<tr>\n")
            WriteTextCell(ihtml,'Z')
            WriteTextCell(ihtml,'H')
            WriteTextCell(ihtml,'&mu;')
            WriteTextCell(ihtml,'&eta;')
            WriteTextCell(ihtml,'H<sub>p</sub>')
            WriteTextCell(ihtml,'n')
            WriteTextCell(ihtml,'V')
            WriteTextCell(ihtml,'&nu;')
            WriteTextCell(ihtml,'L')
            WriteTextCell(ihtml,'&kappa;')
            WriteTextCell(ihtml,'M')
            ihtml.write('</tr>\n')

        if (i==0):
            ihtml.write('<tr>\n') # only printed at top of table
            WriteTextCell(ihtml,'km')
            WriteTextCell(ihtml,'km')
            WriteTextCell(ihtml,'kg m<sup>-1</sup> s<sup>-1</sup>')
            WriteTextCell(ihtml,'m<sup>2</sup>/s')
            WriteTextCell(ihtml,'m')
            WriteTextCell(ihtml,'m<sup>-3</sup>')
            WriteTextCell(ihtml,'m/s')
            WriteTextCell(ihtml,'s<sup>-1</sup>')
            WriteTextCell(ihtml,'m')
            WriteTextCell(ihtml,'W/(m K)')
            WriteBlankCell(ihtml)
            ihtml.write('</tr>\n')

        if (i%10 == 0):
            ihtml.write('<tr class="gray">\n')   # make every other row gray background
        else:
            ihtml.write('<tr>\n')

        altKm=float(i)
        hkm=altKm*REARTH/(altKm+REARTH)      # convert geometric to geopotential
        gRatio=(REARTH/(altKm+REARTH))**2
#        g=gRatio*GZERO
        (sigma,delta,theta) = Atmosphere(altKm)
        (gRatio,dynamicViscosityRatio,kinematicViscosityRatio,numberDensityRatio,
        collisionFrequencyRatio,particleSpeedRatio,meanFreePathRatio,
        thermalConductivityRatio,soundSpeedRatio, molecularWeightRatio) = TransportRatios(altKm)

        dynamicViscosity=MUZERO*dynamicViscosityRatio
        kinematicViscosity=ETAZERO*kinematicViscosityRatio
        pressureScaleHeight=PRESSURE_SCALE_HEIGHT_ZERO*theta/(gRatio*molecularWeightRatio)
        numDensity=NUM_DENSITY_ZERO*numberDensityRatio
        partSpeed=PART_SPEED_ZERO*particleSpeedRatio
        meanFreePath=FREE_PATH_ZERO*meanFreePathRatio
        collisionFrequency=partSpeed/meanFreePath
        kappa=KAPPA_ZERO*thermalConductivityRatio

        WriteIntegerCell(ihtml,i)
        WriteF1Cell(ihtml,hkm)
        WriteES4Cell(ihtml, dynamicViscosity)
        WriteES4Cell(ihtml,kinematicViscosity)
        WriteIntegerCell(ihtml,int(round(pressureScaleHeight)))
        WriteES4Cell(ihtml,numDensity)
        WriteF2Cell(ihtml,partSpeed)     # center this field
        WriteES4Cell(ihtml,collisionFrequency)
        WriteES4Cell(ihtml,meanFreePath)
        WriteES4Cell(ihtml,kappa)
        WriteF3Cell(ihtml,molecularWeightRatio*MOLWT_ZERO)
        ihtml.write('</tr>\n')

    ihtml.write('</table>\n</div>\n')   # ----------- End of Function MakeTable2

def MakeTable3(ihtml):
    ihtml.write( '<div id="table3">\n')
    ihtml.write( "<h2>Table 3 - Atmosphere Properties (non-dimensional)</h2>\n")
    ihtml.write( "<h3>0 to 1000 Km in steps of 5 km</h3>\n")
    ihtml.write( "<p>The units of alt and H are kilometers.</p>\n")
    ihtml.write( "<table>\n")

    for i in range(0,1001,5):
        if i%(5*40) == 0 :    # print every 40 lines
            ihtml.write("<tr>\n")
            WriteTextCell(ihtml,'Z')
            WriteTextCell(ihtml,'H')
            WriteTextCell(ihtml,'&mu;/&mu;<sub>0</sub>')
            WriteTextCell(ihtml,'&eta;/&eta;<sub>0</sub>')
            WriteTextCell(ihtml,'n/n<sub>0</sub>')
            WriteTextCell(ihtml,'V/V<sub>0</sub>')
            WriteTextCell(ihtml,'&nu;/&nu;<sub>0</sub>')
            WriteTextCell(ihtml,'L/L<sub>0</sub>')
            WriteTextCell(ihtml,'&kappa;/&kappa;<sub>0</sub>')
            WriteTextCell(ihtml,'c/c<sub>0</sub>')
            WriteTextCell(ihtml,'M/M<sub>0</sub>')
            ihtml.write('</tr>\n')

        if i%10 == 0:    # Make every other row with light gray background.
            ihtml.write('<tr class="gray">\n')
        else:
            ihtml.write('<tr>\n')

        altKm=float(i)
        hkm=altKm*REARTH/(altKm+REARTH)      # convert geometric to geopotential

        (sigma,delta,theta) = Atmosphere(altKm)
 #       mw=MolecularWeight(altKm)
        (gRatio,dynamicViscosityRatio,kinematicViscosityRatio,
        numberDensityRatio, collisionFrequencyRatio,particleSpeedRatio,
        meanFreePathRatio,thermalConductivityRatio,soundSpeedRatio,
        molecularWeightRatio) = TransportRatios(altKm)

        WriteIntegerCell(ihtml,i)
        WriteF1Cell(ihtml,hkm)
        WriteF4Cell(ihtml,dynamicViscosityRatio)
        WriteES4Cell(ihtml,kinematicViscosityRatio)
        WriteES4Cell(ihtml,numberDensityRatio)
        WriteF4Cell(ihtml,particleSpeedRatio)
        WriteES4Cell(ihtml,collisionFrequencyRatio)
        WriteES4Cell(ihtml,meanFreePathRatio)
        WriteF4Cell(ihtml,thermalConductivityRatio)
        WriteF4Cell(ihtml,soundSpeedRatio)
        WriteF4Cell(ihtml,molecularWeightRatio)
        ihtml.write("</tr>\n")


    ihtml.write('</table>\n</div>\n')   # ----------- End of Function MakeTable3

def MakeTable4(ihtml):
    ihtml.write( '<div id="table4">\n')
    ihtml.write('<h2>Table 4 - Temperature, Pressure, & Density (US units)</h2>\n')
    ihtml.write( '<h3>0 to 3000 Kft in steps of 10 Kft</h3>\n')
    ihtml.write( "<table>\n")

    for i in range(USMIN,USMAX+1,USDEL):
        if (i%(40*USDEL) == 0):
            ihtml.write("<tr>\n")
            WriteTextCell(ihtml,'Z')
            WriteTextCell(ihtml,'H')
            WriteTextCell(ihtml,'T')
            WriteTextCell(ihtml,'T/T<sub>0</sub>')
            WriteTextCell(ihtml,'Pressure')
            WriteTextCell(ihtml,'p/p<sub>0</sub>')
            WriteTextCell(ihtml,'Mass<br />density')
            WriteTextCell(ihtml,'&rho;/&rho;<sub>0</sub>')
            WriteTextCell(ihtml,'Weight<br />density')
            WriteTextCell(ihtml,'Sound<br />speed')
            WriteTextCell(ihtml,'Gravity')
            ihtml.write('</tr>\n')

        if (i==0):
            ihtml.write('<tr>\n') # only printed at top of table
            WriteTextCell(ihtml,'kft')
            WriteTextCell(ihtml,'kft')
            WriteTextCell(ihtml,'R')
            WriteBlankCell(ihtml)
            WriteTextCell(ihtml,'psf')
            WriteBlankCell(ihtml)
            WriteTextCell(ihtml,'slug-ft<sup>-3</sup>')
            WriteBlankCell(ihtml)
            WriteTextCell(ihtml,'lb-ft<sup>-3</sup>')
            WriteTextCell(ihtml,'ft/s')
            WriteTextCell(ihtml,'ft s<sup>-2</sup>')
            ihtml.write('</tr>\n')

        if (i%(2*USDEL) == 0):
            ihtml.write('<tr class="gray">\n')   # make every other row gray
        else:
            ihtml.write('<tr>\n')

        altKft=float(i)
        altKm=FT2METERS*altKft
        hkm=altKm*REARTH/(altKm+REARTH)      # convert geometric to geopotential
        hkft=hkm/FT2METERS
        gRatio=(REARTH/(altKm+REARTH))**2
        gus=(GZERO/FT2METERS)*gRatio
        (sigma, delta, theta) = Atmosphere(altKm)
 #       density=sigma*RHOZERO
        asound=(ASOUNDZERO/FT2METERS)*math.sqrt(theta)

        WriteIntegerCell(ihtml,i)
        WriteF1Cell(ihtml, hkft)
        WriteF3Cell(ihtml, theta*TZERO*KELVIN2RANKINE)
        WriteF4Cell(ihtml, theta)
        WriteES4Cell(ihtml,delta*PZERO/PSF2NSM)
        WriteES4Cell(ihtml,delta)
        WriteES4Cell(ihtml,sigma*RHOZERO/SCF2KCM)
        WriteES4Cell(ihtml,sigma)
        WriteES4Cell(ihtml,gus*sigma*RHOZERO/SCF2KCM)
        WriteF2Cell(ihtml,asound)
        WriteF4Cell(ihtml,gus)
        ihtml.write('</tr>\n')

    ihtml.write('</table>\n</div>\n')   # ----------- End of Function MakeTable4

def MakeTable5(ihtml):
    ihtml.write('<div id="table5">\n')
    ihtml.write('<h2>Table 5 - Transport Properties in US Customary units</h2>\n')
    ihtml.write('<h3>0 to 3000 Kft in steps of 10 Kft</h3>\n')
    ihtml.write('<table>\n')


    ihtml.write("<tr>\n")
    WriteTextCell(ihtml,'alt')
    WriteTextCell(ihtml,'H')
    WriteTextCell(ihtml,'Dynamic<br />viscosity')
    WriteTextCell(ihtml,'Kinematic<br />viscosity')
    WriteTextCell(ihtml,'Pressure<br />scale<br />height')
    WriteTextCell(ihtml,'Number<br />density')
    WriteTextCell(ihtml,'Particle<br />speed')
    WriteTextCell(ihtml,'Collision<br />frequency')
    WriteTextCell(ihtml,'Mean<br />free<br />path')
    WriteTextCell(ihtml,'Thermal<br />conductivity')
    WriteTextCell(ihtml,'Mol.<br />weight')
    ihtml.write('</tr>\n')

    for i in range(USMIN,USMAX+1,USDEL):
        if i%(40*USDEL) == 0:
            ihtml.write("<tr>\n")
            WriteTextCell(ihtml,'Z')
            WriteTextCell(ihtml,'H')
            WriteTextCell(ihtml,'&mu;')
            WriteTextCell(ihtml,'&eta;')
            WriteTextCell(ihtml,'H<sub>p</sub>')
            WriteTextCell(ihtml,'n')
            WriteTextCell(ihtml,'V')
            WriteTextCell(ihtml,'&nu;')
            WriteTextCell(ihtml,'L')
            WriteTextCell(ihtml,'&kappa;')
            WriteTextCell(ihtml,'M')
            ihtml.write('</tr>\n')

        if (i==0):
            ihtml.write('<tr>\n') # only printed at top of table
            WriteTextCell(ihtml,'kft')
            WriteTextCell(ihtml,'kft')
            WriteTextCell(ihtml,'slug ft<sup>-1</sup> s<sup>-1</sup>')
            WriteTextCell(ihtml,'ft<sup>2</sup>/s')
            WriteTextCell(ihtml,'ft')
            WriteTextCell(ihtml,'ft<sup>-3</sup>')
            WriteTextCell(ihtml,'ft/s')
            WriteTextCell(ihtml,'s<sup>-1</sup>')
            WriteTextCell(ihtml,'ft')
            WriteTextCell(ihtml,'BTU/(s ft R)')
            WriteBlankCell(ihtml)
            ihtml.write('</tr>\n')

        if (i%(2*USDEL) == 0):
            ihtml.write('<tr class="gray">\n')   # make every other row gray background
        else:
            ihtml.write('<tr>\n')

        altKft=float(i)
        altKm=FT2METERS*altKft
        hkm=altKm*REARTH/(altKm+REARTH)      # convert geometric to geopotential
        hkft=hkm/FT2METERS

        gRatio=(REARTH/(altKm+REARTH))**2
#        g=gRatio*GZERO/FT2METERS
        (sigma,delta,theta) = Atmosphere(altKm)
        (gRatio,dynamicViscosityRatio,kinematicViscosityRatio,numberDensityRatio,
        collisionFrequencyRatio,particleSpeedRatio,meanFreePathRatio,
        thermalConductivityRatio,soundSpeedRatio, molecularWeightRatio) = TransportRatios(altKm)

        dynamicViscosity=(MUZERO*FT2METERS/SLUG2KG)*dynamicViscosityRatio
        kinematicViscosity=(ETAZERO/FT2METERS**2)*kinematicViscosityRatio
        pressureScaleHeight=theta/(gRatio*molecularWeightRatio)
        pressureScaleHeight=(PRESSURE_SCALE_HEIGHT_ZERO/FT2METERS)*pressureScaleHeight
        numDensity=(NUM_DENSITY_ZERO*FT2METERS**3)*numberDensityRatio
        partSpeed=(PART_SPEED_ZERO/FT2METERS)*particleSpeedRatio
        meanFreePath=(FREE_PATH_ZERO/FT2METERS)*meanFreePathRatio
        collisionFrequency=partSpeed/meanFreePath
        kappa=(KAPPA_ZERO*FT2METERS/(KELVIN2RANKINE*BTU2JOULE))*thermalConductivityRatio

        WriteIntegerCell(ihtml,i)
        WriteF1Cell(ihtml,hkft)
        WriteES4Cell(ihtml, dynamicViscosity)
        WriteES4Cell(ihtml,kinematicViscosity)
        WriteIntegerCell(ihtml,int(round(pressureScaleHeight)))
        WriteES4Cell(ihtml,numDensity)
        WriteF2Cell(ihtml,partSpeed)     # center this field
        WriteES4Cell(ihtml,collisionFrequency)
        WriteES4Cell(ihtml,meanFreePath)
        WriteES4Cell(ihtml,kappa)
        WriteF3Cell(ihtml,molecularWeightRatio*MOLWT_ZERO)
        ihtml.write('</tr>\n')

    ihtml.write('</table>\n</div>\n')   # ----------- End of Function MakeTable5

def MakeTable6(ihtml):
    ihtml.write( '<div id="table6">\n')
    ihtml.write( "<h2>Table 6 - Transport Properties (Non-dimensional)</h2>\n")
    ihtml.write( "<h3>0 to 3000 Kft in steps of 10 Kft</h3>\n")
    ihtml.write( "<p>The units of alt and H are thousands of feet.</p>\n")
    ihtml.write( "<table>\n")


    for i in range(USMIN,USMAX+1,USDEL):
        if i%(40*USDEL) == 0 :
            ihtml.write("<tr>\n")
            WriteTextCell(ihtml,'Z')
            WriteTextCell(ihtml,'H')
            WriteTextCell(ihtml,'&mu;/&mu;<sub>0</sub>')
            WriteTextCell(ihtml,'&eta;/&eta;<sub>0</sub>')
            WriteTextCell(ihtml,'n/n<sub>0</sub>')
            WriteTextCell(ihtml,'V/V<sub>0</sub>')
            WriteTextCell(ihtml,'&nu;/&nu;<sub>0</sub>')
            WriteTextCell(ihtml,'L/L<sub>0</sub>')
            WriteTextCell(ihtml,'&kappa;/&kappa;<sub>0</sub>')
            WriteTextCell(ihtml,'c/c<sub>0</sub>')
            WriteTextCell(ihtml,'M/M<sub>0</sub>')
            ihtml.write('</tr>\n')

        if i%(2*USDEL) == 0:    # Make every other row with light gray background.
            ihtml.write('<tr class="gray">\n')
        else:
            ihtml.write('<tr>\n')

        altKft=float(i)
        altKm=altKft*FT2METERS
        hkm=altKm*REARTH/(altKm+REARTH)      # convert geometric to geopotential
        hkft=hkm/FT2METERS

#        (sigma,delta,theta) = Atmosphere(altKm)
#        mw=MolecularWeight(altKm)
        (gRatio,dynamicViscosityRatio,kinematicViscosityRatio,
        numberDensityRatio, collisionFrequencyRatio,particleSpeedRatio,
        meanFreePathRatio,thermalConductivityRatio,soundSpeedRatio,
        molecularWeightRatio)= TransportRatios(altKm)

        WriteIntegerCell(ihtml,i)
        WriteF1Cell(ihtml,hkft)
        WriteF4Cell(ihtml,dynamicViscosityRatio)
        WriteES4Cell(ihtml,kinematicViscosityRatio)
        WriteES4Cell(ihtml,numberDensityRatio)
        WriteF4Cell(ihtml,particleSpeedRatio)
        WriteES4Cell(ihtml,collisionFrequencyRatio)
        WriteES4Cell(ihtml,meanFreePathRatio)
        WriteF4Cell(ihtml,thermalConductivityRatio)
        WriteF4Cell(ihtml,soundSpeedRatio)
        WriteF4Cell(ihtml,molecularWeightRatio)
        ihtml.write("</tr>\n")

    ihtml.write('</table>\n</div>\n')   # ----------- End of function MakeTable6

# ------------------------------------------------------------------------------
def MakeFooter(ihtml,updater):
    from datetime import date 
    today = date.today()
    dateString = today.strftime(" %d %B %Y ")
    ihtml.write('<footer>\n')
    ihtml.write('Last updated: ' + dateString + ' by ' + updater +',\n')
    ihtml.write('<a href="mailto:pdaerowebmaster@gmail.com">')
    ihtml.write('pdaerowebmaster AT gmail DOT com</a>\n')
    ihtml.write('</footer>\n')
    ihtml.write('<div class="banner">' + PDAS + '</div>\n')
    ihtml.write('<a href="index.html">PDAS home</a>' + 
                    ' &gt; <a href="atmos.html">Standard Atmosphere</a>' +
                    ' &gt; BigTables\n' )

    ihtml.write('</body>\n')
    ihtml.write('</html>\n')   # -------------------- End of function MakeFooter

# ------------------------------------------------------------------------------
def main():
    startTime=time.time()
    print ("Executing ", sys.argv[0])
    print (greeting)
    print ("Version " + version)
    print (author)
    if modifier != "":
        print ("Modified by ", modifier)

    ihtml = open('bigtables.html', 'w')
    MakeMainPage(ihtml)
    MakeTable1(ihtml)
    MakeTable2(ihtml)
    MakeTable3(ihtml)
    MakeTable4(ihtml)
    MakeTable5(ihtml)
    MakeTable6(ihtml)

    if modifier == "":
        MakeFooter(ihtml,author)
    else:
        MakeFooter(ihtml,modifier)
    ihtml.close()
    print (farewell)
    print (finalmess)
    endTime=time.time()

    print(str(round(1000*(endTime-startTime))) + " ms\n")


main()

