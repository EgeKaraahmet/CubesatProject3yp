#!/usr/bin/env python3 
#!/usr/bin/python3

""""

PROGRAM Press;                                           { \atmos\press.py }
(* ----------------------------------------------------------------------- *)
(* PURPOSE - Compute the pressure ratio at the altitudes that are the      *)
(*    boundaries between layers of the 1976 Standard Atmosphere            *)
(* AUTHOR  - Ralph L. Carmichael, Public Domain Aeronautical Software      *)

# Running this program should produce the following output: 

Press - compute pressure at layer boundaries (1976 std)
Ralph L. Carmichael, Public Domain Aeronautical Software
Hydrostatic constant =  34.163194736310366    K/km


   z   temperature        p/p0          density/density0
   0.0   288.15000    1.0000000000000    1.0000000000000
  11.0   216.65000    0.2233611050922    0.2970759401445
  20.0   216.65000    0.0540329501078    0.0718651953546
  32.0   228.65000    0.0085666783593    0.0107959255160
  47.0   270.65000    0.0010945601338    0.0011653334659
  51.0   270.65000    0.0006606353133    0.0007033514337
  71.0   214.65000    0.0000390468337    0.0000524171681
  84.9   186.94600    0.0000036850095    0.0000056799049
Normal termination of press, version 1.5 (19Aug2018)


(* REVISION HISTORY                                                        *)
(*   DATE  VERS PERSON  STATEMENT OF CHANGES                               *)
(* 10Aug95  1.0   RLC   Original coding                                    *)
(*  3Nov95  1.1   RLC   Added column headers                               *)
(* 21Nov96  1.2   RLC   Writes to a file (easier for Delphi)               *)
(* 28Nov00  1.3   RLC   Made it a Delphi Console Application               *)
(* 25Jun12  1.4   RLC   Commented out the uses SysUtils statement          *)
   19Aug18  1.5   RLC   Rewrote in Python
(* ----------------------------------------------------------------------- *)

"""
import math

GREETING  = 'Press - compute pressure at layer boundaries (1976 std)'
AUTHOR    = 'Ralph L. Carmichael, Public Domain Aeronautical Software'
VERSION   = '1.5 (19Aug2018)'
FINALMESSAGE = 'Normal termination of press, version ' + VERSION

GRAVITY = 9.80665      #  accel. of gravity 
MOL_WT  = 28.9644      #  molecular weight of air
R_GAS   = 8314.32      #  gas constant N m / (kmol K)

GMR = 1000 * GRAVITY * MOL_WT / R_GAS
# NOTE - Observe the factor of 1000 in computing GMR. Without it, GMR would have 
#   units of kelvins per meter. But the temperature gradient in LowerAtmosphere
#   is given in kelvins per kilometer.   GMR=34.163195 kelvins/km
 
print(GREETING)
print(AUTHOR)
print('Hydrostatic constant = ',GMR,'   K/km\n\n')
print('   z   temperature        p/p0          density/density0')

t0=288.15
p0=1.0
s=1.0   # density ratio
z0=0.0
print('{:6.1f}'.format(z0),'{:11.5f}'.format(t0),
      '{:18.13f}'.format(p0), '{:18.13f}'.format(s))

z11=11.0
tgrad=-6.5
deltaz=z11-z0
t11 = t0 + tgrad*deltaz
p11 = p0*math.exp(math.log(t0/t11)*GMR/tgrad)
s=p11*t0/t11
print('{:6.1f}'.format(z11),'{:11.5f}'.format(t11),
      '{:18.13f}'.format(p11), '{:18.13}'.format(s))
   
z20=20.0
tgrad=0.0   
deltaz = z20-z11
t20=t11 + tgrad*deltaz
p20=p11*math.exp(-GMR*deltaz/t11)          # special form when tgrad=0 
s=p20*t0/t20
print('{:6.1f}'.format(z20),'{:11.5f}'.format(t20),
      '{:18.13f}'.format(p20), '{:18.13f}'.format(s))

z32=32.0
tgrad=1.0
deltaz = z32-z20
t32=t20+tgrad*deltaz
p32=p20*math.exp(math.log(t20/t32)*GMR/tgrad)
s=p32*t0/t32
print('{:6.1f}'.format(z32),'{:11.5f}'.format(t32),
      '{:18.13f}'.format(p32), '{:18.13f}'.format(s))

z47=47.0
tgrad=2.8
deltaz = z47-z32
t47=t32+tgrad*deltaz
p47=p32*math.exp(math.log(t32/t47)*GMR/tgrad)
s=p47*t0/t47
print('{:6.1f}'.format(z47),'{:11.5f}'.format(t47),
      '{:18.13f}'.format(p47), '{:18.13f}'.format(s))

z51=51.0
tgrad=0.0
deltaz = z51-z47
t51=t47 + tgrad*deltaz
p51=p47*math.exp(-GMR*deltaz/t47)            #  special form when tgrad=0 
s=p51*t0/t51
print('{:6.1f}'.format(z51),'{:11.5f}'.format(t51),
      '{:18.13f}'.format(p51), '{:18.13f}'.format(s))

z71=71.0
tgrad=-2.8
deltaz=z71-z51
t71=t51+tgrad*deltaz
p71=p51*math.exp(math.log(t51/t71)*GMR/tgrad)
s=p71*t0/t71
print('{:6.1f}'.format(z71),'{:11.5f}'.format(t71),
      '{:18.13f}'.format(p71), '{:18.13f}'.format(s))

z85=84.852
tgrad=-2.0
deltaz=z85-z71
t85=t71+tgrad*deltaz
p85=p71*math.exp(math.log(t71/t85)*GMR/tgrad)
s=p85*t0/t85
print('{:6.1f}'.format(z85),'{:11.5f}'.format(t85),
      '{:18.13f}'.format(p85), '{:18.13f}'.format(s))

print(FINALMESSAGE)

