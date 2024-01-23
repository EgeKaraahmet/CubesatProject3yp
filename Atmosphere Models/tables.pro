;+
; NAME:
;   USSA
;
; PURPOSE:
;   Compute the properties of the 1976 US Standard Atmosphere to 86
;   km.
;
; CALLING SEQUENCE:
;   ussa, alt, sigma, delta, theta
;
; INPUTS:
;   alt -> geometric altitude [km]
;
; OUTPUTS:
;   sigma -> density/sea-level standard density (1.347 kg/m3); or
;      density in kg/m3 if absolute keyword is set.
;
;   delta -> pressure/sea-level standard pressure (101325. Pa); or
;      pressure in Pa if absolute keyword is set.
;
;   theta -> temperature/sea-level standard temperature (288.15 K);
;      or temperature in K if absolute keyword is set.
;
; KEYWORDS:
;   absolute (boolean) -> return absolute values for density,
;   pressure, and temperature, based on the tabulated default values
;
; NOTES:
;   (1) The US standard atmosphere is based on measurements and model
;   assumptions valid for northern midlatitudes. 
;   (2) This program produces correct values up to 86 km. Density will
;   be not too far off above 86 km, the terms pressure and temperature
;   are not well defined above this altitude.
;
; AUTHOR:
;   Ralph Carmichael, Public Domain Aeronautical Software
;   Martin Schultz, MPI for Meteorology, IDL translation
;-

PRO ussa, alt, sigma, delta, theta, absolute=absolute

   nalt = n_elements(alt)
   IF nalt EQ 0 THEN BEGIN 
      message, 'Usage: ussa, alt, sigma, delta, theta', /CONTINUE
      print,'    alt: geometric altitude in km'
      print,'  sigma: density/sea-level density'
      print,'  delta: pressure/sea-level pressure'
      print,'  theta: temperature/sea-level temperature'
      return
   ENDIF 

   ;; local constants
   REARTH = 6356.766d0           ; radius of the earth
   GMR = 34.163195d0            ; hydrostatic constant
   
   ;; local arrays (US 1976 standard atmosphere)
   htab = double([ 0.0, 11.0, 20.0, 32.0, 47.0, 51.0, 71.0, 84.852 ])
   ttab = double([288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65, 186.946])
   ptab = double([1.0, 2.233611E-1, 5.403295E-2, 8.5666784E-3, 1.0945601E-3, $
                  6.6063531E-4, 3.9046834E-5, 3.68501E-6])
   gtab = double([-6.5, 0.0, 1.0, 2.8, 0.0, -2.8, -2.0, 0.0])
   
   ;; convert geometric to geopotential altitude
   h = alt*REARTH/(REARTH+alt)

   ;; set up result arrays
   sigma = fltarr(nalt)
   delta = fltarr(nalt)
   theta = fltarr(nalt)

   ;; search of altitude
   FOR i=0L, nalt-1 DO BEGIN 
      w = max(where(htab LE (h[i] > 0.)))
      tgrad = gtab[w]
      tbase = ttab[w]
      deltah = (h[i]-htab[w]) > 0.
      tlocal = tbase+tgrad*deltah
      theta[i] = tlocal/ttab[0]
      IF (abs(tgrad) LT 1.d-6) THEN BEGIN 
         delta[i] = ptab[w]*exp(-GMR*deltah/tbase)
      ENDIF ELSE BEGIN 
         delta[i] = ptab[w]*(tbase/tlocal)^(GMR/tgrad)
      ENDELSE 
      sigma[i] = delta[i]/theta[i]
   ENDFOR 

   IF keyword_set(absolute) THEN BEGIN 
      theta = theta*ttab[0]
      delta = delta*101325.
      sigma = sigma*1.347
   ENDIF 

END 

