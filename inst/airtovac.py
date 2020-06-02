def airtovac(wave_air):
   """
taken from idl astrolib
;+
; NAME:
;       AIRTOVAC
; PURPOSE:
;       Convert air wavelengths to vacuum wavelengths 
; EXPLANATION:
;       Wavelengths are corrected for the index of refraction of air under 
;       standard conditions.  Wavelength values below 2000 A will not be 
;       altered.  Uses relation of Ciddor (1996).
;
; CALLING SEQUENCE:
;       AIRTOVAC, WAVE_AIR, [ WAVE_VAC]
;
; INPUT/OUTPUT:
;       WAVE_AIR - Wavelength in Angstroms, scalar or vector
;               If this is the only parameter supplied, it will be updated on
;               output to contain double precision vacuum wavelength(s). 
; OPTIONAL OUTPUT:
;        WAVE_VAC - Vacuum wavelength in Angstroms, same number of elements as
;                 WAVE_AIR, double precision
;
; EXAMPLE:
;       If the air wavelength is  W = 6056.125 (a Krypton line), then 
;       AIRTOVAC, W yields an vacuum wavelength of W = 6057.8019
;
; METHOD:
;	Formula from Ciddor 1996, Applied Optics 62, 958
;
; NOTES:
;       Take care within 1 A of 2000 A.   Wavelengths below 2000 A *in air* are
;       not altered.
; REVISION HISTORY
;       Written W. Landsman                November 1991
;       Use Ciddor (1996) formula for better accuracy in the infrared 
;           Added optional output vector, W Landsman Mar 2011
;       Iterate for better precision W.L./D. Schlegel  Mar 2011
;-
   """

   wave_vac = wave_air * 1.0
   g = wave_vac > 2000     #Only modify above 2000 A

   if np.sum(g):
      for iter in [0, 1]:
         if isinstance(g, np.ndarray):
            sigma2 = (1e4/wave_vac[g])**2.     #Convert to wavenumber squared
            # Compute conversion factor
            fact = 1. + 5.792105e-2 / (238.0185 - sigma2) + \
                               1.67917e-3 / (57.362 - sigma2)
            wave_vac[g] = wave_air[g] * fact              #Convert Wavelength
         else: # scalar version
            sigma2 = (1e4/wave_vac)**2.     #Convert to wavenumber squared
            # Compute conversion factor
            fact = 1. + 5.792105e-2 / (238.0185 - sigma2) + \
                               1.67917e-3 / (57.362 - sigma2)
            wave_vac = wave_air * fact              #Convert Wavelength

   return wave_vac

