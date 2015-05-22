# -*- coding: utf-8 -*-
import numpy as np

# transits and occultations winn et al

def estimate_duration(R_planet,R_star,period,semi_major,b):

  R_jup = 71492000.0
  R_sun = 69550000.0
  AU = 149597870700.0

  R_planet = R_planet*R_jup
  R_star = R_star*R_sun
  semi_major = semi_major*AU

  i = np.arccos(b*R_star/semi_major)

  k = R_planet/R_star
  
  duration = (period/np.pi)*np.arcsin((R_star/semi_major)*np.sqrt((1.0+k)**2 - b**2)/np.sin(i))
  
  return duration