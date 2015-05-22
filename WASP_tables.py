# -*- coding: utf-8 -*-
import csv
from planet_tables import *
from planet_tables import transit_equations

def add_wasp_planets(p_dict,file_name,search_radius,maglim,declim):

  raw_planets = load_wasp_csv(file_name)

  for p in raw_planets:
    
    ra, dec = ra_dec_from_id(p['obj_id'])
    print ra, dec
    #WHAT IS WASP TIME EPOCH??!

    try:
      try:
	dur = float(p['width'])
	dur_p = float(p['width_p'])
	dur_m = float(p['width_m'])
      except:
	dur = float(p['width_sec'])
	dur_p = float(p['width_sec_p'])
	dur_m = float(p['width_sec_m'])
    except:
      
      try:
	dur = transit_equations.estimate_duration(float(p['p_radius']),float(p['s_radius']),float(p['period']),float(p['separation']),float(p['impact_par']))
	dur_p = 0
	dur_m = 0
      except:
	dur = transit_equations.estimate_duration(float(p['p_radius']),float(p['s_radius']),float(p['period']),float(p['separation']),0)
	dur_p = 0
	dur_m = 0

    t_offset = 2450000.0

    for key in p:
      if p[key] == 'NULL':
	p[key] = 0.0

    vals = [p['nickname'].rstrip('b')+' b',p['nickname'].rstrip('b'),ra,dec,float(p['p_mass']),float(p['p_radius']),float(p['p_temp']),float(p['s_radius']),'',float(p['s_temp']),float(p['separation']),float(p['epoch'])+t_offset,float(p['epoch_p']),float(p['epoch_m']),dur,dur_p,dur_m,float(p['period']),float(p['period_p']),float(p['period_m'])]
    p_dict = add_new_planet(p_dict,vals,search_radius)


  if declim == None:
    is_sane = [(p_dict['st_vj'] < maglim) & (p_dict['st_vj'] > -1)]
  else:
    is_sane = [(p_dict['dec'] < declim) & (p_dict['st_vj'] < maglim) & (p_dict['st_vj'] > -1)]

  for i in range(0,len(p_dict.keys())):
    p_dict[p_dict.keys()[i]] = p_dict[p_dict.keys()[i]][is_sane]

  return p_dict

def ra_dec_from_id(obj_id):
    rastring = obj_id[8:17]
    decstring = obj_id[17:]

    ra_h = float(rastring[:2])
    ra_m = float(rastring[2:4])
    ra_s = float(rastring[4:])

    dec_d = float(decstring[:3])
    dec_m = float(decstring[3:5])
    dec_s = float(decstring[5:])

    ra = 360*(ra_h + (ra_m/60.0) + (ra_s/3600))/24.0


    sign = float(decstring[:1] + '1')

    dec = (sign*dec_d + (dec_m/60.0) + (dec_s/3600))*sign

    return ra, dec

def load_wasp_csv(file_name):

  planets = []
  with open(file_name,mode='rU') as infile:
    dict_reader = csv.DictReader(infile)
    for row in dict_reader:
      planets += [row]
  return planets