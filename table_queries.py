#!/usr/bin/python
# -*- coding: utf-8 -*-
import os
import numpy as np
from pylab import *
import csv
import urllib2
import StringIO
from astroquery.vizier import Vizier
from astroquery.simbad import Simbad
import astropy.units as u
import astropy.coordinates as coord
import pickle
from sys import argv

# Get a file-like object for the Python Web site's home page.

def example_script():

# a (messy) example workflow script

  fname = argv[-1]


  declim = 12
  declim = None
  maglim = 12
  search_radius = 5.0

  try:
    p_dict = pickle.load(open(fname+'.p','rb'))
  except:
    create_dict(fname,declim,maglim,search_radius)
    p_dict = pickle.load(open(fname+'.p','rb'))

  if declim == None:
    is_sane = [(p_dict['st_vj'] < maglim) & (p_dict['st_vj'] > 0)]
  else:
    is_sane = [(p_dict['dec'] < declim) & (p_dict['st_vj'] < maglim) & (p_dict['st_vj'] > 0)]

  for i in range(0,len(p_dict.keys())):
    p_dict[p_dict.keys()[i]] = p_dict[p_dict.keys()[i]][is_sane]

  est_ugriz(p_dict,'WASP-12 b')

  print_comparisons(p_dict)

  p_dict['pl_trandur'][p_dict['pl_name'] == 'GJ 3470 b'] = 0.0673
  p_dict['pl_trandur'][p_dict['pl_name'] == 'CoRoT-1 b'] = 0.09652777777
  p_dict['pl_trandur'][p_dict['pl_name'] == 'CoRoT-4 b'] = 0.18402777777

  p_dict['pl_trandur'][p_dict['pl_name'] == 'WASP-2 b'] = 0.03138
  p_dict['pl_eqt'][p_dict['pl_name'] == 'WASP-2 b'] = 1300

#  kill_list = ['CoRoT-12 b','CoRoT-19 b','CoRoT-18 b','CoRoT-13 b','CoRoT-14 b','CoRoT-20 b','CoRoT-5 b','HAT-P-43 b','CoRoT-7 b','WASP-39 b','WASP-17 b','WASP-52 b','WASP-31 b','GJ 1214 b','SWEEPS-4 b','SWEEPS-11 b','WASP-6 b','HAT-P-26 b','WASP-25 b','WASP-55 b','WASP-34 b','HAT-P-41 b','OGLE-TR-10 b','OGLE-TR-111 b','WASP-103 b','XO-1 b','CoRoT-16 b','OGLE-TR-113 b','WASP-56 b','GJ 436 b','OGLE-TR-056 b','WASP-45 b','HATS-3 b','WASP-26 b','WASP-77 A b','HAT-P-23 b','WASP-43 b','OGLE-TR-132 b','OGLE-TR-211 b','Qatar-2 b','WASP-5 b','HATS-1 b','OGLE-TR-182 b','Lupus-TR-3 b','WASP-47 b','CoRoT-8 b','WASP-37 b','WASP-59 b','OGLE2-TR-L9 b','POTS-1 b','WASP-72 b','WASP-71 b','CoRoT-6 b','CoRoT-10 b','WASP-14 b','CoRoT-17 b','CoRoT-23 b','CoRoT-27 b','55 Cnc e','HD 209458 b','HD 189733 b','WASP-124b','WASP-96b']
#  kill_list = ['WASP-16 b' 'WASP-18 b' 'WASP-24 b' 'WASP-35 b' 'WASP-38 b' 'WASP-41 b' 'WASP-44 b' 'WASP-66 b' 'WASP-79 b' 'CoRoT-1 b' 'CoRoT-3 b' 'CoRoT-4 b' 'CoRoT-9 b' 'CoRoT-11 b' 'HAT-P-25 b' 'HAT-P-30 b' 'HAT-P-42 b' 'HAT-P-49 b' 'HATS-2 b']

# planets that need to be killed outright for having no matching comparison star in distance
  kill_list = ['WASP-17 b','WASP-52 b','Qatar-1 b','TrES-3 b','HATS-3 b','WASP-14 b','WASP-38 b']
  kill_list += ['WASP-31 b', 'WASP-55 b', 'HATS-5 b', 'WASP-72 b', 'WASP-37 b']

# planets that need to be killed for failing to have a comparison star brighter than 14th
  kill_list += ['WASP-79 b','WASP-79 b secondary','WASP-39 b','WASP-25 b','WASP-71 b','WASP-5 b','CoRoT-2 b']

# planets killed for having too poor a S/N


  kill_list += ['CoRoT-7 b','WASP-18 b','WASP-73 b','WASP-8 b','WASP-32 b','WASP-66 b','WASP-61 b','WASP-68 b','WASP-47 b','HAT-P-42 b','WASP-65 b','HATS-1 b','WASP-43 b']
  p_dict = kill_targets(p_dict, kill_list)

  secondary_list = ['WASP-18 b','WASP-103 b','WASP-19 b','WASP-43 b','WASP-4 b','WASP-77 A b','WASP-54 b','WASP-78 b','WASP-46 b','WASP-64 b']
  p_dict = add_secondary(p_dict, secondary_list)

  print ' '
  i = 0
  for planet in p_dict['pl_name']:
    if planet == 'WASP-103 b':
      print 'INVESTIGATING', planet
      for key in p_dict.keys():
	print key,p_dict[key][i]
    i += 1

  print ' '
  i = 0
  for planet in p_dict['pl_name']:
    if planet == 'WASP-103 b secondary':
      print 'INVESTIGATING', planet
      for key in p_dict.keys():
	print key,p_dict[key][i]
    i += 1

  p_dict = calc_eqtemp(p_dict)
  p_dict = calc_signal(p_dict)

  pickle.dump(p_dict,open(fname+'.p','wb'))

#  keep_list = ['WASP-4 b','WASP-19 b','WASP-22 b','WASP-23 b','WASP-36 b','WASP-46 b','WASP-49 b','WASP-50 b','WASP-57 b','WASP-61 b','WASP-62 b','WASP-64 b','WASP-65 b','WASP-78 b','HAT-P-39 b','HAT-P-24 b','HAT-P-27 b','GJ 3470 b']

#  for target in keep_list:
#    print target, p_dict['pl_trandur'][p_dict['pl_name'] == target]

#  kill_list = []

#  for name in p_dict['pl_name']:
#    if any(array(keep_list) == name):
#      name
#    else:
#      p_dict = kill_targets(p_dict,[name])

#  vals = ['pl_hostname','pl_hostname','ra','dec','pl_massj','pl_radj','pl_eqt','st_rad','st_vj','st_teff','pl_orbsmax','pl_tranmid','pl_tranmiderr1','pl_tranmiderr2','pl_trandur','pl_trandurerr1','pl_trandurerr2','pl_orbper','pl_orbpererr1','pl_orbpererr2']
  #vals = ['WASP-83b','WASP-83',190.152125,-19.2842777778,0.260,1.038,0.0,1.095,'',5480,0.053,2453860.8637,0,0,0.12185,0,0,4.971161,0,0]
  #p_dict = add_new_planet(p_dict,vals)

  #vals = ['WASP-88b','WASP-88',309.51125,-48.462,0.560,1.700,0.0,2.08,'',6430,0.064310,2454626.8780,0,0,0.2603,0,0,4.954114,0,0]
  #p_dict = add_new_planet(p_dict,vals)

  #vals = ['WASP-96b','WASP-96',1.04641666667,-47.3606111111,0.48,1.20,0.0,1.05,'',5540,0.045300,2455353.7909,0,0,0.07937,0,0,3.425317,0,0]
  #p_dict = add_new_planet(p_dict,vals)

  #vals = ['WASP-101b','WASP-101',98.3510833333,-23.4861666667,0.50,1.40,0.0,1.29,'',6400,0.050600,2454848.7294,0,0,0.11318,0,0,3.585725,0,0]
  #p_dict = add_new_planet(p_dict,vals)

  #vals = ['WASP-110b','WASP-110',305.873125,-44.0584166667,0.646,1.311,0.0,0.949,'',5690,0.047000,2453861.6195,0,0,0.1106,0,0,3.778422,0,0]
  #p_dict = add_new_planet(p_dict,vals)

  #vals = ['WASP-118b','WASP-118',19.5505,2.70283333333,0.380,1.249,0.0,1.513,'',6350,0.050000,2454663.6427,0,0,0.20169,0,0,4.046042,0,0]
  #p_dict = add_new_planet(p_dict,vals)

  #vals = ['WASP-124b','WASP-124',332.714291667,-30.7495277778,0.690,1.715,0.0,1.230,'',6100,0.047700,2453865.0280,0,0,0.07557,0,0,3.372657,0,0]
  #p_dict = add_new_planet(p_dict,vals)

  run_ephem(p_dict)

  print_result(p_dict)
  
  quit()

  run_ephem(p_dict)

  print_result(p_dict)

  quit()
  is_sane = [(p_dict['eqt'] != 0) & (p_dict['pl_massj'] != 0) & (p_dict['pl_radj'] != 0) & (p_dict['st_rad'] != 0) & (p_dict['dec'] < 12.0) & (p_dict['ra'] > 50.0) & (p_dict['ra'] < 145.0)]

  for i in range(0,len(p_dict.keys())):
    p_dict[p_dict.keys()[i]] = array(p_dict[p_dict.keys()[i]])[is_sane]

  targets = ['WASP-2 b','GJ 1214 b','CoRoT-3 b','SWEEPS-4 b','SWEEPS-11 b']
  targets = ['WASP-23 b','WASP-22 b','CoRoT-1 b']
  for target in targets:
    print_target(target,p_dict)


#  print_positions(p_dict)
#  print_comparisons(p_dict)

def est_ugriz(p_dict,p_name):
  print p_dict

def add_new_planet(p_dict,vals,search_radius):

  # if you don't know the planets vband magnitude, put v_j as '' and the code will find one using vizier.

  keys = array(['pl_name','pl_hostname','ra','dec','pl_massj','pl_radj','pl_eqt','st_rad','st_vj','st_teff','pl_orbsmax','pl_tranmid','pl_tranmiderr1','pl_tranmiderr2','pl_trandur','pl_trandurerr1','pl_trandurerr2','pl_orbper','pl_orbpererr1','pl_orbpererr2'])
  new_dict = {}
  for i in range(0,len(keys)):
    new_dict[keys[i]] = array([vals[i]])
    print keys[i], new_dict[keys[i]]
  new_dict = calc_eqtemp(new_dict)
  for i in range(0,len(new_dict.keys())):
    new_dict[new_dict.keys()[i]] = new_dict[new_dict.keys()[i]]
  new_dict = vizier_find_mags(new_dict,search_radius)

  try:
    if new_dict['st_vj'] == '':
      new_dict['st_vj'] = new_dict['nearby_objects'][0]['Vmag'][argmin(new_dict['nearby_objects'][0]['sep'])]
  except:
    #err... not really sure what to do if that fails... set it to an obviously wrong value that doesn't break the code, for now.
    new_dict['st_vj'] = 0
 
  new_dict = calc_signal(new_dict)
  for key in new_dict.keys():
    p_dict[key] = append(p_dict[key],new_dict[key])

  return p_dict

def print_target(target,p_dict):
  match = target == p_dict['pl_name']
  print match
  for key in p_dict.keys():
    print key, array(p_dict[key])[match]


def run_ephem(p_dict):

  out_start = open('planets_start','wb')
  out_mid = open('planets_mid','wb')
  out_end = open('planets_end','wb')

  tlen = 0.75

  is_sane = [(p_dict['pl_tranmid'] != 0)]

  for i in range(0,len(p_dict.keys())):
    p_dict[p_dict.keys()[i]] = p_dict[p_dict.keys()[i]][is_sane]


  for i in range(0,len(p_dict['pl_name'])):
    out_start.write(p_dict['pl_name'][i] + '\n')
    out_mid.write(p_dict['pl_name'][i] + '\n')
    out_end.write(p_dict['pl_name'][i] + '\n')

    out_start.write(ra_to_sex(p_dict['ra'][i]) + ' ' + dec_to_sex(p_dict['dec'][i])+ '\n')
    out_mid.write(ra_to_sex(p_dict['ra'][i]) + ' ' + dec_to_sex(p_dict['dec'][i])+ '\n')
    out_end.write(ra_to_sex(p_dict['ra'][i]) + ' ' + dec_to_sex(p_dict['dec'][i])+ '\n')

    dur = p_dict['pl_trandur'][i]
    
    out_start.write('HJD linear' + ' ' + str(p_dict['pl_tranmid'][i] - tlen*dur) + ' ' + str((abs(p_dict['pl_tranmiderr1'][i])+abs(p_dict['pl_tranmiderr2'][i]))/2.0) + ' ' + str(p_dict['pl_orbper'][i]) + ' ' + str((abs(p_dict['pl_orbpererr1'][i])+abs(p_dict['pl_orbpererr2'][i]))/2.0)+ '\n')
    out_mid.write('HJD linear' + ' ' + str(p_dict['pl_tranmid'][i]) + ' ' + str((abs(p_dict['pl_tranmiderr1'][i])+abs(p_dict['pl_tranmiderr2'][i]))/2.0) + ' ' + str(p_dict['pl_orbper'][i]) + ' ' + str((abs(p_dict['pl_orbpererr1'][i])+abs(p_dict['pl_orbpererr2'][i]))/2.0)+ '\n')
    out_end.write('HJD linear' + ' ' + str(p_dict['pl_tranmid'][i] + tlen*dur) + ' ' + str((abs(p_dict['pl_tranmiderr1'][i])+abs(p_dict['pl_tranmiderr2'][i]))/2.0) + ' ' + str(p_dict['pl_orbper'][i]) + ' ' + str((abs(p_dict['pl_orbpererr1'][i])+abs(p_dict['pl_orbpererr2'][i]))/2.0)+ '\n')

def ra_to_sex(deg):
  hours = 24*(deg / 360.0)
  minutes = (hours - floor(hours))*60
  seconds = (minutes - floor(minutes))*60

  hours = str(int(floor(hours)))
  minutes = str(int(floor(minutes)))

  if len(hours) < 2:
    hours = '0' + hours

  if len(minutes) < 2:
    minutes = '0' + minutes

  sex = hours +' '+ minutes +' '+ str(seconds)
  return sex

def dec_to_sex(deg):
  minutes = (deg - floor(deg))*60
  seconds = (minutes - floor(minutes))*60

  deg = str(int(floor(deg)))
  minutes = str(int(floor(minutes)))

  if deg[0] != '-':
    deg = '+'+deg

  if len(deg) < 3:
    deg = deg[0] + '0' + deg[1]

  if len(minutes) < 2:
    minutes = '0' + minutes

  sex = deg +' '+ minutes +' '+ str(seconds)
  return sex

def add_columns(p_dict,add_cols):
  Columns='pl_name,'+add_cols
  p_dict_new = grab_table(Columns,Transit=True)
  return p_dict

def kill_targets(p_dict, kill_list):
  
  for string in kill_list:
    i = 0
    for planet in p_dict['pl_name']:
      if string == planet:
	print 'KILLING', planet
	for key in p_dict.keys():
	  p_dict[key] = delete(p_dict[key],i)
	i -= 1
      i += 1
  return p_dict

def keep_targets(p_dict, keep_list):
  
  i = 0
  for planet in p_dict['pl_name']:
    if planet in keep_list:
      print 'KEEPING', planet
    else:
      print 'KILLING', planet
      for key in p_dict.keys():
	p_dict[key] = delete(p_dict[key],i)
      i -= 1
    i += 1
  return p_dict

def add_secondary(p_dict, target_list):
  
  for string in target_list:
    newname = string + ' secondary' 
    if newname not in p_dict['pl_name']:
      i = 0
      for planet in p_dict['pl_name']:
	if string == planet:
	  print 'adding secondary', planet
	  for key in p_dict.keys():
	    if key == 'pl_name':
	      p_dict[key] = append(p_dict[key],(p_dict[key][i] + ' secondary'))
	    elif key == 'pl_tranmid':
	      p_dict[key] = append(p_dict[key],(p_dict[key][i] + 0.5*p_dict['pl_orbper'][i]))
	    else:
	      p_dict[key] = append(p_dict[key],p_dict[key][i])
	i += 1
  return p_dict

def print_comparisons(p_dict):

  for i in range(0,len(p_dict['pl_name'])):
    print p_dict['pl_name'][argsort(p_dict['trans_SN'])][::-1][i]
    target = array(p_dict['nearby_objects'])[argsort(p_dict['trans_SN'])][::-1][i]
    try:
      for x in range (0,len(array(p_dict['nearby_objects'])[argsort(p_dict['trans_SN'])][::-1][i]['Vmag'])):
	targ_ra = target['_RAJ2000'][x]*(pi/180.0)
	targ_dec = target['_DEJ2000'][x]*(pi/180.0)      
	ra = p_dict['ra'][argsort(p_dict['trans_SN'])][::-1][i]*(pi/180.0)
	dec = p_dict['dec'][argsort(p_dict['trans_SN'])][::-1][i]*(pi/180.0)

	sep = 60*(arccos(sin(targ_dec)*sin(dec) + cos(dec)*cos(targ_dec)*cos(ra - targ_ra)))*180.0/pi

	print sep, target['Bmag'][x], target['Vmag'][x], target['B-V'][x]
    except:
	print 'NO MATCHES!'

def print_comparisons(p_dict):

  for i in range(0,len(p_dict['pl_name'])):
    print p_dict['pl_name'][argsort(p_dict['trans_SN'])][::-1][i]
    target = array(p_dict['nearby_objects'])[argsort(p_dict['trans_SN'])][::-1][i]
    try:
      for x in range (0,len(array(p_dict['nearby_objects'])[argsort(p_dict['trans_SN'])][::-1][i]['Vmag'])):
	targ_ra = target['_RAJ2000'][x]*(pi/180.0)
	targ_dec = target['_DEJ2000'][x]*(pi/180.0)      
	ra = p_dict['ra'][argsort(p_dict['trans_SN'])][::-1][i]*(pi/180.0)
	dec = p_dict['dec'][argsort(p_dict['trans_SN'])][::-1][i]*(pi/180.0)

	sep = 60*(arccos(sin(targ_dec)*sin(dec) + cos(dec)*cos(targ_dec)*cos(ra - targ_ra)))*180.0/pi

	print sep, target['Bmag'][x], target['Vmag'][x], target['B-V'][x]
    except:
	print 'NO MATCHES!'


def print_result(p_dict):

  print ' '
  print 'Planet name | Transit signal | Reflection signal | est S/N | est reflection S/N'

  sort_key = 'signal'

  sort_index = argsort(p_dict[sort_key])[::-1]

  for i in range(0,len(p_dict['pl_name'])):
    print p_dict['pl_name'][sort_index[i]], round(p_dict['signal'][sort_index[i]]*1e4,2), 'x 1e-4', round(p_dict['reflect_signal'][sort_index[i]]*1e4,2), 'x 1e-4',p_dict['trans_SN'][sort_index][i],p_dict['reflect_SN'][sort_index][i],p_dict['st_vj'][sort_index][i]

def print_positions(p_dict):

  for i in range(0,len(p_dict['pl_name'])):
    print '['+p_dict['pl_name'][argsort(p_dict['trans_SN'])][::-1][i].replace(' ','')+']',p_dict['ra'][argsort(p_dict['trans_SN'])][::-1][i], p_dict['dec'][argsort(p_dict['trans_SN'])][::-1][i], ';'

def grab_table(Columns,Table='exoplanets',Order ='dec',Format = 'csv',Transit=False):

  Base = "http://exoplanetarchive.ipac.caltech.edu/cgi-bin/nstedAPI/nph-nstedAPI?"

  if Transit == True:
    Columns += ',pl_tranflag'

  Query = Base+'table='+Table+'&select='+Columns+'&order='+Order+'&format='+Format

  response = urllib2.urlopen(Query)

  cr = csv.reader(response)

  name_array = Columns.split(',')

  p_dict = {}

  for i in range(0,len(Columns.split(','))):
    p_dict[name_array[i]] = np.array([])

  l = 0

  for line in cr:
    if l != 0:
      for i in range(0,len(Columns.split(','))):
	val = line[i]
	try:
	  val = float(val)
	except:
	  val = val
	if val == '':
	  val = 0.0
	p_dict[name_array[i]] = np.append(p_dict[name_array[i]],[val])
    l +=1

  if Transit==True:
    is_transit = [p_dict['pl_tranflag'] == 1]
    for i in range(0,len(Columns.split(','))):
      p_dict[name_array[i]] = p_dict[name_array[i]][is_transit]
    p_dict.pop('pl_tranflag', None)

  return p_dict

def calc_eqtemp(p_dict):

  sol_rad = 695500.0e3
  AU = 149597870700.0

  p_dict['eqt'] = p_dict['pl_eqt'].copy()

  for i in range(0,len(p_dict['pl_name'])):
    eqt = p_dict['st_teff'][i]*sqrt((sol_rad*p_dict['st_rad'][i])/(2.0*AU*p_dict['pl_orbsmax'][i]))
    if eqt != inf:
      p_dict['eqt'][i] = eqt

  return p_dict

def try_calc_seps(p_dict,i,search_radius,v):
  try:
    result = v.query_object(p_dict['pl_hostname'][i],catalog=['UCAC'])
    result = v.query_region(coord.ICRS(ra=p_dict['ra'][i],dec=p_dict['dec'][i], unit=(u.deg,u.deg)),radius=search_radius,catalog=['UCAC'])
    viz_result = result['I/322A/out']
    out_dict = {}
    for key in viz_result.keys():
      out_dict[key] = array(viz_result[key].filled(999).tolist())
#      print out_dict['Bmag'], out_dict['Vmag']
    out_dict['B-V'] = out_dict['Bmag'] - out_dict['Vmag']

    search_ra = p_dict['ra'][i]*(pi/180.0)
    search_dec = p_dict['dec'][i]*(pi/180.0)

    sep_list = []
    for x in range(0,len(out_dict['RAJ2000'])):
      targ_ra = out_dict['RAJ2000'][x]*(pi/180.0)
      targ_dec = out_dict['DEJ2000'][x]*(pi/180.0)      
      sep = 60*(arccos(sin(targ_dec)*sin(search_dec) + cos(search_dec)*cos(targ_dec)*cos(search_ra - targ_ra)))*180.0/pi
      sep_list += [sep]
    out_dict['sep'] = sep_list

    p_dict['nearby_objects'][i] = out_dict
  except:
    p_dict['nearby_objects'][i] = 'NO MATCHES!'
  return p_dict 

def find_IR_mags(p_dict,search_radius=0.05,verbose=False):

  search_radius = str(search_radius)+"m"
  maxjmag = 19
  maxj = '<'+str(maxjmag)

  nplanets = len(p_dict['pl_hostname'])

  p_dict['st_jj'] = (p_dict['st_vj'].copy())*0.0
  p_dict['st_hj'] = (p_dict['st_vj'].copy())*0.0
  p_dict['st_kj'] = (p_dict['st_vj'].copy())*0.0

  for i in range(0,nplanets):
    v = Vizier(columns=['+_r','Jmag','Hmag','Kmag'],column_filters={"Jmag":maxj})
    result = v.query_region(coord.ICRS(ra=p_dict['ra'][i],dec=p_dict['dec'][i], unit=(u.deg,u.deg)),radius=search_radius,catalog=['2MASS'])
    viz_result = result['II/246/out']

    jmag = viz_result['Jmag'][np.argsort(viz_result['_r'])[0]]
    hmag = viz_result['Hmag'][np.argsort(viz_result['_r'])[0]]
    kmag = viz_result['Kmag'][np.argsort(viz_result['_r'])[0]]

    vmag = p_dict['st_vj'][i]

    p_dict['st_jj'][i] = jmag
    p_dict['st_hj'][i] = hmag
    p_dict['st_kj'][i] = kmag

    if verbose == True:
      print p_dict['pl_name'][i], p_dict['ra'][i], p_dict['dec'][i]
      print 100.0*i/len(p_dict['pl_hostname']),'%'
      print 'best match', np.sort(viz_result['_r'])[0]
      print 'J MAG:',jmag,' (vmag ',vmag,')'
      print ''

  return p_dict

def vizier_find_mags(p_dict,search_radius):

  max_delta_m = 1

  search_radius = str(search_radius)+"m"

  p_dict['nearby_objects'] = [[]]*len((p_dict['pl_hostname']))

  for i in range(0,len(p_dict['pl_hostname'])):
    if p_dict['st_vj'][i] == 0:
      maxv = '<'+str(19)
    else:
      maxv = '<'+str(p_dict['st_vj'][i] + max_delta_m)
    v = Vizier(columns=['+_r','RAJ2000', 'DEJ2000','B-V', 'Vmag','Bmag'],column_filters={"Vmag":maxv})

    print p_dict['pl_name'][i], p_dict['ra'][i], p_dict['dec'][i]
    print 100.0*i/len(p_dict['pl_hostname']),'%'

    p_dict = try_calc_seps(p_dict,i,search_radius,v)
    print 'best match', np.sort(p_dict['nearby_objects'][i]['sep'])[0]

#    print p_dict['nearby_objects'][i]
#    if (p_dict['st_vj'][i] == 0) & (np.sort(p_dict['nearby_objects'][i]['sep'])[0] < 1):
    newmag = p_dict['nearby_objects'][i]['Vmag'][np.argsort(p_dict['nearby_objects'][i]['sep'])[0]]
    oldmag = p_dict['st_vj'][i]
    p_dict['st_vj'][i] = newmag
    print 'NEW MAG:',newmag,' (was ',oldmag,')'
    maxv = '<'+str(p_dict['st_vj'][i] + max_delta_m)
    v = Vizier(columns=['+_r','RAJ2000', 'DEJ2000','B-V', 'Vmag','Bmag'],column_filters={"Vmag":maxv})
    p_dict = try_calc_seps(p_dict,i,search_radius,v)

    print ''
  return p_dict

def simbad_find_mags(p_dict):


  Simbad.get_votable_fields()
  Simbad.add_votable_fields('Vmag')

  for i in range(0,len(p_dict['pl_hostname'])):
    result = Simbad.query_object(p_dict['pl_hostname'][i])
    print p_dict['pl_hostname'][i]
    try:
      print result['I/322A/out']
    except:
      'No result found!'

def calc_signal(p_dict):

  mol_mass = 3.65e-27

  Kb = 1.3806488e-23

  surf_grav_jup = 24.79

  jup_rad = 69911e3

  sol_rad = 695500e3

  AU = 149597870700.0

  albedo = 0.3

  print p_dict['pl_trandur']

  mirror = pi*(420.0/2.0)**2.0

  p_dict["p_grav"] = surf_grav_jup*p_dict["pl_massj"]/(p_dict["pl_radj"])**2
  print p_dict["eqt"], p_dict["p_grav"]
  p_dict["H"] = Kb*p_dict["eqt"]/(p_dict["p_grav"]*mol_mass)
  p_dict["signal"] = 2.0*p_dict['H']*p_dict['pl_radj']*jup_rad/((p_dict['st_rad']*sol_rad)**2)

  p_dict["reflect_signal"] = albedo*((p_dict['pl_radj']*jup_rad)/(p_dict['pl_orbsmax']*AU))**2.0

  transmission = 1.0

  width = 1000.0

  len_transit = p_dict['pl_trandur']*24*3600.0

  flux = width*1000*10**(-1.0*(p_dict['st_vj'] - 0.03)/2.512)

  nphotons = len_transit*flux*mirror*transmission

  noise = (1.0/sqrt(nphotons))*sqrt(2)

  reflect_SN = p_dict["reflect_signal"] / noise
  trans_SN = p_dict["signal"] / noise

  p_dict['reflect_SN'] = reflect_SN
  p_dict['trans_SN'] = trans_SN

  return p_dict

def create_dict(fname,declim,maglim,search_radius):

  Columns= 'pl_name,pl_hostname,ra,dec,pl_massj,pl_radj,pl_eqt,st_rad,st_vj,st_teff,pl_orbsmax,pl_tranmid,pl_tranmiderr1,pl_tranmiderr2,pl_trandur,pl_trandurerr1,pl_trandurerr2,pl_orbper,pl_orbpererr1,pl_orbpererr2'
  p_dict = grab_table(Columns,Transit=True)
  p_dict = calc_eqtemp(p_dict)
  print len(p_dict['st_vj'])

  if declim == None:
    is_sane = [(p_dict['st_vj'] < maglim) & (p_dict['st_vj'] > -1)]
  else:
    is_sane = [(p_dict['dec'] < declim) & (p_dict['st_vj'] < maglim) & (p_dict['st_vj'] > -1)]

  print len(p_dict['st_vj'][is_sane])

  for i in range(0,len(p_dict.keys())):
    p_dict[p_dict.keys()[i]] = p_dict[p_dict.keys()[i]][is_sane]
  p_dict = vizier_find_mags(p_dict,search_radius)
  p_dict = find_IR_mags(p_dict)
  p_dict = calc_signal(p_dict)
  
  pickle.dump(p_dict, open( fname+".p", "wb" ) )