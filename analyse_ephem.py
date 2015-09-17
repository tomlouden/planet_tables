#!/usr/bin/python
# -*- coding: utf-8 -*-
from pylab import *
from astropy.time import Time
import re

def main():

  min_len = 3

  name = 'planets'

  start_dict=load_data(name+'_start_out')
  mid_dict=load_data(name+'_mid_out')
  end_dict=load_data(name+'_end_out')

  mid_dict = remove_bad_nights(start_dict,mid_dict,end_dict)

  chunks(mid_dict,min_len)

def remove_bad_nights(start_dict,mid_dict,end_dict):

  mask = []

  chosen = 'WASP-39 b'

  for i in range(0,len(mid_dict['name'])):
    mid_transit_times = mid_dict['jd'][mid_dict['name'][i] == mid_dict['name']]
    end_transit_times = end_dict['jd'][mid_dict['name'][i] == end_dict['name']]
    start_transit_times = start_dict['jd'][mid_dict['name'][i] == start_dict['name']]
    this_date = mid_dict['jd'][i]
    dif2 = end_transit_times - this_date
    dif = start_transit_times - this_date
    if any(abs(dif2) < 0.5) & any(abs(dif) < 0.5):
#      if mid_dict['name'][i] == chosen:
#	print mid_dict['date'][i], mid_dict['name'][i]
      mask += [True]
    else:
      mask += [False]
  mask = array(mask)
  for key in mid_dict.keys():
    mid_dict[key] = mid_dict[key][mask]
  return mid_dict


def chunks(e_dict,min_len):
  length_chunk = 0
  dates = array([])
  names = array([])
  day = array([])

  for i in range(1,len(e_dict['name'])):
    delta = e_dict['jd'][i] - e_dict['jd'][i-1]
    if abs(delta) > 1.5:
      try:
	stretch = (max(dates) - min(dates))
	if len(dates) > min_len:
	  print 'A succesful run!'
	  for x in range(0,len(dates)):
	    print dates[x], names[x], day[x]
	  print ''
	dates=[]
	names=[]
	day=[]
      except:
	dates=[]
	names=[]
	day=[]

    else:
      dates = append(dates,e_dict['jd'][i])
      names = append(names,e_dict['name'][i])
      day = append(day,e_dict['date'][i])

def load_data(fname):
  e_dict = {}
  e_dict['name']= array([])
  e_dict['date']= array([])
  e_dict['time']= array([])
  e_dict['airmass']= array([])
  e_dict['jd']= array([])

  for line in open(fname,'rb'):
    data = line.strip('\n')
    print data
    data = re.sub('  +', ',', data)
    print data
    data = data.split(',')
    print data
    e_dict['name'] = append(e_dict['name'],data[0])
    e_dict['date'] = append(e_dict['date'],data[1])
    e_dict['time'] = append(e_dict['time'],data[2])
    e_dict['airmass'] = append(e_dict['airmass'],float(data[-2]))
    split_date = data[1].lstrip(' ').rstrip(' ').split(' ')
    iso_date = split_date[2] + '-' + mon_replace(split_date[1]) + '-' + split_date[0]
    iso_time = iso_date + data[2]
    t = Time(iso_time, format='iso', scale='utc')
    e_dict['jd'] = append(e_dict['jd'],t.jd)

  return e_dict

def mon_replace(mon):
  if mon == 'Jan':
    return '01'
  if mon == 'Feb':
    return '02'
  if mon == 'Mar':
    return '03'
  if mon == 'Apr':
    return '04'
  if mon == 'May':
    return '05'
  if mon == 'Jun':
    return '06'
  if mon == 'Jul':
    return '07'
  if mon == 'Aug':
    return '08'
  if mon == 'Sep':
    return '09'
  if mon == 'Oct':
    return '10'
  if mon == 'Nov':
    return '11'
  if mon == 'Dec':
    return '12'


main()
