#!/usr/bin/env python2.7 

""" 
DESCRIPTION:
  Convert text data files containing Florida current 
  observations into netcdf data files

AUTHOR:

  Chris Roberts     ECMWF     2017

"""

import argparse
from netCDF4 import Dataset, date2num
import numpy as np
import datetime


def get_args():
    """   Get arguments from command line.  """
    parser = argparse.ArgumentParser(
        description='Convert Florida current text data to netcdf')
    parser.add_argument(
        'dataf', type=str, help='Name of data file')
    args = parser.parse_args()

    return args


def read_txt(dataf):
    """ Read Florida Current data from text file """
    dts, dat  = [], []
    
    with open(dataf) as txtf:
        lines = txtf.readlines()
        
    lines = [l.split() for l in lines]
        
    for l in lines:
        if l[0] != '%':
            yy, mm, dd = np.int(l[0]), np.int(l[1]), np.int(l[2])
            dts.append(datetime.datetime(yy,mm,dd))
            
            if l[3].lower() == 'nan':
                dat.append(1.e20)
            else:
                dat.append(l[3])
                
    return np.array(dts), np.array(dat)


def create_netcdf(savef, dts, dat):
    """ Write Florida current data to netcdf file """
    dataset = Dataset(savef, 'w', format='NETCDF4_CLASSIC') 

    # Create time coordinate
    tdim = dataset.createDimension('time', None)
    time = dataset.createVariable('time',np.float64,(tdim.name,))
    time.units = 'hours since 0001-01-01 00:00:00.0'
    time.calendar = 'gregorian'
    time[:] = date2num(dts, time.units, calendar=time.calendar)

    # Create data variable
    fc = dataset.createVariable('florida_current_transport',np.float64,(tdim.name),fill_value=1.e20)
    fc.units = 'Sv'
    fc[:] = dat
    
    # close file
    print 'SAVING: %s' % savef 
    dataset.close()
    

def dat2nc(f):
    """ Convert Florida current text data to netcdf """ 
    savef = '%s.nc' % f
    dts, dat = read_txt(f)
    create_netcdf(savef, dts, dat)

    
if __name__ == '__main__':
    args = get_args()
    dat2nc(args.dataf)
    
