"""
Module containing utility functions

"""


import numpy as np
import os
from netCDF4 import num2date

class ShapeError(Exception):
    pass



def get_ncdates(nc, tvar='time'):
    """ Return dates from nercdf time coordinate """
    t = nc.variables[tvar]
    dts = num2date(t[:], t.units, calendar=t.calendar)
    
    return dts


def get_datestr(dates, date_format):
    """ Return string of form 'mindate-maxdate' for specified format """
    datestr = '%s-%s' % (dates.min().strftime(date_format), 
                        dates.max().strftime(date_format))
    
    return datestr


def get_savename(outdir, name, dates, date_format, suffix=''):
    """ Return savename for output file """
    datestr = get_datestr(dates, date_format)
    savename = os.path.join(outdir,'%s_%s%s' % (name, datestr, suffix))
    
    return savename


def get_daterange(dates1, dates2):
    """ Return min and max date range for overlapping period """
    mindt = max([dates1.min(),dates2.min()])
    maxdt = min([dates1.max(),dates2.max()])
    
    return mindt, maxdt


def get_dateind(dates, mindt, maxdt):
    """ Return index to extract data over specified date range """
    return (dates >= mindt) & (dates <= maxdt)


def get_indrange(vals,minval,maxval):
    """ 
    Return max and min indices to use for simple slicing when
    updating multi-dim np.array that avoids creation of copies.
        
    """ 
    
    if vals.ndim != 1:
        raise ShapeError('get_inds: expected 1-d numpy array')
        
    inds = np.where((vals >= minval) & (vals < maxval))[0]
    
    if len(inds) != 0:
        minind = inds.min()
        maxind = inds.max() + 1
    else:
        raise ValueError('get_inds: no matching data')
        
    return minind, maxind
    
    

def find_nearest(array, val, min=False):
    """ Returns index for value in array that is closest to val. """
    if not isinstance(array, np.ndarray):
        raise TypeError('Array must be a np.ndarray object')

    abs_array = np.absolute(array - val)
    minval = abs_array.min()
    inds = np.where(abs_array == minval)[0]

    if len(inds) > 1:
        if min:
            inds = np.where(array == array[inds].min())[0]
        else:
            inds = np.where(array == array[inds].max())[0]
    
    return inds
