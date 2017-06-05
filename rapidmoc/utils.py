"""
Module containing utility functions

"""


import numpy as np
import os

class ShapeError(Exception):
    pass


def get_savename(config, dates, suffix):
    """ Return savename for output file """
    mindt = dates.min().strftime('%Y%m%d')
    maxdt = dates.max().strftime('%Y%m%d')
    outdir = config.get('options', 'outdir')
    name = config.get('options', 'name')
    savename = os.path.join(outdir,'%s_%s-%s%s' % (name, mindt, maxdt, suffix))
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
    """
    Returs index for value in array that is closest to val.

    Args:

    * array:
        np.ndarray containing data

    * val:
        Value to search for in array.

    Kwargs:

    * min - boolean:
        If True and multiple values are equally close,
        return minimum value (def=False).

    Returns:
        Index for array.

    """
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
