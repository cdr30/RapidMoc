"""
Module containing utility functions

"""


import numpy as np


class ShapeError(Exception):
    pass

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
