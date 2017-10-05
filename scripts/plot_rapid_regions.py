#!/usr/bin/env python2.7 

""" 
DESCRIPTION:
  Plot meridional currents at RAPID section and overly with specified 
  values for fc_minlon, fc_maxlon and wbw_maxlon.

AUTHOR:

  Chris Roberts     ECMWF     2017

"""

import argparse
import numpy as np
import ConfigParser
from netCDF4 import Dataset
import matplotlib.pyplot as plt

import rapidmoc.sections as sections
import rapidmoc.plotdiag as rapidplot


def get_args():
    """   Get arguments from command line.  """
    parser = argparse.ArgumentParser(
        description=' Plot meridional currents at RAPID section and region boundaries.')
    parser.add_argument(
        'configf', type=str, help='RapidMoc config file')
    parser.add_argument(
        'vfile', type=str, help='Netcdf file containing meridional velocity data')
    parser.add_argument(
        '--name', type=str, help='Name to use on plot', default='model data')
    args = parser.parse_args()

    return args


def main(configf, vfile, name):
    """ Plot meridional currents at RAPID section and region boundaries. """

    # Read config file
    config = ConfigParser.ConfigParser()
    config.read(configf)

    # Read boundaries of RAPID regions
    fc_minlon = config.getfloat('options', 'fc_minlon')
    fc_maxlon = config.getfloat('options', 'fc_maxlon')
    wbw_maxlon = config.getfloat('options', 'wbw_maxlon')

    # Load velocity data
    v = sections.ZonalSections(vfile, config, 'meridional_velocity')
    lons = v.x
    z_scaled = -np.sqrt(np.abs(v.z))
    
    if len(v.data.shape) == 3:
        vdat = v.data[0]
    else:
        vdat = v.data
        
    # Plot velocity data and overly boundaries
    fig = plt.figure(figsize=(10,6))
    zticks, zlabels = rapidplot.get_zscale_ticks()
    plt.pcolormesh(lons, z_scaled, vdat, vmin=vdat.min(), vmax=np.abs(vdat.min()),
                   cmap=plt.cm.RdBu_r)
    
    for lon in [fc_minlon, fc_maxlon, wbw_maxlon]:
        plt.plot([lon, lon], [z_scaled.min(), z_scaled.max()], '-k', linewidth=2)

    # Annotate figure
    plt.xlim([fc_minlon-2, wbw_maxlon+2])
    plt.xlabel('Longitude')
    plt.ylabel('Depth (m)')
    plt.yticks(zticks, zlabels)
    plt.title('Meridional velocity in %s\n Black lines indicate position of fc_minlon, fc_maxlon, wbw_maxlon' % 
              name)
    plt.colorbar()
    plt.show()
    
if __name__ == '__main__':
    args = get_args()
    main(args.configf, args.vfile, args.name)
    
