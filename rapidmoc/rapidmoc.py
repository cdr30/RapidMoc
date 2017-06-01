"""
Module containing main routines to execute RapidMoc

"""

import argparse
import ConfigParser
import matplotlib.pyplot as plt
import copy

import utils    
import sections
import transports

def get_args():
    """   Get arguments from command line.  """
    parser = argparse.ArgumentParser(
        description='Calculate RAPID AMOC diagnostics using ocean model data')
    parser.add_argument(
        'config_file', type=str, help='Path to configuration file.')
    parser.add_argument(
        'tfile', type=str, help='Path for netcdf file(s) containing temperature data.')
    parser.add_argument(
        'sfile', type=str, help='Path for netcdf file(s) containing salinity data.')
    parser.add_argument(
        'ufile', type=str, help='Path for netcdf file(s) containing zonal wind stress data.')
    parser.add_argument(
        'vfile', type=str, help='Path for netcdf file(s) containing meridional velocity data.')
    args = parser.parse_args()

    return args

 
def get_config(args):
    """ Return configuration options as <ConfigParser> object. """
    config = ConfigParser.ConfigParser()
    config.read(args.config_file)

    return config
    

def main():
    """ Parse options and run RapidMoc. """
    args = get_args()
    config = get_config(args)

    # Read data
    t = sections.ZonalSections(args.tfile, config, 'temperature')
    s = sections.ZonalSections(args.sfile, config, 'salinity')
    tau = sections.ZonalSections(args.ufile, config, 'taux')
    v = sections.ZonalSections(args.vfile, config, 'meridional_velocity')

    # Interpolate T & S data onto v-grid
    t_on_v = sections.interpolate(t, v)
    s_on_v = sections.interpolate(s, v)
  
    # Calculate dynamic heights
    dh = transports.dynamic_heights(t_on_v, s_on_v)
    
    # Calculate geostrophic transports
    if config.getboolean('rapid_options', 'geo_approx'):
        georef = config.getfloat('rapid_options', 'georef_level')
        vgeo = transports.geostrophic(v, dh, georef=georef)
        
        # Reference geostrophic transports to model velocities
        if config.has_option('rapid_options', 'vref_level'):
            vref_level = config.getfloat('rapid_options', 'vref_level') 
            vgeo = transports.update_reference(vgeo, v, vref_level)
    else:
        vgeo = copy.deepcopy(v)

    # Calculate Ekman velocities
    ek = transports.ekman(tau, v, config)

    # Use model velocities in FS and WBW regions
    vgeo = transports.merge(vgeo, v, config)

    # Apply mass-balance constraints to section
    vgeo = transports.rapid_mass_balance(vgeo, ek, config)

    # Add ekman to geostrophic transports for combined rapid velocities
    vrapid = copy.deepcopy(vgeo)
    vrapid.data = vgeo.data + ek.data
    
    # Get sub-section indices
    fs_minlon = config.getfloat('rapid_options','fs_minlon')
    fs_maxlon = config.getfloat('rapid_options','fs_maxlon')
    wbw_maxlon = config.getfloat('rapid_options','wbw_maxlon')
    int_maxlon = config.getfloat('rapid_options','int_maxlon')
    fsmin, fsmax = utils.get_indrange(vgeo.x, fs_minlon, fs_maxlon)
    wbwmin, wbwmax = utils.get_indrange(vgeo.x, fs_maxlon, wbw_maxlon)
    intmin, intmax = utils.get_indrange(vgeo.x, wbw_maxlon, int_maxlon)

    # Calculate volume and heat transports on each (sub-)section
    fs_trans = transports.Transports(vgeo, t_on_v, fsmin, fsmax)
    wbw_trans = transports.Transports(vgeo, t_on_v, wbwmin, wbwmax)
    int_trans = transports.Transports(vgeo, t_on_v, intmin, intmax)
    ek_trans = transports.Transports(ek, t_on_v, intmin, intmax)
    model_trans = transports.Transports(v, t_on_v, fsmin, intmax)
    rapid_trans = transports.Transports(vrapid, t_on_v, fsmin, intmax)

    # !!! Check CMIP6 integrations !!!

    # Check OHT data
    # Save/write data
    # Optional plot data
 
