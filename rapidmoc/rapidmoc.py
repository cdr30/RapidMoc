"""
Module containing main routines to execute RapidMoc

"""

import argparse
import ConfigParser
import matplotlib.pyplot as plt
import copy

import utils    
import sections
import plotdiag
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
  
    # Return integrated transports on RAPID section as netcdf object
    trans = transports.calc_transports_from_sections(
        config, v, tau, t_on_v, s_on_v)
    
    # Plot data
    
    # Write data





    # Plot diagnostics
#    if config.getboolean('options','plot'):
#        plotdiag.plot_diagnostics(config, model_trans, rapid_trans, fs_trans,
#                              wbw_trans, int_trans, ek_trans)
        
    
    # Write data to netcdf file
    import pdb; pdb.set_trace()
    

    import pdb; pdb.set_trace()

    # Check OHT data
    # Save/write data
    # Optional plot data
 
