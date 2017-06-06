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
import observations
import plotdiag

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
        'taufile', type=str, help='Path for netcdf file(s) containing zonal wind stress data.')
    parser.add_argument(
        'vfile', type=str, help='Path for netcdf file(s) containing meridional velocity data.')
    args = parser.parse_args()

    return args

 
def get_config(args):
    """ Return configuration options as <ConfigParser> object. """
    config = ConfigParser.ConfigParser()
    config.read(args.config_file)

    return config
    

def load_observations(config):
    """ Load observational data sets """
    if config.has_option('observations', 'time_avg'):
        time_avg = config.get('observations', 'time_avg')
    else:
        time_avg = None

    sf = observations.StreamFunctionObs(
        config.get('observations', 'streamfunctions'), time_avg=time_avg)
    vol = observations.VolumeTransportObs(
        config.get('observations', 'volume_transports'), time_avg=time_avg)
    oht = observations.HeatTransportObs(
        config.get('observations', 'heat_transports'), time_avg=time_avg)    

    return sf, vol, oht


def call_plotdiag(config, trans):
    """ Call plotting routines to compare against RAPID observations """ 
    obs_sf, obs_vol, obs_oht = load_observations(config)        
    outdir = config.get('output', 'outdir')
    date_format = config.get('output', 'date_format')
    name = config.get('output', 'name')
    plotdiag.plot_diagnostics(trans, obs_sf, obs_vol, obs_oht, name=name,
                              outdir=outdir, date_format=date_format)


def main():
    """ Parse options and run RapidMoc. """
    args = get_args()
    config = get_config(args)

    # Read data
    t = sections.ZonalSections(args.tfile, config, 'temperature')
    s = sections.ZonalSections(args.sfile, config, 'salinity')
    tau = sections.ZonalSections(args.taufile, config, 'taux')
    v = sections.ZonalSections(args.vfile, config, 'meridional_velocity')

    # Interpolate T & S data onto v-grid
    t_on_v = sections.interpolate(t, v)
    s_on_v = sections.interpolate(s, v)
  
    # Return integrated transports on RAPID section as netcdf object
    trans = transports.calc_transports_from_sections(
        config, v, tau, t_on_v, s_on_v)
    
    ## Plot diagnostics
    if config.getboolean('output','plot'):
        call_plotdiag(config, trans)
        
    # Write data
    print 'SAVING: %s' % trans.filepath()
    trans.close()

