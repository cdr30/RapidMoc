#!/usr/bin/env python3


import argparse
import configparser
from netCDF4 import Dataset

import rapidmoc.observations as rapidobs
import rapidmoc.plotdiag as rapidplot


def get_args():
    """   Get arguments from command line.  """
    parser = argparse.ArgumentParser(
        description=' Plot RapidMoc time series data against observations.')
    parser.add_argument(
        'configf', help='RapidMoc config file containing paths to observational data')
    parser.add_argument(
        'datafile', type=str, help='Output netcdf file created by RapidMoc')
    parser.add_argument(
        '--name', type=str, help='Name to use in plots and output files', default='model')

    args = parser.parse_args()

    return args


def main(args):
    """ Plot meridional currents at RAPID section and region boundaries. """

    # Read config file
    config = configparser.ConfigParser()
    config.read(args.configf)
    obs_oht_f = config.get('observations', 'heat_transports')
    obs_vol_f = config.get('observations', 'volume_transports')
    obs_sf_f = config.get('observations', 'streamfunctions')
    obs_fc_f = config.get('observations', 'florida_current')
    time_avg = config.get('observations', 'time_avg')
    outdir = config.get('output', 'outdir')
   
    # Read data 
    obs_oht = rapidobs.HeatTransportObs(obs_oht_f, time_avg=time_avg)
    obs_fc = rapidobs.FloridaCurrentObs(obs_fc_f, time_avg=time_avg)
    obs_sf = rapidobs.StreamFunctionObs(obs_sf_f, time_avg=time_avg)
    obs_vol = rapidobs.VolumeTransportObs(obs_vol_f, time_avg=time_avg)
    trans = Dataset(args.datafile)

    # Plot data
    rapidplot.plot_diagnostics(trans, name=args.name, outdir=outdir, 
				      obs_sf=obs_sf, obs_oht=obs_oht,
				      obs_vol=obs_vol, obs_fc=obs_fc)

if __name__ == '__main__':
    args = get_args()
    main(args)
    
