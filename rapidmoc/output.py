"""
Module containing code to generate netcdf object for output.
 
"""

 
from netCDF4 import Dataset, date2num
import numpy as np
import os

import utils

               
def open_ncfile(config, dates):
    """ Return handle for netcdf file """
    outdir = config.get('output', 'outdir')
    name = config.get('output', 'name')
    date_format = config.get('output', 'date_format')
    
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
        
    savef = utils.get_savename(
        outdir, name, dates, date_format,
        suffix='_natl_meridional_transports_at_26N.nc')
    dataset = Dataset(savef, 'w', format='NETCDF4_CLASSIC')
    
    return dataset

def create_netcdf(config, rapid_trans, model_trans, fc_trans,
                  wbw_trans, int_trans, ek_trans):
    """ 
    Return diagnostics for comparison with RAPID-MOCHA observations
    in netcdf object for plotting/output.
    
    """
    
    # Configuration options
    zind = rapid_trans.streamfunction.mean(axis=0).argmax()
    fc_minlon = config.getfloat('options','fc_minlon') 
    fc_maxlon = config.getfloat('options','fc_maxlon')  
    wbw_maxlon = config.getfloat('options','wbw_maxlon')
    int_maxlon = config.getfloat('options','int_maxlon')
    georef = config.getfloat('options','georef_level')
    ek_level = config.getfloat('options','ekman_depth')
    
    try:
        vref_level =  config.getfloat('options','vref_level')
    except:
        vref_level = 'None'
    
    # Create netcdf file and add dimensions
    dataset = open_ncfile(config, rapid_trans.dates)
    zdim = dataset.createDimension('depth', rapid_trans.z.size)
    tdim = dataset.createDimension('time', None)
    
    # Create time coordinate
    time = dataset.createVariable('time',np.float64,(tdim.name,))
    time.units = 'hours since 0001-01-01 00:00:00.0'
    time.calendar = 'gregorian'
    time[:] = date2num(rapid_trans.dates, time.units, calendar=time.calendar)

    # Create depth coordinate 
    z = dataset.createVariable('depth',np.float64,(zdim.name,))
    z.units = 'm'
    z[:] = rapid_trans.z
    
    # Create depth coordinate 
    dz = dataset.createVariable('level_thickness',np.float64,(zdim.name,))
    dz.units = 'm'
    dz[:] = rapid_trans.dz
    
    # Global attributes
    dataset.geostrophic_reference_level = georef
    dataset.rhocp = rapid_trans.rhocp
    dataset.ekman_level = ek_level
    dataset.reference_to_model_velocity = vref_level
    dataset.contact = 'chris.roberts@ecmwf.int'
    dataset.code_reference = 'Roberts, C. D., et al. (2013), Atmosphere drives recent interannual variability of the Atlantic meridional overturning circulation at 26.5N, Geophys. Res. Lett., 40, 5164-5170 doi:10.1002/grl.50930.'
    dataset.method_references = '(1) McCarthy, G. D., and Coauthors, 2015: Measuring the Atlantic Meridional Overturning Circulation at 26 degrees N. Progress in Oceanography, 130, 91-111. (2) Johns, W.E., M.O. Baringer, L.M. Beal, S.A. Cunningham, T. Kanzow, H.L. Bryden, J.J. Hirschi, J. Marotzke, C.S. Meinen, B. Shaw, and R. Curry, 2011: Continuous, Array-Based Estimates of Atlantic Ocean Heat Transport at 26.5N. J. Climate, 24, 2429-2449, doi: 10.1175/2010JCLI3997.1.'
    
    # Basinwide potential temperature profile
    t_basin = dataset.createVariable('t_basin',np.float64,(tdim.name,zdim.name))
    t_basin.units = 'degC'
    t_basin.minimum_longitude = fc_minlon
    t_basin.maximum_longitude = int_maxlon
    t_basin.comment = 'Basinwide zonal mean potential temperature profile'
    t_basin[:] = rapid_trans.zonal_avg_t

    # Florida current flow-weighted potential temperature
    t_fc_fwt = dataset.createVariable('t_fc_fwt',np.float64,(tdim.name))
    t_fc_fwt.units = 'degC'
    t_fc_fwt.minimum_longitude = fc_minlon
    t_fc_fwt.maximum_longitude = fc_maxlon
    t_fc_fwt.comment = 'Florida current flow-weighted potential temperature'
    t_fc_fwt[:] = fc_trans.oht_total / (fc_trans.rhocp * fc_trans.net_transport)

    # Basinwide salinity profile
    s_basin = dataset.createVariable('s_basin',np.float64,(tdim.name,zdim.name))
    s_basin.units = 'PSU'
    s_basin.minimum_longitude = fc_minlon
    s_basin.maximum_longitude = int_maxlon
    s_basin.comment = 'Basinwide zonal mean salinity profile'
    s_basin[:] = rapid_trans.zonal_avg_s

    # Florida current flow-weighted salinity
    s_fc_fwt = dataset.createVariable('s_fc_fwt',np.float64,(tdim.name))
    s_fc_fwt.units = 'PSU'
    s_fc_fwt.minimum_longitude = fc_minlon
    s_fc_fwt.maximum_longitude = fc_maxlon
    s_fc_fwt.comment = 'Florida current flow-weighted salinity'
    s_fc_fwt[:] = fc_trans.oft_total / ((-1.0/fc_trans.sref) * fc_trans.net_transport)
        
    # Basinwide transport profile - RAPID approx
    v_basin_rapid = dataset.createVariable('v_basin_rapid',np.float64,(tdim.name,zdim.name))
    v_basin_rapid.units = 'Sv/m'
    v_basin_rapid.minimum_longitude = fc_minlon
    v_basin_rapid.maximum_longitude = int_maxlon
    v_basin_rapid.comment = 'Basinwide transport profile using RAPID approximations'
    v_basin_rapid[:] = rapid_trans.zonal_sum_v / 1e6
    
    # Basinwide transport profile - model v
    v_basin_model = dataset.createVariable('v_basin_model',np.float64,(tdim.name,zdim.name))
    v_basin_model.units = 'Sv/m'
    v_basin_model.minimum_longitude = fc_minlon
    v_basin_model.maximum_longitude = int_maxlon
    v_basin_model.comment = 'Basinwide transport profile using model velocities'
    v_basin_model[:] = model_trans.zonal_sum_v / 1e6
    
    # Florida current transport profile
    v_fc = dataset.createVariable('v_fc',np.float64,(tdim.name,zdim.name))
    v_fc.units = 'Sv/m'
    v_fc.minimum_longitude = fc_minlon
    v_fc.maximum_longitude = fc_maxlon
    v_fc.comment = 'Florida current transport profile'
    v_fc[:] = fc_trans.zonal_sum_v / 1e6
    
    # Ekman transport time series
    ekman = dataset.createVariable('ekman',np.float64,(tdim.name))
    ekman.units = 'Sv'
    ekman.minimum_longitude = wbw_maxlon
    ekman.maximum_longitude = int_maxlon
    ekman.comment = 'Ekman transport time series (streamfunction at 1000m)'
    ekman[:] = ek_trans.streamfunction[:,zind] / 1e6
    
    # Gyre interior transport time series
    geoint = dataset.createVariable('geoint',np.float64,(tdim.name))
    geoint.units = 'Sv'
    geoint.minimum_longitude = wbw_maxlon
    geoint.maximum_longitude = int_maxlon
    geoint.comment = 'Geostrophic interior transport time series (streamfunction at 1000m).'
    geoint[:] = int_trans.streamfunction[:,zind] / 1e6

    # Western-boundary wedge transport time series
    wbw = dataset.createVariable('wbw',np.float64,(tdim.name))
    wbw.units = 'Sv'
    wbw.minimum_longitude = fc_maxlon
    wbw.maximum_longitude = wbw_maxlon
    wbw.comment = 'Western boundary wedge transport time series (streamfunction at 1000m).'
    wbw[:] = wbw_trans.streamfunction[:,zind] / 1e6

    # Florida current transport time series
    fc = dataset.createVariable('fc',np.float64,(tdim.name))
    fc.units = 'Sv'
    fc.minimum_longitude = fc_minlon
    fc.maximum_longitude = fc_maxlon
    fc.comment = 'Florida current transport time series (streamfunction at 1000m).'
    fc[:] = fc_trans.streamfunction[:,zind] / 1e6

    # Upper mid ocean transport time series
    umo = dataset.createVariable('umo',np.float64,(tdim.name))
    umo.units = 'Sv'
    umo.minimum_longitude = fc_maxlon
    umo.maximum_longitude = int_maxlon
    umo.comment = 'Upper mid-ocean transport time series (streamfunction at 1000m). umo = wbw + geoint'
    umo[:] = wbw[:] + geoint[:]
    
    # Meridional overturning transport time series - RAPID approx
    moc_rapid = dataset.createVariable('moc_rapid',np.float64,(tdim.name))
    moc_rapid.units = 'Sv'
    moc_rapid.minimum_longitude = fc_minlon
    moc_rapid.maximum_longitude = int_maxlon
    moc_rapid.comment = 'Time series of meridional overturning transport using RAPID approximation (streamfunction at 1000m)'
    moc_rapid[:] =  rapid_trans.streamfunction[:,zind] / 1e6

    # Meridional overturning transport time series - model v
    moc_model = dataset.createVariable('moc_model',np.float64,(tdim.name))
    moc_model.units = 'Sv'
    moc_model.minimum_longitude = fc_minlon
    moc_model.maximum_longitude = int_maxlon
    moc_model.comment = 'Time series of meridional overturning transport using model velocities (streamfunction at 1000m)'
    moc_model[:] =  model_trans.streamfunction[:,zind] / 1e6
    
    # Meridional overturning transport maxima time series - RAPID approx
    mocmax_rapid = dataset.createVariable('mocmax_rapid',np.float64,(tdim.name))
    mocmax_rapid.units = 'Sv'
    mocmax_rapid.minimum_longitude = fc_minlon
    mocmax_rapid.maximum_longitude = int_maxlon
    mocmax_rapid.comment = 'Time series of meridional overturning transport using RAPID approximation (streamfunction maxima)'
    mocmax_rapid[:] =  rapid_trans.streamfunction.max(axis=1) / 1e6

    # Meridional overturning transport maxima time series - model v
    mocmax_model = dataset.createVariable('mocmax_model',np.float64,(tdim.name))
    mocmax_model.units = 'Sv'
    mocmax_model.minimum_longitude = fc_minlon
    mocmax_model.maximum_longitude = int_maxlon
    mocmax_model.comment = 'Time series of meridional overturning transport using model velocities (streamfunction maxima)'
    mocmax_model[:] =  model_trans.streamfunction.max(axis=1) / 1e6
        
    # Overturning streamfunctions - RAPID approx
    sf_rapid = dataset.createVariable('sf_rapid',np.float64,(tdim.name,zdim.name))
    sf_rapid.units = 'Sv'
    sf_rapid.minimum_longitude = fc_minlon
    sf_rapid.maximum_longitude = int_maxlon
    sf_rapid.comment = 'Overturning streamfunctions using RAPID approximation.'
    sf_rapid[:] =  rapid_trans.streamfunction/ 1e6

    # Meridional overturning transport time series - model v
    sf_model = dataset.createVariable('sf_model',np.float64,(tdim.name,zdim.name))
    sf_model.units = 'Sv'
    sf_model.minimum_longitude = fc_minlon
    sf_model.maximum_longitude = int_maxlon
    sf_model.comment = 'Overturning streamfunctions using model velocities.'
    sf_model[:] =  model_trans.streamfunction / 1e6
    
    # Florida current stream function 
    sf_fc = dataset.createVariable('sf_fc',np.float64,(tdim.name,zdim.name))
    sf_fc.units = 'Sv'
    sf_fc.minimum_longitude = fc_minlon
    sf_fc.maximum_longitude = fc_maxlon
    sf_fc.comment = 'Florida current overturning streamfunction.'
    sf_fc[:] =  fc_trans.streamfunction/ 1e6
    
    # Ekman stream function 
    sf_ek = dataset.createVariable('sf_ek',np.float64,(tdim.name,zdim.name))
    sf_ek.units = 'Sv'
    sf_ek.minimum_longitude = wbw_maxlon
    sf_ek.maximum_longitude = int_maxlon
    sf_ek.comment = 'Ekman overturning streamfunction.'
    sf_ek[:] =  ek_trans.streamfunction/ 1e6
    
    # Wbw stream function 
    sf_wbw = dataset.createVariable('sf_wbw',np.float64,(tdim.name,zdim.name))
    sf_wbw.units = 'Sv'
    sf_wbw.minimum_longitude = fc_minlon
    sf_wbw.maximum_longitude = wbw_maxlon
    sf_wbw.comment = 'Western boundary wedge overturning streamfunction.'
    sf_wbw[:] =  wbw_trans.streamfunction/ 1e6
    
    # Geostrophic interior stream function 
    sf_geoint = dataset.createVariable('sf_geoint',np.float64,(tdim.name,zdim.name))
    sf_geoint.units = 'Sv'
    sf_geoint.minimum_longitude = wbw_maxlon
    sf_geoint.maximum_longitude = int_maxlon
    sf_geoint.comment = 'Geostrophic interior overturning streamfunction.'
    sf_geoint[:] =  int_trans.streamfunction/ 1e6
    
    # mid ocean stream function 
    sf_mo = dataset.createVariable('sf_mo',np.float64,(tdim.name,zdim.name))
    sf_mo.units = 'Sv'
    sf_mo.minimum_longitude = fc_maxlon
    sf_mo.maximum_longitude = int_maxlon
    sf_mo.comment = 'Mid ocean overturning streamfunction (sf_mo = sf_wbw + sf_int).'
    sf_mo[:] =  sf_geoint[:] + sf_wbw[:]
    
    # Total heat transport - RAPID approx
    q_sum_rapid = dataset.createVariable('q_sum_rapid',np.float64,(tdim.name))
    q_sum_rapid.units = 'PW'
    q_sum_rapid.minimum_longitude = fc_minlon
    q_sum_rapid.maximum_longitude = int_maxlon
    q_sum_rapid.comment = 'Total heat transport across section calculated using RAPID approximations (q_sum_rapid = q_fc + q_ek + q_mo = q_ot_rapid + q_gyre_rapid + q_net_rapid)'
    q_sum_rapid[:] = rapid_trans.oht_total / 1e15
        
    # Gyre heat transport - RAPID approx
    q_gyre_rapid = dataset.createVariable('q_gyre_rapid',np.float64,(tdim.name))
    q_gyre_rapid.units = 'PW'
    q_gyre_rapid.minimum_longitude = fc_minlon
    q_gyre_rapid.maximum_longitude = int_maxlon
    q_gyre_rapid.comment = 'Heat transport by the horizontal circulation calculated using RAPID approximations '
    q_gyre_rapid[:] = rapid_trans.oht_by_horizontal / 1e15
    
    # Overturning heat transport - RAPID approx
    q_ot_rapid = dataset.createVariable('q_ot_rapid',np.float64,(tdim.name))
    q_ot_rapid.units = 'PW'
    q_ot_rapid.minimum_longitude = fc_minlon
    q_ot_rapid.maximum_longitude = int_maxlon
    q_ot_rapid.comment = 'Heat transport by the overturning circulation calculated using RAPID approximations'
    q_ot_rapid[:] = rapid_trans.oht_by_overturning / 1e15
    
    # Heat transport by net throughflow - RAPID approx
    q_net_rapid = dataset.createVariable('q_net_rapid',np.float64,(tdim.name))
    q_net_rapid.units = 'PW'
    q_net_rapid.minimum_longitude = fc_minlon
    q_net_rapid.maximum_longitude = int_maxlon
    q_net_rapid.comment = 'Heat transport referenced to 0C by the net flow through the section using RAPID approximations'
    q_net_rapid[:] = rapid_trans.oht_by_net / 1e15
    
    # Total heat transport - model v
    q_sum_model = dataset.createVariable('q_sum_model',np.float64,(tdim.name))
    q_sum_model.units = 'PW'
    q_sum_model.minimum_longitude = fc_minlon
    q_sum_model.maximum_longitude = int_maxlon
    q_sum_model.comment = 'Total heat transport across section calculated using model velocities (q_sum_model = q_gyre_model + q_ot_model + q_net_model)'
    q_sum_model[:] = model_trans.oht_total / 1e15
    
    # Gyre heat transport -model v
    q_gyre_model = dataset.createVariable('q_gyre_model',np.float64,(tdim.name))
    q_gyre_model.units = 'PW'
    q_gyre_model.minimum_longitude = fc_minlon
    q_gyre_model.maximum_longitude = int_maxlon
    q_gyre_model.comment = 'Heat transport by the horizontal circulation calculated using model velocities'
    q_gyre_model[:] = model_trans.oht_by_horizontal / 1e15
    
    # Overturning heat transport - model v
    q_ot_model = dataset.createVariable('q_ot_model',np.float64,(tdim.name))
    q_ot_model.units = 'PW'
    q_ot_model.minimum_longitude = fc_minlon
    q_ot_model.maximum_longitude = int_maxlon
    q_ot_model.comment = 'Heat transport by the overturning circulation calculated using model velocities'
    q_ot_model[:] = model_trans.oht_by_overturning / 1e15
    
    # Heat transport by net throughflow - model v
    q_net_model = dataset.createVariable('q_net_model',np.float64,(tdim.name))
    q_net_model.units = 'PW'
    q_net_model.minimum_longitude = fc_minlon
    q_net_model.maximum_longitude = int_maxlon
    q_net_model.comment = 'Heat transport referenced to 0C by the net flow through the section using model velocities'
    q_net_model[:] = model_trans.oht_by_net / 1e15
        
    # Heat transport by florida current
    q_fc = dataset.createVariable('q_fc',np.float64,(tdim.name))
    q_fc.units = 'PW'
    q_fc.minimum_longitude = fc_minlon
    q_fc.maximum_longitude = fc_maxlon
    q_fc.comment = 'Heat transport referenced to 0C by the Florida current'
    q_fc[:] = fc_trans.oht_total / 1e15
    
    # Heat transport by ekman
    q_ek = dataset.createVariable('q_ek',np.float64,(tdim.name))
    q_ek.units = 'PW'
    q_ek.minimum_longitude = wbw_maxlon
    q_ek.maximum_longitude = int_maxlon
    q_ek.comment = 'Heat transport referenced to 0C by Ekman transport'
    q_ek[:] = ek_trans.oht_total / 1e15
    
    # Heat transport by wbw
    q_wbw = dataset.createVariable('q_wbw',np.float64,(tdim.name))
    q_wbw.units = 'PW'
    q_wbw.minimum_longitude = fc_maxlon
    q_wbw.maximum_longitude = wbw_maxlon
    q_wbw.comment = 'Heat transport referenced to 0C by western boundary wedge transport'
    q_wbw[:] = wbw_trans.oht_total / 1e15
    
    # Heat transport by zonal mean geostrophic interior
    q_geoint = dataset.createVariable('q_geoint',np.float64,(tdim.name))
    q_geoint.units = 'PW'
    q_geoint.minimum_longitude = wbw_maxlon
    q_geoint.maximum_longitude = int_maxlon
    q_geoint.comment = 'Heat transport referenced to 0C by zonal mean of geostrophic interior transport'
    q_geoint[:] = (int_trans.oht_total - int_trans.oht_by_horizontal )/ 1e15
    
    # Heat transport by standing "eddy" component of geostrophic interior
    q_eddy = dataset.createVariable('q_eddy',np.float64,(tdim.name))
    q_eddy.units = 'PW'
    q_eddy.minimum_longitude = wbw_maxlon
    q_eddy.maximum_longitude = int_maxlon
    q_eddy.comment = 'Heat transport referenced to 0C by standing eddy component of geostrophic interior transport'
    q_eddy[:] = (int_trans.oht_by_horizontal )/ 1e15
        
    # Heat transport by mid ocean
    q_mo = dataset.createVariable('q_mo',np.float64,(tdim.name))
    q_mo.units = 'PW'
    q_mo.minimum_longitude = wbw_maxlon
    q_mo.maximum_longitude = int_maxlon
    q_mo.comment = 'Heat transport referenced to 0C by mid-ocean transport (q_mo = q_geoint + q_wbw + q_eddy)'
    q_mo[:] = q_geoint[:] + q_wbw[:] + q_eddy[:]

    # Total freshwater transport - RAPID approx
    fw_sum_rapid = dataset.createVariable('fw_sum_rapid',np.float64,(tdim.name))
    fw_sum_rapid.units = 'Sv'
    fw_sum_rapid.minimum_longitude = fc_minlon
    fw_sum_rapid.maximum_longitude = int_maxlon
    fw_sum_rapid.reference_salinity = rapid_trans.sref
    fw_sum_rapid.comment = 'Total equivalent freshwater transport across section calculated using RAPID approximations (fw_sum_rapid = fw_fc + fw_ek + fw_mo = fw_ot_rapid + fw_gyre_rapid + fw_net_rapid)'
    fw_sum_rapid[:] = rapid_trans.oft_total /1.0e6
        
    # Gyre freshwater transport - RAPID approx
    fw_gyre_rapid = dataset.createVariable('fw_gyre_rapid',np.float64,(tdim.name))
    fw_gyre_rapid.units = 'Sv'
    fw_gyre_rapid.minimum_longitude = fc_minlon
    fw_gyre_rapid.maximum_longitude = int_maxlon
    fw_gyre_rapid.comment = 'freshwater transport by the horizontal circulation calculated using RAPID approximations '
    fw_gyre_rapid[:] = rapid_trans.oft_by_horizontal/1.0e6 
    
    # Overturning freshwater transport - RAPID approx
    fw_ot_rapid = dataset.createVariable('fw_ot_rapid',np.float64,(tdim.name))
    fw_ot_rapid.units = 'Sv'
    fw_ot_rapid.minimum_longitude = fc_minlon
    fw_ot_rapid.maximum_longitude = int_maxlon
    fw_ot_rapid.comment = 'freshwater transport by the overturning circulation calculated using RAPID approximations'
    fw_ot_rapid[:] = rapid_trans.oft_by_overturning /1.0e6
    
    # freshwater transport by net throughflow - RAPID approx
    fw_net_rapid = dataset.createVariable('fw_net_rapid',np.float64,(tdim.name))
    fw_net_rapid.units = 'Sv'
    fw_net_rapid.minimum_longitude = fc_minlon
    fw_net_rapid.maximum_longitude = int_maxlon
    fw_net_rapid.reference_salinity = rapid_trans.sref
    fw_net_rapid.comment = 'equivalent freshwater transport by the net flow through the section using RAPID approximations'
    fw_net_rapid[:] = rapid_trans.oft_by_net /1.0e6
    
    # Total freshwater transport - model v
    fw_sum_model = dataset.createVariable('fw_sum_model',np.float64,(tdim.name))
    fw_sum_model.units = 'Sv'
    fw_sum_model.minimum_longitude = fc_minlon
    fw_sum_model.maximum_longitude = int_maxlon
    fw_sum_model.reference_salinity = model_trans.sref
    fw_sum_model.comment = 'Total freshwater transport across section calculated using model velocities (fw_sum_model = fw_gyre_model + fw_ot_model + fw_net_model)'
    fw_sum_model[:] = model_trans.oft_total /1.0e6
    
    # Gyre freshwater transport -model v
    fw_gyre_model = dataset.createVariable('fw_gyre_model',np.float64,(tdim.name))
    fw_gyre_model.units = 'Sv'
    fw_gyre_model.minimum_longitude = fc_minlon
    fw_gyre_model.maximum_longitude = int_maxlon
    fw_gyre_model.comment = 'freshwater transport by the horizontal circulation calculated using model velocities'
    fw_gyre_model[:] = model_trans.oft_by_horizontal /1.0e6
    
    # Overturning freshwater transport - model v
    fw_ot_model = dataset.createVariable('fw_ot_model',np.float64,(tdim.name))
    fw_ot_model.units = 'Sv'
    fw_ot_model.minimum_longitude = fc_minlon
    fw_ot_model.maximum_longitude = int_maxlon
    fw_ot_model.comment = 'freshwater transport by the overturning circulation calculated using model velocities'
    fw_ot_model[:] = model_trans.oft_by_overturning /1.0e6
    
    # freshwater transport by net throughflow - model v
    fw_net_model = dataset.createVariable('fw_net_model',np.float64,(tdim.name))
    fw_net_model.units = 'Sv'
    fw_net_model.minimum_longitude = fc_minlon
    fw_net_model.maximum_longitude = int_maxlon
    fw_net_model.reference_salinity = model_trans.sref
    fw_net_model.comment = 'equivalent freshwater transport by the net flow through the section using model velocities'
    fw_net_model[:] = model_trans.oft_by_net /1.0e6
        
    # freshwater transport by florida current
    fw_fc = dataset.createVariable('fw_fc',np.float64,(tdim.name))
    fw_fc.units = 'Sv'
    fw_fc.minimum_longitude = fc_minlon
    fw_fc.maximum_longitude = fc_maxlon
    fw_fc.reference_salinity = fc_trans.sref
    fw_fc.comment = 'equivalent freshwater transport by the Florida current'
    fw_fc[:] = fc_trans.oft_total /1.0e6
    
    # freshwater transport by ekman
    fw_ek = dataset.createVariable('fw_ek',np.float64,(tdim.name))
    fw_ek.units = 'Sv'
    fw_ek.minimum_longitude = wbw_maxlon
    fw_ek.maximum_longitude = int_maxlon
    fw_ek.reference_salinity = ek_trans.sref
    fw_ek.comment = 'equivalent freshwater transport by Ekman transport'
    fw_ek[:] = ek_trans.oft_total/1.0e6
    
    # freshwater transport by wbw
    fw_wbw = dataset.createVariable('fw_wbw',np.float64,(tdim.name))
    fw_wbw.units = 'Sv'
    fw_wbw.minimum_longitude = fc_maxlon
    fw_wbw.maximum_longitude = wbw_maxlon
    fw_wbw.reference_salinity = wbw_trans.sref
    fw_wbw.comment = 'equivalent freshwater transport by western boundary wedge transport'
    fw_wbw[:] = wbw_trans.oft_total /1.0e6
    
    # freshwater transport by zonal mean geostrophic interior
    fw_geoint = dataset.createVariable('fw_geoint',np.float64,(tdim.name))
    fw_geoint.units = 'Sv'
    fw_geoint.minimum_longitude = wbw_maxlon
    fw_geoint.maximum_longitude = int_maxlon
    fw_geoint.reference_salinity = int_trans.sref
    fw_geoint.comment = 'equivalent freshwater transport by zonal mean of geostrophic interior transport'
    fw_geoint[:] = (int_trans.oft_total - int_trans.oft_by_horizontal )/1.0e6
    
    # freshwater transport by standing "eddy" component of geostrophic interior
    fw_eddy = dataset.createVariable('fw_eddy',np.float64,(tdim.name))
    fw_eddy.units = 'Sv'
    fw_eddy.minimum_longitude = wbw_maxlon
    fw_eddy.maximum_longitude = int_maxlon
    fw_eddy.reference_salinity = int_trans.sref
    fw_eddy.comment = 'equivalent freshwater transport by standing eddy component of geostrophic interior transport'
    fw_eddy[:] = (int_trans.oft_by_horizontal )/1.0e6
        
    # freshwater transport by mid ocean
    fw_mo = dataset.createVariable('fw_mo',np.float64,(tdim.name))
    fw_mo.units = 'Sv'
    fw_mo.minimum_longitude = wbw_maxlon
    fw_mo.maximum_longitude = int_maxlon
    fw_mo.reference_salinity = int_trans.sref
    fw_mo.comment = 'equivalent freshwater transport by mid-ocean transport (fw_mo = fw_geoint + fw_wbw + fw_eddy)'
    fw_mo[:] = fw_geoint[:] + fw_wbw[:] + fw_eddy[:]

    return dataset
    
            
