"""
Module containing code to plot volume and heat transport diagnostics

"""

import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt

import observations
import utils


# COLORS
c1='#d95f02'
c2='#1b9e77'
c3='#7570b3'


def plot_streamfunctions(config, rapid_trans, model_trans, obs_moc):
    """ Plot time mean overturning stream functions """

    fig = plt.figure(figsize=(6,8))
        
    # Plot observations
    plt.plot(obs_moc.streamfunction.mean(axis=0), -obs_moc.z, '-k', linewidth=4,
             label='RAPID observations')
    for nt in range(len(obs_moc.dates)):
        plt.plot(obs_moc.streamfunction[nt,:], -obs_moc.z, '-k', linewidth=0.5,
                 alpha=0.4)
   
   # Plot model data
    plt.plot(model_trans.streamfunction.mean(axis=0) / 1e6,
             -model_trans.z, '-', color=c1, linewidth=4, 
             label=model_trans.name)
    for nt in range(len(model_trans.dates)):
        plt.plot(model_trans.streamfunction[nt,:] / 1e6 ,
                 -model_trans.z, '-',
                 color=c1, linewidth=0.5, alpha=0.4)
    
    # Plot model data using geostrophic assumption
    plt.plot(rapid_trans.streamfunction.mean(axis=0) / 1e6,
             -rapid_trans.z, '-', linewidth=4, color=c2,
             label="%s (geostrophic approx)" % rapid_trans.name)
    for nt in range(len(rapid_trans.dates)):
        plt.plot(rapid_trans.streamfunction[nt,:] / 1e6,
                 -rapid_trans.z, '-',
                 color=c2, linewidth=0.5, alpha=0.4)
        
    # Annotate plot
    plt.title('Atlantic overturning streamfunction at 26N')
    plt.xlabel('Sverdrups')
    plt.ylabel('Depth (m)')
    plt.legend(loc='best', fontsize=10)   
    
    # Save plot
    plt.tight_layout()
    suffix = '_overturning_streamfunctions.png'
    savef = utils.get_savename(config, model_trans.dates, suffix)
    print 'SAVING: %s' % savef
    fig.savefig(savef, resolution=300)
    plt.close()


def plot_transport_ts(config, obs_vol, ek_trans, int_trans, 
                      wbw_trans, fs_trans, depth=1000):
    """ Plot volume transport time series """
    fig = plt.figure(figsize=(8,12))
    
    ### Plot observed time series
    fig.add_subplot(2,1,1)
    
    # Labels 
    fs_label = 'Florida straits (%6.1f Sv)' % (obs_vol.fs.mean())
    ek_label = 'Ekman transport (%6.1f Sv)' % (obs_vol.ekman.mean())
    umo_label = 'Upper-mid ocean (%4.1f Sv)' % (obs_vol.umo.mean())
    moc_label = 'MOC (%6.1f Sv)' % (obs_vol.moc.mean())
    
    # Add to axis
    plt.plot(obs_vol.dates, obs_vol.fs, linewidth=4, color=c1, label=fs_label)
    plt.plot(obs_vol.dates, obs_vol.ekman, linewidth=4, color=c2, label=ek_label)
    plt.plot(obs_vol.dates, obs_vol.umo, linewidth=4, color=c3, label=umo_label)
    plt.plot(obs_vol.dates, obs_vol.moc, linewidth=4, color='k', label=moc_label)
    
    # Annotate 
    plt.xlabel('Date')
    plt.ylim([-40,50])
    plt.ylabel('Sverdrups')
    plt.title('Overturning components at 26N in RAPID observations')
    plt.legend(loc=8, fontsize=8, ncol=2)
    
    ### Plot model time series
    fig.add_subplot(2,1,2)
    zind = utils.find_nearest(fs_trans.z, depth)
    dates = fs_trans.dates

    # Extract model data at specified depth
    fs = fs_trans.streamfunction[:,zind] / 1e6
    ek = ek_trans.streamfunction[:,zind] / 1e6
    interior = int_trans.streamfunction[:,zind] / 1e6
    wbw = wbw_trans.streamfunction[:,zind] / 1e6
    
    # Combine to get RAPID components
    umo = interior + wbw
    moc = fs + ek + umo
    
    # Model labels
    fs_label = 'Florida straits (%6.1f Sv)' % (fs.mean())
    ek_label = 'Ekman transport (%6.1f Sv)' % (ek.mean())
    umo_label = 'Upper-mid ocean (%4.1f Sv)' % (umo.mean())
    moc_label = 'MOC (%6.1f Sv)' % (moc.mean())
    
    # Add to axis
    plt.plot(dates, fs, linewidth=4, color=c1, label=fs_label)
    plt.plot(dates, ek, linewidth=4, color=c2, label=ek_label)
    plt.plot(dates, umo, linewidth=4, color=c3, label=umo_label)
    plt.plot(dates, moc, linewidth=4, color='k', label=moc_label)
    
    # Annotate 
    name = config.get('options', 'name')
    plt.xlabel('Date')
    plt.ylim([-40,50])
    plt.ylabel('Sverdrups')
    plt.title('Overturning components at 26N in %s' % name)
    plt.legend(loc=8, fontsize=8, ncol=2)
   
    # Save plot
    plt.tight_layout()
    suffix = '_volume_transport_time_series.png'
    savef = utils.get_savename(config, dates, suffix)
    print 'SAVING: %s' % savef
    fig.savefig(savef, resolution=300)
    plt.close()


def load_observations(config):
    """ Load observational data sets """
    if config.has_option('observations', 'time_avg'):
        time_avg = config.get('observations', 'time_avg')
    else:
        time_avg = None

    moc = observations.StreamFunctionObs(
        config.get('observations', 'streamfunctions'), time_avg=time_avg)
    vol = observations.VolumeTransportObs(
        config.get('observations', 'volume_transports'), time_avg=time_avg)
    oht = observations.HeatTransportObs(
        config.get('observations', 'heat_transports'), time_avg=time_avg)    

    return moc, vol, oht


def plot_diagnostics(config, model_trans, rapid_trans, fs_trans,
                     wbw_trans, int_trans, ek_trans):
    """ Plot volume and heat transport diagnostics against RAPID observations """
    
    obs_moc, obs_vol, obs_oht = load_observations(config)
    plot_streamfunctions(config, rapid_trans, model_trans, obs_moc)


    # Plot transport time series
    plot_transport_ts(config, obs_vol, ek_trans, int_trans, wbw_trans, fs_trans)
    
    
    # Plot heat transport time series
    
    # Plot heat transport decomposition
    
    # Plot heat transport vs volume transports 
    
    import pdb; pdb.set_trace()
    
    

