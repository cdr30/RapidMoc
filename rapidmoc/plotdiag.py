"""
Module containing code to plot volume and heat transport diagnostics

"""


import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm

import numpy as np
from scipy import stats

import observations
import utils


# COLORS 
c1='#a6cee3'
c2='#1f78b4'
c3='#b2df8a'
c4='#33a02c'


def plot_streamfunctions(trans, name='simulated', basename='', obs=None,lw=4):
    """ Plot time mean overturning stream functions"""

    # Extract variables from data objects
    z = trans.variables['depth'][:]
    sf_rapid = trans.variables['sf_rapid'][:].mean(axis=0)
    sf_model = trans.variables['sf_model'][:].mean(axis=0)
    sfmax_rapid = sf_rapid.max()
    zmax_rapid = z[np.argmax(sf_rapid)]
    sfmax_model = sf_model.max()
    zmax_model = z[np.argmax(sf_model)]

    # Create labels
    model_label = ('%s (model velocities) (max=%4.1f Sv, depth=%6i m)' %
                 (name, sfmax_model, zmax_model))
    rapid_label = ('%s (RAPID approx) (max=%4.1f Sv, depth=%6i m)' %
                 (name, sfmax_rapid, zmax_rapid))

    # Add data to axis
    fig = plt.figure(figsize=(6,8))
    plt.plot(sf_model, -z,'-', color=c1, linewidth=lw, label=model_label)
    plt.plot(sf_rapid, -z,'-', linewidth=lw, color=c2, label=rapid_label)

    # Plot optional observational data
    if obs is not None:
        z_obs = obs.z
        sf_obs = obs.sf.mean(axis=0)
        sfmax_obs = sf_obs.max()
        zmax_obs = z_obs[np.argmax(sf_obs)]
        obs_label = ('RAPID observations (max=%4.1f Sv, depth=%6i m)' %
                     (sfmax_obs, zmax_obs))
        plt.plot(sf_obs, -z_obs, '-k', linewidth=lw, label=obs_label)
        
    # Annotate plot
    plt.title('Atlantic overturning streamfunction at 26N')
    plt.xlabel('Sverdrups')
    plt.ylabel('Depth (m)')
    plt.legend(loc='best', fontsize=8)   

    # Save plot
    plt.tight_layout()
    savef = basename + 'overturning_streamfunctions_at_26n.png'
    print 'SAVING: %s' % savef
    fig.savefig(savef, resolution=300)
    plt.close()
    

def plot_streamfunction_hovmollers(trans, name='simulated', basename='', obs=None):
    """ Plot overturning stream function hovmoller diagrams"""

    # Extract variables from data objects
    dts = utils.get_ncdates(trans)
    z = trans.variables['depth'][:]
    sf_rapid = trans.variables['sf_rapid'][:]
    sf_model = trans.variables['sf_model'][:]
    
   # Set up figure
    fig = plt.figure(figsize=(8,12))
    cmap=plt.cm.viridis
    levels = np.arange(15) * 2 - 4
    norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
    cmin,cmax=-5,30
    
    # Add model data to axis
    fig.add_subplot(3,1,1)
    plt.pcolormesh(dts, -z, sf_model.transpose(), 
                   vmin=cmin,vmax=cmax,cmap=cmap, norm=norm)
    plt.colorbar(orientation='vertical')
    plt.title('Overturning streamfunction at 26N in %s (model velocities)' % name)
    plt.xlabel('Dates')
    plt.ylabel('Depth (m)')
    
    # Add model data to axis (RAPID approx)
    fig.add_subplot(3,1,2)
    plt.pcolormesh(dts, -z, sf_rapid.transpose(), 
                   vmin=cmin,vmax=cmax,cmap=cmap, norm=norm)
    plt.colorbar(orientation='vertical')
    plt.title('Overturning streamfunction at 26N in %s (RAPID approx)' % name)
    plt.xlabel('Dates')
    plt.ylabel('Depth (m)')

    # Add optional observed data to axis
    if obs is not None:
        fig.add_subplot(3,1,3)
        plt.pcolormesh(obs.dates, -obs.z, obs.sf.transpose(), 
                       vmin=cmin,vmax=cmax, cmap=cmap, norm=norm)
        plt.colorbar(orientation='vertical')
        plt.title('Overturning streamfunction at 26N from RAPID array')
        plt.xlabel('Dates')
        plt.ylabel('Depth (m)')
            
    # Save plot
    plt.tight_layout()
    savef = basename + 'overturning_streamfunction_at_26n_hovmoller.png'
    print 'SAVING: %s' % savef
    fig.savefig(savef, resolution=300)
    plt.close() 
    
    
def plot_zonal_mean_temperature(trans, name='simulated', basename='', obs=None,lw=4):
    """ Plot basin-wide zonal mean otential temperature profile """

    # Extract variables from data objects
    z = trans.variables['depth'][:]
    t_basin = trans.variables['t_basin'][:].mean(axis=0)
    z_scaled = np.sqrt(z)
    
    # Add data to axis
    fig = plt.figure(figsize=(6,8))
    plt.plot(t_basin, -z_scaled, color=c1, linewidth=lw, label=name)

    # Add optional observed data
    if obs is not None:
        z_obs = obs.z
        z_obs_scaled = np.sqrt(z_obs)
        t_basin_obs = obs.t_basin.mean(axis=0)
        plt.plot(t_basin_obs, -z_obs_scaled, '-k', linewidth=lw, label='RAPID observations')
        
    # Annotate plot
    zticks, zlabels = get_zscale_ticks()
    plt.title('zonal mean temperature profile at 26N')
    plt.xlabel('degC')
    plt.ylabel('Depth (m)')
    plt.yticks(zticks, zlabels)
    plt.legend(loc='best', fontsize=8)   
    
    # Save plot
    plt.tight_layout()
    savef = basename + 'zonal_mean_temperature_at_26n.png'
    print 'SAVING: %s' % savef
    fig.savefig(savef, resolution=300)
    plt.close()


def get_zscale_ticks():
    """ Return zticks and zlabels """
    z = (np.arange(16)* 5)**2
    zticks = - np.sqrt(np.abs(z))
    zlabels = [str(zlab) for zlab in z]
    
    return zticks, zlabels


def plot_fc_transport_profile(trans, name='simulated', basename='', obs=None,lw=4):
    """ Plot Florida Current transport profile """

    # Add model data to axis
    fig = plt.figure(figsize=(6,8))
    z = trans.variables['depth'][:]
    z_scaled = np.sqrt(z)
    v_fc = trans.variables['v_fc'][:].mean(axis=0)
    plt.plot(v_fc, -z_scaled, color=c2, linewidth=lw, label=name)

    # add optional obs data to axis
    if obs is not None:
        z_obs = obs.z
        z_obs_scaled = np.sqrt(z_obs)
        v_fc_obs = obs.v_fc.mean(axis=0)
        plt.plot(v_fc_obs, -z_obs_scaled, '-k', linewidth=lw, 
                 label='RAPID observations')
 
    # Annotate plot
    zticks, zlabels = get_zscale_ticks()
    plt.title('Florida current transport profile at 26N')
    plt.xlabel('Sv/m')
    plt.ylabel('Depth (m)')
    plt.yticks(zticks, zlabels)
    plt.legend(loc='best', fontsize=8)   
    
    # Save plot
    plt.tight_layout()
    savef = basename + 'florida_current_transport_profile_at_26n.png'
    print 'SAVING: %s' % savef
    fig.savefig(savef, resolution=300)
    plt.close()
    
    
def plot_transport_profile(trans, name='simulated', basename='', obs=None,lw=4):
    """ Plot basin-wide zonal transport profile """

    # Add model data to axis
    fig = plt.figure(figsize=(6,8))
    z = trans.variables['depth'][:]
    z_scaled = np.sqrt(z)
    v_basin_rapid = trans.variables['v_basin_rapid'][:].mean(axis=0)
    v_basin_model = trans.variables['v_basin_model'][:].mean(axis=0)
    
    plt.plot(v_basin_model, -z_scaled, color=c1, linewidth=lw, 
             label='%s (model velocities)' %name)
    plt.plot(v_basin_rapid, -z_scaled, color=c2, linewidth=lw, 
             label='%s (RAPID approx)' %name)

    # Add optional obs data to axis
    if obs is not None:
        z_obs = obs.z
        v_basin_obs = obs.v_basin.mean(axis=0)
        z_obs_scaled = np.sqrt(z_obs)
        plt.plot(v_basin_obs, -z_obs_scaled, '-k', linewidth=lw, 
                 label='RAPID observations')
   
    # Annotate plot
    zticks, zlabels = get_zscale_ticks()
    plt.title('Basin-wide transport profile')
    plt.xlabel('Sv/m')
    plt.ylabel('Depth (m)')
    plt.yticks(zticks, zlabels)
    plt.legend(loc='best', fontsize=8)   
    
    # Save plot
    plt.tight_layout()
    savef = basename + 'basinwide_transport_profile_at_26n.png'
    print 'SAVING: %s' % savef
    fig.savefig(savef, resolution=300)
    plt.close()
    
    
def plot_volume_components(trans, basename='', name='simulated', obs_vol=None, obs_fc=None,lw=4):
    """ Plot volume transport component time series """
    
    # Add model data to sub-axis
    fig = plt.figure(figsize=(8,12))
    fig.add_subplot(2,1,1)

    dts = utils.get_ncdates(trans)
    fc = trans.variables['fc'][:]
    ekman = trans.variables['ekman'][:]
    umo = trans.variables['umo'][:]
    moc = trans.variables['moc_rapid'][:]

    fc_label = 'Florida current (%6.1f Sv)' % (fc.mean())
    ek_label = 'Ekman transport (%6.1f Sv)' % (ekman.mean())
    umo_label = 'Upper-mid ocean (%4.1f Sv)' % (umo.mean())
    moc_label = 'MOC (%6.1f Sv)' % (moc.mean())

    plt.plot(dts, fc, linewidth=lw, color=c1, label=fc_label)
    plt.plot(dts, ekman, linewidth=lw, color=c2, label=ek_label)
    plt.plot(dts, umo, linewidth=lw, color=c3, label=umo_label)
    plt.plot(dts, moc, linewidth=lw, color='k', label=moc_label)

    plt.xlabel('Date')
    plt.ylim([-50,50])
    plt.ylabel('Sverdrups')
    plt.title('Overturning components at 26N in %s' % name)
    plt.legend(loc=8, fontsize=8, ncol=2)


    # Add optional observational data to sub-axis
    if (obs_vol is not None) or (obs_fc is not None):
        fig.add_subplot(2,1,2)
        
        if obs_vol is not None:
            fc_obs_label = 'Florida current (%6.1f Sv)' % (obs_vol.fc.mean())
            ek_obs_label = 'Ekman transport (%6.1f Sv)' % (obs_vol.ekman.mean())
            umo_obs_label = 'Upper-mid ocean (%4.1f Sv)' % (obs_vol.umo.mean())
            moc_obs_label = 'MOC (%6.1f Sv)' % (obs_vol.moc.mean())
        
            plt.plot(obs_vol.dates, obs_vol.fc, linewidth=lw, color=c1, label=fc_obs_label)
            plt.plot(obs_vol.dates, obs_vol.ekman, linewidth=lw, color=c2, label=ek_obs_label)
            plt.plot(obs_vol.dates, obs_vol.umo, linewidth=lw, color=c3, label=umo_obs_label)
            plt.plot(obs_vol.dates, obs_vol.moc, linewidth=lw, color='k', label=moc_obs_label)

        if obs_fc is not None:
            fc_obs_ext_label = 'FC - extended (%6.1f Sv)' % (np.nanmean(obs_fc.fc))
            plt.plot(obs_fc.dates, obs_fc.fc, linewidth=lw, color=c4, label=fc_obs_ext_label)
            

        plt.xlabel('Date')
        plt.ylim([-50,50])
        plt.ylabel('Sverdrups')
        plt.title('Overturning components at 26N in RAPID observations')
        plt.legend(loc=8, fontsize=8, ncol=2)

    # Save plot
    plt.tight_layout()
    savef = basename + 'volume_transport_components_at_26n.png'
    print 'SAVING: %s' % savef
    fig.savefig(savef, resolution=300)
    plt.close()
    
        
def linreg_vol_vs_oht(vol,oht):
    """ Return linear regression of volume vs heat transports """
    if len(vol) > 1:
        slope, intercept, r_value, p_value, std_err =  stats.linregress(vol,oht)
        oht_model = vol * slope + intercept
        label = '(%5.3f PW/Sv)' % slope
    else:
        oht_model = None
        label = ''

    return oht_model, label

        
def plot_moc_vs_oht(trans, basename='', name='simulated', obs_vol=None, obs_oht=None):
    """ Plot total overturning vs geometric heat transports """
    
    # Add model data to axis (model v)
    fig = plt.figure(figsize=(15,5))
    fig.add_subplot(1,3,1)

    dts = utils.get_ncdates(trans)
    moc_model = trans.variables['moc_rapid'][:]
    q_sum_model = trans.variables['q_sum_model'][:]
    q_gyre_model = trans.variables['q_gyre_model'][:]
    q_ot_model = trans.variables['q_ot_model'][:]
    
    q_sum_model_lin, q_sum_model_label = linreg_vol_vs_oht(moc_model, q_sum_model)
    q_gyre_model_lin, q_gyre_model_label = linreg_vol_vs_oht(moc_model, q_gyre_model)
    q_ot_model_lin, q_ot_model_label = linreg_vol_vs_oht(moc_model, q_ot_model)
    
    plt.plot(moc_model, q_sum_model,'x', color='k', label='total %s' % q_sum_model_label)
    plt.plot(moc_model, q_ot_model,'x', color=c1, label='overturning %s' % q_ot_model_label)
    plt.plot(moc_model, q_gyre_model,'x', color=c2, label='gyre %s' % q_gyre_model_label)

    if q_sum_model_lin is not None:
        plt.plot(moc_model, q_sum_model_lin,'-', color='k')
        plt.plot(moc_model, q_ot_model_lin,'-', color=c1)
        plt.plot(moc_model, q_gyre_model_lin,'-', color=c2)
    
    plt.xlabel('MOC (Sv)')
    plt.ylabel('Heat transport (PW)')
    plt.title('MOC vs OHT in %s (model velocities)' % name)
    plt.legend(loc='best', fontsize=8, )


    # Add model data to axis (RAPID approx)
    fig.add_subplot(1,3,2)
    
    moc_rapid = trans.variables['moc_rapid'][:]
    q_sum_rapid = trans.variables['q_sum_rapid'][:]
    q_gyre_rapid = trans.variables['q_gyre_rapid'][:]
    q_ot_rapid = trans.variables['q_ot_rapid'][:]

    q_sum_rapid_lin, q_sum_rapid_label = linreg_vol_vs_oht(moc_rapid, q_sum_rapid)
    q_gyre_rapid_lin, q_gyre_rapid_label = linreg_vol_vs_oht(moc_rapid, q_gyre_rapid)
    q_ot_rapid_lin, q_ot_rapid_label = linreg_vol_vs_oht(moc_rapid, q_ot_rapid)

    plt.plot(moc_rapid, q_sum_rapid,'x', color='k', label='total %s' % q_sum_rapid_label)
    plt.plot(moc_rapid, q_ot_rapid,'x', color=c1, label='overturning %s' % q_ot_rapid_label)
    plt.plot(moc_rapid, q_gyre_rapid,'x', color=c2, label='gyre %s' % q_gyre_rapid_label)

    if q_sum_rapid_lin is not None:
        plt.plot(moc_rapid, q_sum_rapid_lin,'-', color='k')
        plt.plot(moc_rapid, q_ot_rapid_lin,'-', color=c1)
        plt.plot(moc_rapid, q_gyre_rapid_lin,'-', color=c2)

    plt.xlabel('MOC (Sv)')
    plt.ylabel('Heat transport (PW)')
    plt.title('MOC vs OHT in %s (RAPID approx)' % name)
    plt.legend(loc='best', fontsize=8, )


    # Add optional obs data to axis
    if (obs_oht is not None) and (obs_vol is not None):
        fig.add_subplot(1,3,3)

        mindt, maxdt = utils.get_daterange(obs_vol.dates, obs_oht.dates)
        vol_ind = utils.get_dateind(obs_vol.dates, mindt, maxdt)
        oht_ind = utils.get_dateind(obs_oht.dates, mindt, maxdt)
        q_sum_obs = obs_oht.q_sum[oht_ind]
        q_gyre_obs = obs_oht.q_gyre[oht_ind]
        q_ot_obs = obs_oht.q_ot[oht_ind]
        moc_obs = obs_vol.moc[vol_ind]
        
        q_sum_obs_lin, q_sum_obs_label = linreg_vol_vs_oht(moc_obs, q_sum_obs)
        q_gyre_obs_lin, q_gyre_obs_label = linreg_vol_vs_oht(moc_obs, q_gyre_obs)
        q_ot_obs_lin, q_ot_obs_label = linreg_vol_vs_oht(moc_obs, q_ot_obs)
        
        plt.plot(moc_obs, q_sum_obs,'x', color='k', label='total %s' % q_sum_obs_label)
        plt.plot(moc_obs, q_ot_obs,'x', color=c1, label='overturning %s' % q_ot_obs_label)
        plt.plot(moc_obs, q_gyre_obs,'x', color=c2, label='gyre %s' % q_gyre_obs_label)
        
        if q_sum_obs_lin is not None:
            plt.plot(moc_obs, q_sum_obs_lin,'-', color='k')       
            plt.plot(moc_obs, q_ot_obs_lin,'-', color=c1)
            plt.plot(moc_obs, q_gyre_obs_lin,'-', color=c2)

        plt.xlabel('MOC (Sv)')
        plt.ylabel('Heat transport (PW)')
        plt.title('MOC vs OHT in RAPID observations')
        plt.legend(loc='best', fontsize=8, )

    # Save plot
    plt.tight_layout()
    savef = basename + 'moc_vs_heat_transports_at_26n.png'
    print 'SAVING: %s' % savef
    fig.savefig(savef, resolution=300)
    plt.close()


def plot_vol_vs_heat_transports(trans, basename='', name='simulated', obs_vol=None, obs_oht=None):
    """ Plot total overturning vs geometric heat transports """
    
    # Extract model variables from data objects
    dts = utils.get_ncdates(trans)
    fc = trans.variables['fc'][:]
    ek = trans.variables['ekman'][:]
    umo = trans.variables['umo'][:]
    q_ek = trans.variables['q_ek'][:]
    q_fc = trans.variables['q_fc'][:]
    q_mo = trans.variables['q_mo'][:]

    q_fc_lin, q_fc_label = linreg_vol_vs_oht(fc, q_fc)
    q_ek_lin, q_ek_label = linreg_vol_vs_oht(ek, q_ek)
    q_mo_lin, q_mo_label = linreg_vol_vs_oht(umo, q_mo)

    
    # Extract obs variables from data objects
    if (obs_vol is not None) and (obs_oht is not None):
        mindt, maxdt = utils.get_daterange(obs_vol.dates, obs_oht.dates)
        vol_ind = utils.get_dateind(obs_vol.dates, mindt, maxdt)
        oht_ind = utils.get_dateind(obs_oht.dates, mindt, maxdt)
        ek_obs = obs_vol.ekman[vol_ind]
        fc_obs = obs_vol.fc[vol_ind]
        umo_obs = obs_vol.umo[vol_ind]
        q_ek_obs = obs_oht.q_ek[oht_ind]
        q_fc_obs = obs_oht.q_fc[oht_ind]
        q_mo_obs = obs_oht.q_mo[oht_ind]
        
        q_fc_obs_lin, q_fc_obs_label = linreg_vol_vs_oht(fc_obs, q_fc_obs)
        q_ek_obs_lin, q_ek_obs_label = linreg_vol_vs_oht(ek_obs, q_ek_obs)
        q_mo_obs_lin, q_mo_obs_label = linreg_vol_vs_oht(umo_obs, q_mo_obs)
    
    # Plot ekman 
    fig = plt.figure(figsize=(15,5))
    fig.add_subplot(1,3,1)

    plt.plot(ek, q_ek,'x', color=c1, label='%s %s' % (name, q_ek_label))
    if q_ek_lin is not None:
        plt.plot(ek, q_ek_lin,'-', color=c1)
   
    if (obs_vol is not None) and (obs_oht is not None):
        plt.plot(ek_obs, q_ek_obs,'x', color='k', label='RAPID observations %s' % q_ek_obs_label)
        if q_ek_obs_lin is not None:
            plt.plot(ek_obs, q_ek_obs_lin,'-', color='k')
        
    plt.xlabel('Volume transport (Sv)')
    plt.ylabel('Heat transport (PW)')
    plt.title('Ekman volume vs heat transport')
    plt.legend(loc='best', fontsize=8)
    
    # Plot florida current 
    fig.add_subplot(1,3,2)

    plt.plot(fc, q_fc,'x', color=c1, label='%s %s' % (name, q_fc_label))
    if q_fc_lin is not None:
        plt.plot(fc, q_fc_lin,'-', color=c1)

    if (obs_vol is not None) and (obs_oht is not None):
        plt.plot(fc_obs, q_fc_obs,'x', color='k', label='RAPID observations %s' % q_fc_obs_label)
        if q_fc_obs_lin is not None:
            plt.plot(fc_obs, q_fc_obs_lin,'-', color='k')
   
    plt.xlabel('Volume transport (Sv)')
    plt.ylabel('Heat transport (PW)')
    plt.title('Florida current volume vs heat transport')
    plt.legend(loc='best', fontsize=8)
    
    # Plot mid-ocean 
    fig.add_subplot(1,3,3)
    
    plt.plot(umo, q_mo,'x', color=c1, label='%s %s' % (name, q_mo_label))
    if q_mo_lin is not None:
        plt.plot(umo, q_mo_lin,'-', color=c1)
   
    if (obs_vol is not None) and (obs_oht is not None):
        plt.plot(umo_obs, q_mo_obs,'x', color='k', label='RAPID observations %s' % q_mo_obs_label)
        if q_mo_obs_lin is not None:
            plt.plot(umo_obs, q_mo_obs_lin,'-', color='k')

    plt.xlabel('Volume transport (Sv)')
    plt.ylabel('Heat transport (PW)')
    plt.title('Mid-ocean volume vs heat transport')
    plt.legend(loc='best', fontsize=8)

    # Save plot
    plt.tight_layout()
    savef = basename + 'volume_vs_heat_transports_at_26n.png'
    print 'SAVING: %s' % savef
    fig.savefig(savef, resolution=300)
    plt.close()

    
def plot_geometric_heat_components(trans, basename='', name='simulated', obs=None,lw=4):
    """ Plot geometric heat transport components """
      
    # Add model data to sub-axis (model v)
    fig = plt.figure(figsize=(8,12))
    fig.add_subplot(3,1,1)

    dts = utils.get_ncdates(trans)
    q_sum_model = trans.variables['q_sum_model'][:]
    q_gyre_model = trans.variables['q_gyre_model'][:]
    q_ot_model = trans.variables['q_ot_model'][:]
    q_net_model = trans.variables['q_net_model'][:]

    q_sum_model_label = 'Total  (%4.2f PW)' % (q_sum_model.mean())
    q_gyre_model_label = 'Gyre (%4.2f PW)' % (q_gyre_model.mean())
    q_ot_model_label = 'Overturning  (%4.2f PW)' % (q_ot_model.mean())
    q_net_model_label = 'Net (%4.2f PW)' % (q_net_model.mean())
    
    plt.plot(dts, q_sum_model, linewidth=lw, color='k', label=q_sum_model_label)
    plt.plot(dts, q_ot_model, linewidth=lw, color=c1, label=q_ot_model_label)
    plt.plot(dts, q_gyre_model, linewidth=lw, color=c2, label=q_gyre_model_label)
    plt.plot(dts, q_net_model, linewidth=lw, color=c3, label=q_net_model_label)

    plt.xlabel('Date')
    plt.ylim([-.5,2.5])
    plt.ylabel('PW')
    plt.title('Heat transports at 26N in %s (model velocities)'  % name)
    plt.legend(loc=8, fontsize=8, ncol=2)


    # Add model data to sub-axis (RAPID approx)
    fig.add_subplot(3,1,2)
       
    q_sum_rapid = trans.variables['q_sum_rapid'][:]
    q_gyre_rapid = trans.variables['q_gyre_rapid'][:]
    q_ot_rapid = trans.variables['q_ot_rapid'][:]
    q_net_rapid = trans.variables['q_net_rapid'][:]
    
    q_sum_rapid_label = 'Total (%4.2f PW)' % (q_sum_rapid.mean())
    q_gyre_rapid_label = 'Gyre (%4.2f PW)' % (q_gyre_rapid.mean())
    q_ot_rapid_label = 'Overturning (%4.2f PW)' % (q_ot_rapid.mean())
    q_net_rapid_label = 'Net (%4.2f PW)' % (q_net_rapid.mean())

    plt.plot(dts, q_sum_rapid, linewidth=lw, color='k', label=q_sum_rapid_label)
    plt.plot(dts, q_ot_rapid, linewidth=lw, color=c1, label=q_ot_rapid_label)
    plt.plot(dts, q_gyre_rapid, linewidth=lw, color=c2, label=q_gyre_rapid_label)
    plt.plot(dts, q_net_rapid, linewidth=lw, color=c3, label=q_net_rapid_label)
        
    plt.xlabel('Date')
    plt.ylim([-.5,2.5])
    plt.ylabel('PW')
    plt.title('Heat transports at 26N in %s (RAPID approx)'  % name)
    plt.legend(loc=8, fontsize=8, ncol=2)


    # Add optional observational data to sub-axis
    if obs is not None:
        fig.add_subplot(3,1,3)
        q_sum_obs_label = 'Total (%4.2f PW)' % (obs.q_sum.mean())
        q_gyre_obs_label = 'Gyre (%4.2f PW)' % (obs.q_gyre.mean())
        q_ot_obs_label = 'Overturning (%4.2f PW)' % (obs.q_ot.mean())
        
        plt.plot(obs.dates, obs.q_sum, linewidth=lw, color='k',
                 label=q_sum_obs_label)
        plt.plot(obs.dates, obs.q_ot, linewidth=lw, color=c1, 
                 label=q_ot_obs_label)
        plt.plot(obs.dates, obs.q_gyre, linewidth=lw, color=c2, 
                 label=q_gyre_obs_label)
    
        plt.xlabel('Date')
        plt.ylim([-.5,2.5])
        plt.ylabel('PW')
        plt.title('Heat transports at 26N in RAPID observations')
        plt.legend(loc=8, fontsize=8, ncol=2)
   

    # Save plot
    plt.tight_layout()
    savef = basename + 'heat_transports_geometric_decomposition_at_26n.png'
    print 'SAVING: %s' % savef
    fig.savefig(savef, resolution=300)
    plt.close()
    
    
def plot_rapid_heat_components(trans, basename='', name='simulated', obs=None,lw=4):
    """ Plot RAPID heat transport components """
    
    # Add model data to sub-axis
    fig = plt.figure(figsize=(8,12))
    fig.add_subplot(2,1,1)

    dts = utils.get_ncdates(trans)
    q_sum = trans.variables['q_sum_rapid'][:]
    q_ek = trans.variables['q_ek'][:]
    q_fc = trans.variables['q_fc'][:]
    q_geoint = trans.variables['q_geoint'][:]
    q_eddy = trans.variables['q_eddy'][:]
    q_wbw = trans.variables['q_wbw'][:]

    q_sum_label = 'Total (%4.2f PW)' % (q_sum.mean())
    q_ek_label = 'Ekman (%4.2f PW)' % (q_ek.mean())
    q_fc_label = 'Florida current (%4.2f PW)' % (q_fc.mean())
    q_geoint_label = 'Geostrophic interior (%4.2f PW)' % (q_geoint.mean())
    q_wbw_label = 'WBW (%4.2f PW)' % (q_wbw.mean())
    q_eddy_label = 'Eddies (%4.2f PW)' % (q_eddy.mean())

    plt.plot(dts, q_sum, linewidth=lw, color='k', label=q_sum_label)
    plt.plot(dts, q_ek, linewidth=lw, color=c1, label=q_ek_label)
    plt.plot(dts, q_fc, linewidth=lw, color=c2, label=q_fc_label)
    plt.plot(dts, q_wbw, linewidth=lw, color=c3, label=q_wbw_label)
    plt.plot(dts, q_geoint, linewidth=lw, color=c4, label=q_geoint_label)
    plt.plot(dts, q_eddy, linewidth=lw, color='0.5', label=q_eddy_label)
    plt.xlabel('Date')
    plt.ylim([-4,4])
    plt.ylabel('PW')
    plt.title('Heat transports relative to 0C at 26N in %s'  % name)
    plt.legend(loc=8, fontsize=8, ncol=2)


    # Add observational data to sub-axis
    if obs is not None:
        fig.add_subplot(2,1,2)
    
        q_sum_obs_label = 'Total (%4.2f PW)' % (obs.q_sum.mean())
        q_ek_obs_label = 'Ekman (%4.2f PW)' % (obs.q_ek.mean())
        q_fc_obs_label = 'Florida current (%4.2f PW)' % (obs.q_fc.mean())
        q_geoint_obs_label = 'Geostrophic interior (%4.2f PW)' % (obs.q_geoint.mean())
        q_wbw_obs_label = 'WBW (%4.2f PW)' % (obs.q_wbw.mean())
        q_eddy_obs_label = 'Eddies (%4.2f PW)' % (obs.q_eddy.mean())
    
        plt.plot(obs.dates, obs.q_sum, linewidth=lw, color='k',
                 label=q_sum_obs_label)
        plt.plot(obs.dates, obs.q_ek, linewidth=lw, color=c1, 
                 label=q_ek_obs_label)
        plt.plot(obs.dates, obs.q_fc, linewidth=lw, color=c2, 
                 label=q_fc_obs_label)
        plt.plot(obs.dates, obs.q_wbw, linewidth=lw, color=c3, 
                 label=q_wbw_obs_label)
        plt.plot(obs.dates, obs.q_geoint, linewidth=lw, color=c4,
                 label=q_geoint_obs_label)
        plt.plot(obs.dates, obs.q_eddy, linewidth=lw, color='0.5', 
                 label=q_eddy_obs_label)

        plt.xlabel('Date')
        plt.ylim([-4,4])
        plt.ylabel('PW')
        plt.title('Heat transports relative to 0C at 26N in RAPID observations')
        plt.legend(loc=8, fontsize=8, ncol=2)
     
    # Save plot
    plt.tight_layout()
    savef = basename + 'heat_transports_rapid_decomposition_at_26n.png'
    print 'SAVING: %s' % savef
    fig.savefig(savef, resolution=300)
    plt.close()
    

def plot_diagnostics(trans, name='simulated', outdir='./', date_format='%Y%m%d',
                     obs_vol=None, obs_oht=None, obs_sf=None, obs_fc=None,lw=4):
    """ Plot volume and heat transport diagnostics against RAPID observations """
    
    # Create basename for output files
    dts = utils.get_ncdates(trans)
    basename = utils.get_savename(outdir, name, dts, date_format,suffix='_')
    
    # Plot data
    plot_streamfunctions(trans, basename=basename, name=name, obs=obs_sf,lw=lw)
    plot_streamfunction_hovmollers(trans, basename=basename, name=name, obs=obs_sf)
    plot_volume_components(trans, basename=basename, name=name, obs_vol=obs_vol, obs_fc=obs_fc,lw=lw)
    plot_rapid_heat_components(trans, basename=basename, name=name, obs=obs_oht,lw=lw)
    plot_geometric_heat_components(trans, basename=basename, name=name, obs=obs_oht,lw=lw)
    plot_zonal_mean_temperature(trans, basename=basename, name=name, obs=obs_oht,lw=lw)
    plot_transport_profile(trans, basename=basename, name=name, obs=obs_oht,lw=lw)
    plot_fc_transport_profile(trans, basename=basename, name=name, obs=obs_oht,lw=lw)
    plot_moc_vs_oht(trans, basename=basename, name=name, obs_vol=obs_vol, obs_oht=obs_oht)
    plot_vol_vs_heat_transports(trans, basename=basename, name=name, obs_vol=obs_vol, obs_oht=obs_oht)


