"""
Module containing code to work with ocean transports

"""

import numpy as np
import utils
import copy

import output

# Constants
G = 9.81          # Gravitational acceleration (m/s2)
ROT = 7.292116E-5 # Rotation rate of earth (rad/s) 
RHO_REF = 1025.   # Reference density for sea water (kg/m3)
CP = 3985         # Specific heat capacity of sea water (J/kg/K)

    
class Transports(object):
    """ Class to interface with volume and heat transport diagnostics """
    
    def __init__(self, v, t_on_v,minind,maxind):
        """ Initialize with velocity and temperature sections """
        
        # Initialize data
        self.name = v.name
        self.v = v.data[:,:,minind:maxind]
        self.t = t_on_v.data[:,:,minind:maxind]
        self.rhocp = RHO_REF * CP
        self.x = v.x[minind:maxind]        
        self.y = v.y[minind:maxind]
        self.z = v.z
        self.dates = v.dates
        self.dz = v.dz
        self.dz_as_data = v.dz_as_data[:,:,minind:maxind]
        self.dx = v.cell_widths
        self.dx_as_data = v.cell_widths_as_data[:,:,minind:maxind]
        self.da = self.dx_as_data * self.dz_as_data
        
        # Set null values for property attributes
        self._avg_t = None
        self._avg_v = None
        self._v_no_net = None
        self._net_transport = None
        self._zonal_avg_v = None
        self._zonal_avg_v_no_net = None
        self._zonal_anom_v = None
        self._zonal_sum_v_no_net = None
        self._zonal_sum_v = None
        self._zonal_avg_t = None
        self._zonal_anom_t = None
        self._streamfunction = None
        self._oht_total = None
        self._oht_by_net = None
        self._oht_by_horizontal = None
        self._oht_by_overturning = None
            
    def section_avg(self, data, total=False):
        """ Return avg across whole section """
        if total:
            return (data * self.da).sum(axis=(1,2)) 
        else:
            return (data * self.da).sum(axis=(1,2)) / self.da.sum(axis=(1,2))
    
    def zonal_avg(self, data, total=False):
        """ Return zonal avg across section """
        if total:
            return (data * self.dx_as_data).sum(axis=2) 
        else:
            return (data * self.dx_as_data).sum(axis=2) / self.dx_as_data.sum(axis=2)
    
    @property
    def streamfunction(self):
        """ Return contribution to basin-wide overturning streamfunction """
        if self._streamfunction is None:
            self._streamfunction = np.cumsum(
                (self.v.filled(0) * self.da.filled(0)).sum(axis=2), axis=1)
        return self._streamfunction
    
    @property 
    def avg_t(self):
        """ Return section average temperature """
        if self._avg_t is None:
            self._avg_t = self.section_avg(self.t)
        return self._avg_t
    
    @property 
    def avg_v(self):
        """ Return section average velocity """
        if self._avg_v is None:
            self._avg_v = self.section_avg(self.v) 
        return self._avg_v 
    
    @property
    def v_no_net(self):
        """ Return velocity after removing net transport through section """
        if self._v_no_net is None:
            self._v_no_net = self.v - self.avg_v[:,np.newaxis,np.newaxis]
        return self._v_no_net
    
    @property 
    def net_transport(self):
        """ Return net transport through section """
        if self._net_transport is None:
            self._net_transport = self.section_avg(self.v, total=True) 
        return self._net_transport
       
    @property
    def zonal_avg_v_no_net(self):
        """ Return zonal mean of v_no_net """
        if self._zonal_avg_v_no_net is None:
            self._zonal_avg_v_no_net = self.zonal_avg(self.v_no_net)
        return self._zonal_avg_v_no_net
    
    @property
    def zonal_avg_v(self):
        """ Return zonal mean of v """
        if self._zonal_avg_v is None:
            self._zonal_avg_v = self.zonal_avg(self.v)
        return self._zonal_avg_v
    
    @property
    def zonal_avg_t(self):
        """ Return zonal mean temperature profile """
        if self._zonal_avg_t is None:
            self._zonal_avg_t = self.zonal_avg(self.t)
        return self._zonal_avg_t
    
    @property
    def zonal_anom_v(self):
        """ Return velocity anomalies relative to zonal mean profile """
        if self._zonal_anom_v is None:
            self._zonal_anom_v = self.v_no_net - self.zonal_avg_v_no_net[:,:,np.newaxis]
        return self._zonal_anom_v
        
    @property
    def zonal_anom_t(self):
        """ Return temperature anomalies relative to zonal mean profile """
        if self._zonal_anom_t is None:
            self._zonal_anom_t = self.t - self.zonal_avg_t[:,:,np.newaxis]
        return self._zonal_anom_t
    
    @property
    def zonal_sum_v(self):
        """ Return zonal integral of v """
        if self._zonal_sum_v is None:
            self._zonal_sum_v = self.zonal_avg(self.v, total=True)
        return self._zonal_sum_v
              
    @property
    def zonal_sum_v_no_net(self):
        """ Return zonal integral of v_no_net """
        if self._zonal_sum_v_no_net is None:
            self._zonal_sum_v_no_net = self.zonal_avg(self.v_no_net, total=True)
        return self._zonal_sum_v_no_net
                
    @property
    def oht_by_net(self):
        """ Return heat transport by net transport through section """
        if self._oht_by_net is None:
            self._oht_by_net = self.net_transport * self.avg_t * self.rhocp
        return self._oht_by_net
    
    @property
    def oht_total(self):
        """ Return total heat transport through section """
        if self._oht_total is None:
            self._oht_total = self.section_avg(self.v * self.t, total=True) * self.rhocp
        return self._oht_total    
    
    @property
    def oht_by_horizontal(self):
        """ Return heat transport by horizontal circulation """
        if self._oht_by_horizontal is None:
            self._oht_by_horizontal = self.section_avg(self.zonal_anom_v * self.zonal_anom_t,
                                                    total=True) * self.rhocp            
        return self._oht_by_horizontal   

    @property
    def oht_by_overturning(self):
        """ Return heat transport by local overturning circulation """
        if self._oht_by_overturning is None:
            self._oht_by_overturning = (self.zonal_sum_v_no_net * self.zonal_avg_t *
                                        self.dz[np.newaxis,:]).sum(axis=1) * self.rhocp
        return self._oht_by_overturning   


def calc_transports_from_sections(config, v, tau, t_on_v, s_on_v):
    """
    High-level routine to call transport calculations and return
    integrated transports on RAPID section as netcdf object
    
    """ 
    # Extract sub-section boundaries
    fc_minlon = config.getfloat('options','fc_minlon')   # Minimum longitude for Florida current
    fc_maxlon = config.getfloat('options','fc_maxlon')   # Longitude of Florida current/WBW boundary
    wbw_maxlon = config.getfloat('options','wbw_maxlon') # Longitude of WBW/gyre boundary
    int_maxlon = config.getfloat('options','int_maxlon') # Maximum longitude of gyre interior
    
    # Get indices for sub-sections 
    fcmin, fcmax = utils.get_indrange(v.x, fc_minlon, fc_maxlon)     # Florida current
    wbwmin, wbwmax = utils.get_indrange(v.x, fc_maxlon, wbw_maxlon)  # WBW
    intmin, intmax = utils.get_indrange(v.x, wbw_maxlon, int_maxlon) # Gyre interior
    
    # Calculate dynamic heights
    dh = calc_dh(t_on_v, s_on_v)
    
    # Calculate geostrophic transports
    georef = config.getfloat('options', 'georef_level')
    vgeo = calc_vgeo(v, dh, georef=georef)
        
    # Optionally reference geostrophic transports to model velocities
    if config.has_option('options', 'vref_level'):
        vref_level = config.getfloat('options', 'vref_level') 
        vgeo = update_georef(vgeo, v, vref_level)
   
    # Calculate Ekman velocities
    ek_level = config.getfloat('options','ekman_depth')
    ek = calc_ek(v, tau, wbw_maxlon, int_maxlon, ek_level)

    # Use model velocities in fc and WBW regions
    vgeo = merge_vgeo_and_v(vgeo, v, fc_minlon, wbw_maxlon)

    # Apply mass-balance constraints to section
    vgeo = rapid_mass_balance(vgeo, ek, fc_minlon, wbw_maxlon, int_maxlon)

    # Add ekman to geostrophic transports for combined rapid velocities
    vrapid = copy.deepcopy(vgeo)
    vrapid.data = vgeo.data + ek.data
    
    # Get volume and heat transports on each (sub-)section
    fc_trans = Transports(vgeo, t_on_v, fcmin, fcmax)        # Florida current transports
    wbw_trans = Transports(vgeo, t_on_v, wbwmin, wbwmax)     # Western-boundary wedge transports
    int_trans = Transports(vgeo, t_on_v, intmin, intmax)     # Gyre interior transports
    ek_trans = Transports(ek, t_on_v, intmin, intmax)        # Ekman transports
    model_trans = Transports(v, t_on_v, fcmin, intmax)       # Total section transports using model velocities
    rapid_trans = Transports(vrapid, t_on_v, fcmin, intmax)  # Total section transports using RAPID approximation

    # Create netcdf object for output/plotting
    trans = output.create_netcdf(config,rapid_trans, model_trans, fc_trans, 
                                 wbw_trans, int_trans, ek_trans)
    
    return trans
    

def calc_dh(t_on_v, s_on_v):
    """ 
    Return ZonalSections containing dynamic heights calculated from 
    from temperature and salinity interpolated onto velocity boundaries. 
   
    """
    # Calculate in situ density at bounds
    rho = copy.deepcopy(t_on_v)
    rho.data = None # Density not needed at v mid-points
    rho.bounds_data = eos_insitu(t_on_v.bounds_data, s_on_v.bounds_data,
                                 t_on_v.z_as_bounds_data)

    # Calculate dynamic height relative to a reference level
    dh = copy.deepcopy(rho)
    rho_anom = (rho.bounds_data - RHO_REF) / RHO_REF
    # Depth axis reversed for integral from sea-floor. 
    dh.bounds_data = np.cumsum((rho_anom * rho.dz_as_bounds_data)[:,::-1,:],
                                axis=1)[:,::-1,:]

    return dh


def calc_vgeo(v, dh, georef=4750.):
    """ 
    Return ZonalSections containing geostrophic velocities 
    relative to specified reference level. 
    
    """
    vgeo = copy.deepcopy(v) # Copy velocity data structure
     
    for nprof in range(len(vgeo.x)): # Loop through profiles
        if not v.mask[:,:,nprof].all():
            
            # Extract depth and dynamic height profiles at bounds
            z = dh.z
            z1 = dh.z_as_bounds_data[:,:,nprof]
            z2 = dh.z_as_bounds_data[:,:,nprof+1]
            dh1 = dh.bounds_data[:,:,nprof]
            dh2 = dh.bounds_data[:,:,nprof+1]
        
            # Coriolis parameter at profile location
            corf = 2 * ROT * np.sin(np.pi * (vgeo.y[nprof]/180.) ) 
            
            # cell width along section
            dx = vgeo.cell_widths[nprof]
            
            # Clip reference depth using ocean floor.
            maxz = np.min([z1.max(),z2.max()])
            zref = min(georef, maxz)
        
            # Adjust dh to new reference level
            zind = utils.find_nearest(z,zref)
            dh1 -= dh1[:,zind]
            dh2 -= dh2[:,zind]

            # Calculate geostrophic velocity
            vgeo_profile = (-1. * (G / corf) * ( (dh2 - dh1) / dx))
            vgeo.data[:,:,nprof] = vgeo_profile
            
    return vgeo


def update_georef(vgeo, v, vref_level):
    """ 
    Return vgeo after updating geostrophic reference depth by constraining
    velocities in vgeo to match those in v at the specified depth.
    
    """
    vgeodat = vgeo.data.filled(0)
    vdat = v.data.filled(0)
    zind = utils.find_nearest(v.z,vref_level)
    vadj = np.ones_like(vgeodat) * (vdat[:,zind,:] - vgeodat[:,zind,:])
    vgeo.data = np.ma.MaskedArray(vgeo.data + vadj, mask=vgeo.mask)
    
    return vgeo


def calc_ek(v, tau, minlon, maxlon, ek_level):
    """ Return ZonalSections containing Ekman velocities """

    # Copy velocity data structure
    ek = copy.deepcopy(v)
    ek.data = np.zeros_like(v.data)
    
    # Get indices for gyre interior
    intmin, intmax = utils.get_indrange(tau.x, minlon, maxlon)

    # Calculate depth-integrated Ekman transports
    dx = tau.cell_widths_as_data[:,intmin:intmax]
    lats = tau.y[intmin:intmax]
    taux = tau.data[:,intmin:intmax]
    corf = 2 * ROT * np.sin(np.pi * (lats / 180.) ) 
    ek_trans = ((-1. *  taux / (corf * RHO_REF)) * dx ).sum(axis=1)

    # Calculate average velocity over ekman layer
    ek_minind, ek_maxind = utils.get_indrange(v.z, 0, ek_level)
    dz = v.dz_as_data[0,ek_minind:ek_maxind,intmin:intmax]
    dx = v.cell_widths_as_data[0,ek_minind:ek_maxind,intmin:intmax]
    ek_area = (dx * dz).sum()
    ek.data[:,ek_minind:ek_maxind,intmin:intmax] = ek_trans[:,np.newaxis, np.newaxis] / ek_area
    ek.data = np.ma.MaskedArray(ek.data, mask=v.mask)

    return ek    


def merge_vgeo_and_v(vgeo, v, minlon, maxlon):
    """ Return vgeo with velocities from v west of lonbnd """
    minind, maxind = utils.get_indrange(vgeo.x, minlon, maxlon)
    vgeo.data[:,:,minind:maxind] = v.data[:,:,minind:maxind]
    
    return vgeo

    
def section_integral(v, xmin, xmax):
    """ Section integral between x values """
    minind, maxind = utils.get_indrange(v.x, xmin, xmax)
    if v.surface_field:
        dx = v.cell_widths_as_data[:,minind:maxind]
        return np.sum(v.data[:,minind:maxind] * dx, axis=1)
    else:
        da = v.cell_widths_as_data[:,:,minind:maxind] * v.dz_as_data[:,:,minind:maxind]
        return np.sum(v.data[:,:,minind:maxind] * da , axis=(1,2))
        

def rapid_mass_balance(vgeo, ek, minlon, midlon, maxlon):
    """ 
    Return vgeo after applying RAPID-style mass-balance constraint
    as a barotropic velocity over geostrophic interior
    """
    
    # Calculate net transports
    fcwbw_tot = section_integral(vgeo, minlon, midlon)
    ek_tot = section_integral(ek, midlon, maxlon)
    int_tot = section_integral(vgeo, midlon, maxlon)
    net = int_tot + ek_tot + fcwbw_tot

    # Get cell dimensions in gyre interior
    minind,maxind = utils.get_indrange(vgeo.x, midlon, maxlon) 
    dz = vgeo.dz_as_data[0]
    dx = vgeo.cell_widths_as_data[0]
    da = (dx[:,minind:maxind] * dz[:,minind:maxind])

    # Correct geostrophic transports in gyre interior
    corr = net / da.sum()
    vgeo.data[:,:,minind:maxind] = (vgeo.data[:,:,minind:maxind] -
                                    corr[:,np.newaxis,np.newaxis])
    
    return vgeo


def total_mass_balance(v):
    """ Apply mass balance evenly across entire section """
    da =  v.cell_widths_as_data * v.dz_as_data
    v.data = (v.data - ((v.data * da).sum(axis=(1,2)) / 
                        da.sum(axis=(1,2)))[:,np.newaxis,np.newaxis])
    
    return v


def eos_insitu(t, s, p):
    """
    Returns in situ density of seawater as calculated by the NEMO
    routine eos_insitu.f90. Computes the density referenced to
    a specified depth/pressure from potential temperature and salinity
    using the Jackett and McDougall (1994) equation of state.
    
    """
    # Convert to double precision
    ptem = np.double(t)    # potential temperature (celcius)
    psal = np.double(s)    # salintiy (psu)
    depth = np.double(p)   # pressure (decibar) = depth (m)
    rau0 = np.double(1035) # volumic mass of reference (kg/m3)

    # Read into eos_insitu.f90 varnames  
    zrau0r = 1.e0 / rau0
    zt = ptem
    zs = psal
    zh = depth            
    zsr= np.sqrt(np.abs(psal))   # square root salinity

    # compute volumic mass pure water at atm pressure
    zr1 = ( ( ( ( 6.536332e-9*zt-1.120083e-6 )*zt+1.001685e-4)*zt-9.095290e-3 )*zt+6.793952e-2 )*zt+999.842594

    # seawater volumic mass atm pressure
    zr2 = ( ( ( 5.3875e-9*zt-8.2467e-7 ) *zt+7.6438e-5 ) *zt-4.0899e-3 ) *zt+0.824493
    zr3 = ( -1.6546e-6*zt+1.0227e-4 ) *zt-5.72466e-3
    zr4 = 4.8314e-4

    #  potential volumic mass (reference to the surface)
    zrhop = ( zr4*zs + zr3*zsr + zr2 ) *zs + zr1

    # add the compression terms
    ze = ( -3.508914e-8*zt-1.248266e-8 ) *zt-2.595994e-6
    zbw = (  1.296821e-6*zt-5.782165e-9 ) *zt+1.045941e-4
    zb = zbw + ze * zs
    zd = -2.042967e-2
    zc =   (-7.267926e-5*zt+2.598241e-3 ) *zt+0.1571896
    zaw = ( ( 5.939910e-6*zt+2.512549e-3 ) *zt-0.1028859 ) *zt - 4.721788
    za = ( zd*zsr + zc ) *zs + zaw
    zb1 =   (-0.1909078*zt+7.390729 ) *zt-55.87545
    za1 = ( ( 2.326469e-3*zt+1.553190)*zt-65.00517 ) *zt+1044.077
    zkw = ( ( (-1.361629e-4*zt-1.852732e-2 ) *zt-30.41638 ) *zt + 2098.925 ) *zt+190925.6
    zk0 = ( zb1*zsr + za1 )*zs + zkw

    # Caculate density
    prd = (  zrhop / (  1.0 - zh / ( zk0 - zh * ( za - zh * zb ) )  ) - rau0  ) * zrau0r
    rho = (prd*rau0) + rau0

    return rho
