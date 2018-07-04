"""
Module holding code to work with ocean sections 

"""

from math import radians, cos, sin, asin, sqrt
from netCDF4 import MFTime, MFDataset, Dataset, num2date
import numpy as np
import copy
import glob
import pandas as pd

class NetCDF4Error(Exception):
    pass

class MaskError(Exception):
    pass

class ShapeError(Exception):
    pass

class ZonalSections(object):
    """ 
    Class to interface with zonal section data.
    
    ** Assumptions **
    
    Data arrays must be shaped as follows:        
    
    3D fields - data[t,z,y,x]
    2D fields - data[t,y,x]
        
    where the x-coordinate is assumed to be approximately
    zonal, and the y-coordinate is assumed to be approximately
    meridional. 
    
    Some deviations from this assumption are accounted for
    when estimating cell bounds and for along-section
    interpolation of data onto the velocity grid. These 
    operations are designed to work with generalized
    curvilinear coordinate grids (e.g. the NEMO ORCA grid)
    and different grid stencils. 
    
    However, it is assumed throughout the code that 
    sub-sections can be selected using longitude pairs.
    Extracted sections that do not have an x-coordinate of 
    monotonically increasing longitude data will probably
    not behave as expected.
    
    """
    def __init__(self, f, config, section):
        """ Initialize class and read data """
        self.f = f
        self.name = config.get('output', 'name')
        self.section = section.lower()
        self.var = config.get(section, 'var')
        self.xcoord = config.get(section, 'xcoord')
        self.ycoord = config.get(section, 'ycoord')
        self.tcoord = config.get(section, 'tcoord')
        self.i1 = config.getint(section, 'i1')
        self.i2 = config.getint(section, 'i2')
        self.j1 = config.getint(section, 'j1')
        self.j2 = config.getint(section, 'j2')
        self.bounds_data = None

        if config.has_option(section, 'zcoord'):
            self.surface_field = False
            self.zcoord = config.get(section, 'zcoord')
        else:
            self.surface_field = True
            self.zcoord = None
            
        self._read_data()
        self._read_xcoord()
        self._read_ycoord()
        self._read_zcoord()
        self._read_tcoord()

        if config.has_option(section, 'maskf'):
            self.maskvar = config.get(section,'maskvar')
            self.maskf = config.get(section, 'maskf')
            self.maskmdi = config.getfloat(section, 'maskmdi')
            mask = self._read_mask()
            self._update_data_mask(mask)
        elif hasattr(self.data, 'mask'):
            print "%s: using mask information from %s." % (self.var, self.f)
        else:
            raise MaskError('Mask information must be specified')
        
        if config.has_option(section, 'fill_missing_coords'):
            self.fill_missing_coords = config.getboolean(section, 'fill_missing_coords')
        else:
            self.fill_missing_coords = False

        if self.fill_missing_coords:
            self._infill_coords()
                    
        self._apply_meridional_average()


    @property
    def cell_widths(self):
        """ Return width of cells along section """
        
        x = self.xbounds
        y = self.ybounds
        nmax = len(x) - 1
        dists = []
        
        for n in range(nmax):
            dists.append(self._haversine(x[n],y[n],x[n+1],y[n+1]))
            
        return np.array(dists)
    
    @property
    def dz(self):
        """ Return width of cells in z direction """
        if self.surface_field:
            return None
        else:
            return np.abs(np.diff(self.zbounds))
                
    @property
    def xbounds(self):
        """ Return estimated y coordinate bounds """
        return self._get_bounds(self.x)
    
    @property
    def ybounds(self):
        """ Return estimated y coordinate bounds """
        return self._get_bounds(self.y)
    
    @property    
    def zbounds(self):
        """ Return estimated z coordinate bounds """
        if self.surface_field:
            return None
        else:
            return self._get_bounds(self.z)

    @property 
    def mask(self):
        """ Return mask data """
        return self.data.mask
        
    @property 
    def bounds_mask(self):
        """ Return mask at cell bounds """
        if self.surface_field:
            nt,nx = self.data.shape
            bounds_mask = np.ones((nt,len(self.xbounds))) == 1
            bounds_mask[:,:-1] = bounds_mask[:,:-1] & self.mask
            bounds_mask[:,1:] = bounds_mask[:,1:] & self.mask
        else:
            nt,nz,nx = self.data.shape
            bounds_mask = np.ones((nt,nz,len(self.xbounds))) == 1
            bounds_mask[:,:,:-1] = bounds_mask[:,:,:-1] & self.mask
            bounds_mask[:,:,1:] = bounds_mask[:,:,1:] & self.mask

        return bounds_mask
    
    @property
    def z_as_data(self):
        """ Return z with same shape as self.data """
        if (self.surface_field) or (self.data is None):
            return None
        else:
            return (np.ones_like(self.data) * self.z[None,:,None])

    @property
    def z_as_bounds_data(self):
        """ Return z with same shape as self.bounds_data  """
        if (self.surface_field) or (self.bounds_data is None):
            return None
        else:
            return (np.ones_like(self.bounds_data) * self.z[None,:,None])

    @property
    def dz_as_data(self):
        """ Return dz with same shape as self.data  """
        if (self.surface_field) or (self.data is None):
            return None
        else:
            return (np.ones_like(self.data) * self.dz[None,:,None])

    @property
    def dz_as_bounds_data(self):
        """ Return dz with same shape as self.bounds_data  """
        if (self.surface_field) or (self.bounds_data is None):
            return None
        else:
            return (np.ones_like(self.bounds_data) * self.dz[None,:,None])
            
    @property
    def cell_widths_as_data(self):
        """ Return cell_widths with same shape as self.data  """
        if (self.surface_field):
            return (np.ones_like(self.data) * self.cell_widths[None,:])
        else:
            return (np.ones_like(self.data) * self.cell_widths[None,None,:])

    def distance_along_section(self, x0, y0, bounds=False):
        """ Return distance along section relative to (x0,y0) """
        if bounds:
            x = self.xbounds
            y = self.ybounds
        else:
            x = self.x
            y = self.y
        
        dists = []
        
        for n in range(len(x)):
            if n == 0:
                dists.append(self._haversine(x0,y0,x[0],y[0]))
            else:
                dists.append(self._haversine(x[n-1],y[n-1],x[n],y[n]))
            
        return np.cumsum(dists)

    def _interp_surface(self, distances, x0, y0):
        """ Interpolate surface data along section for each time """
        nt,nx = self.data.shape
        idata = np.ma.MaskedArray(np.zeros((nt,len(distances))), mask=True)
        orig_dist = self.distance_along_section(x0,y0)
        
        for tind in range(nt):
            data = self.data[tind,:]
            ind = data.mask == False
            idata[tind,:] = np.interp(distances, orig_dist[ind], data[ind])

        return idata
        
    def _interp_section(self, distances, x0, y0):
        """ Interpolate data along section for each time and depth """
        nt,nz,nx = self.data.shape
        idata = np.ma.MaskedArray(np.zeros((nt,nz,len(distances))), mask=True)
        orig_dist = self.distance_along_section(x0,y0)

        for tind in range(nt):
            for zind in range(nz):
                data = self.data[tind,zind,:]
                ind = data.mask == False
                if ind.any():
                    idata[tind,zind,:] = np.interp(distances,orig_dist[ind],data[ind])

        return idata

    def interp_along_section(self, distances, x0, y0):
        """ Return data interpolated to specified distances relative to (x0,y0) """
        if self.surface_field:
            idata = self._interp_surface(distances, x0, y0)
        else:
            idata = self._interp_section(distances, x0, y0)

        return idata

    def _haversine(self, lon1, lat1, lon2, lat2):
        """ Return distance between lon/lat pairs using haversine formula """
        lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
        dlon = lon2 - lon1 
        dlat = lat2 - lat1 
        a = sin(dlat/2.)**2 + cos(lat1) * cos(lat2) * sin(dlon/2.)**2
        c = 2 * asin(sqrt(a)) 
        m = 6371229. * c # Earth radius as used in NEMO

        return m

    def _get_bounds(self, data):
        """ Estimate bounds using mid-points """
        nmax = len(data) + 1
        bounds = np.zeros(nmax)       
        
        for nbound in range(nmax):
            if nbound == 0:
                bounds[nbound] = data[0] - 0.5 * (data[1] - data[0])
            elif nbound == nmax-1:
                bounds[nbound] = data[nbound-1] + 0.5 * (data[nbound-1] - data[nbound-2])
            else:
                bounds[nbound] = 0.5 * (data[nbound] + data[nbound-1]) 
        
        return bounds

    def _opennc(self, f):
        """ Open netcdf data set either Dataset or MFDataset. """
        if len(glob.glob(f)) > 1:
            try:
                nc = MFDataset(f)
            except ValueError as err:
                print 'netcdf4.MFDataset incompatible with NETCDF4. Try concatenating data into a single file.'
                raise NetCDF4ERROR(err)
        else:
            nc = Dataset(f)

        return nc

    def _apply_meridional_average(self):
        """ Average data across j-indices """
        if self.surface_field:
            self.data = self.data.mean(axis=1)
        else:
            self.data = self.data.mean(axis=2)
        self.x = self.x.mean(axis=0)
        self.y = self.y.mean(axis=0)

    def _update_data_mask(self, mask):
        """ Update mask information for section data """
        # Dummy time indices added again for array broadcasting 
        if self.surface_field:
            mask = np.ones(self.data.shape) * mask[np.newaxis,0]
        else:
            mask = np.ones(self.data.shape) * mask[np.newaxis]
        self.data = np.ma.MaskedArray(self.data, mask=mask)

    def _read_mask(self):
        """ Read mask data from netcdf file(s) """
        nc = self._opennc(self.maskf)
        ncvar = nc.variables[self.maskvar]        
        if ncvar.ndim == 4: # If time dim present in mask field
            mask = ncvar[0,:,self.j1:self.j2+1, self.i1:self.i2+1] == self.maskmdi
        else: 
            mask = ncvar[:,self.j1:self.j2+1, self.i1:self.i2+1] == self.maskmdi
        nc.close()
        return mask 
 
    def _infill_coords(self):
        """ Infill missing coordinate info assuming cells are equally spaced """
        # First mask coordinsates using data mask
        if self.surface_field:
            coordmask = self.mask[0,:,:]
        else:
            coordmask = self.mask[0,0,:,:]
        self.x = np.ma.MaskedArray(self.x, mask=coordmask)    
        self.y = np.ma.MaskedArray(self.y, mask=coordmask) 
        
        # Linearly interpolate coordinates assuming equal spacing.
        fillx = np.array(self.x) 
        filly = np.array(self.y)
        idx=np.arange(len(fillx[0]))
        for n in range(fillx.shape[0]):
            missing = self.x.mask[n]
            not_missing = missing == False
            fillx[n,missing] = np.interp(idx[missing], idx[not_missing],self.x[n,not_missing])
            filly[n,missing] = np.interp(idx[missing], idx[not_missing],self.y[n,not_missing])
        self.x = np.array(fillx)
        self.y = np.array(filly)
        

    def _read_data(self):
        """ Read data from netcdf file(s) """
        nc = self._opennc(self.f)
        ncvar = nc.variables[self.var]
        if self.surface_field:
            self.data = ncvar[:, self.j1:self.j2+1, self.i1:self.i2+1]
            if self.data.ndim != 3:
                raise ShapeError('Surface fields must have dimensions [t,y,x]')
        else:
            self.data = ncvar[:, :, self.j1:self.j2+1, self.i1:self.i2+1]       
            if self.data.ndim != 4:
                raise ShapeError('Section fields must have dimensions [t,z,y,x]')
        nc.close()        

    def _read_xcoord(self):
        """ Read longitude data in decimal degrees from netcdf file(s). """
        nc = self._opennc(self.f)
        ncvar = nc.variables[self.xcoord]
        
        if ncvar.ndim == 2:
            self.x = ncvar[self.j1:self.j2+1,self.i1:self.i2+1]
            
        else:
            # Dummy j-index for averaging
            self.x = ncvar[self.i1:self.i2+1][np.newaxis] 
    
        nc.close()

        # Set range -180 to 180
        if self.x.max() > 180:
            self.x[self.x > 180] = self.x[self.x > 180] - 360
                
    def _read_ycoord(self):
        """ Read latitude data in decimal degrees from netcdf file(s) """
        nc = self._opennc(self.f)
        ncvar = nc.variables[self.ycoord]
        if ncvar.ndim == 2:
            self.y = ncvar[self.j1:self.j2+1,self.i1:self.i2+1]
        else:
             # Dummy j-index for averaging
            self.y = ncvar[self.j1:self.j2+1][:,np.newaxis]
        nc.close()
        
    def _read_tcoord(self):
        """ Read time coordinate information from netcdf file(s) """
        nc = self._opennc(self.f)
        t = nc.variables[self.tcoord]
        
        if len(glob.glob(self.f)) > 1:
            try:
                self.dates = num2date(MFTime(t)[:], calendar=t.calendar, units=t.units)
            except:
                print 'netcdf4.MFTime incompatible with NETCDF4. Try concatenating data into a single file.'
                raise NetCDF4ERROR(err)
        else:
            self.dates = num2date(t[:], calendar=t.calendar, units=t.units)            
            
    def _read_zcoord(self):
        """ Read z coordinate data from netcdf file(s) """
        if self.surface_field:
            self.z = None
        else:
            nc = self._opennc(self.f)
            self.z = np.abs(nc.variables[self.zcoord][:])
            nc.close()


def interpolate(s1, s2):
    """ Return data in s1 interpolated onto coordinates in s2 """
    sinterp = copy.deepcopy(s2)
    x0 = np.min([s1.xbounds[0],s2.xbounds[0]])
    y0 = np.min([s1.ybounds[0],s2.ybounds[0]])
    dist = s2.distance_along_section(x0,y0)
    dist_bounds = s2.distance_along_section(x0,y0, bounds=True)

    if s1.surface_field:
        bounds_mask=s2.bounds_mask[:,0,:]
        mask=s2.mask[:,0,:]
    else:
        bounds_mask=s2.bounds_mask
        mask=s2.mask
        
    sinterp.data = np.ma.MaskedArray(s1.interp_along_section(dist, x0, y0), mask=mask)    
    sinterp.bounds_data = np.ma.MaskedArray(s1.interp_along_section(dist_bounds, x0, y0), mask=bounds_mask)

    return sinterp
