"""
Module containing code to work with Rapid observational data

"""

from netCDF4 import Dataset, num2date, date2num
import datetime
import numpy as np
import utils

class TransportObs(object):
    """ Template class to interface with observed ocean transports """

    def __init__(self, f, time_avg=None, mindt=None, maxdt=None):
        """ Create instance holding ocean transport data """
        self.f = f
        self.time_avg = time_avg
        self.mindt = mindt
        self.maxdt = maxdt
        self._read_data()

    def _read_data(self):
        """ Abstract method to read data and apply time averaging """
        pass
    
    def _read_dates(self):
        """ Abstract method to initialized dates """
        pass
    
    def _ym_dates(self):
        """ Return yearly mean date time objects """
        ym_dates = []
        
        for yr in range(self.yy.min(), self.yy.max()+1):
            ind = (self.yy == yr)
            if ind.any():
                ym_dates.append(datetime.datetime(yr, 7, 1))
                    
        return np.array(ym_dates)
    
    def _mm_dates(self):
        """ Return monthly mean date time objects """
        mm_dates = []
        
        for yr in range(self.yy.min(), self.yy.max()+1):
            for mon in range(1,12+1):
                ind = (self.yy == yr) & (self.mm == mon)
                if ind.any():
                    mm_dates.append(datetime.datetime(yr, mon, 15))
                    
        return np.array(mm_dates)
    
    def _calc_ym(self, data, profile=False):
        """ Return yearly mean values """
        ym_data = []
        
        for yr in range(self.yy.min(), self.yy.max()+1):
            ind = (self.yy == yr)
            if ind.any():
                if profile:
                    ym_data.append(np.mean(data[ind,:],axis=0))
                else:
                    ym_data.append(np.mean(data[ind]))
                    
        return np.array(ym_data)
    
    def _calc_mm(self, data, profile=False):
        """ Return monthly mean values """
        mm_data = []
        
        for yr in range(self.yy.min(), self.yy.max()+1):
            for mon in range(1,12+1):
                ind = (self.yy == yr) & (self.mm == mon)
                if ind.any():
                    if profile:
                        mm_data.append(np.mean(data[ind,:],axis=0))
                    else:
                        mm_data.append(np.mean(data[ind]))
        
        return np.array(mm_data)
         
    def _readnc(self, ncvar):
        """ Read variable from netcdf file """
        nc = Dataset(self.f)
        data = nc.variables[ncvar][:]
        nc.close()
        
        return data
    

class StreamFunctionObs(TransportObs):
    """ 
    Sub-class to hold overturning streamfunction observations
    from the RAPID-MOCHA-WBTS array at 26N.
    
    Data source:
    https://www.bodc.ac.uk/data/published_data_library/catalogue/
    10.5285/35784047-9b82-2160-e053-6c86abc0c91b/
    
    
    Data reference:
    Smeed D.; McCarthy G.; Rayner D.; Moat B.I.; Johns W.E.;
    Baringer M.O.; Meinen C.S. (2016). Atlantic meridional 
    overturning circulation observed by the RAPID-MOCHA-WBTS
    (RAPID-Meridional Overturning Circulation and Heatflux 
    Array-Western Boundary Time Series) array at 26N from 
    2004 to 2015. British Oceanographic Data Centre - Natural
    Environment Research Council, UK. doi:10/bkzc.
 
    """
    def _read_data(self):
        """ Read data and apply time averaging """
        self._read_dates()
        self.z = self._readnc('depth')
        
        if self.time_avg is None:
            self.dates = self.original_dates
            self.sf = self._readnc('stream_function_mar').transpose()
        elif self.time_avg == 'monthly':
            self.dates = self._mm_dates()
            self.sf = self._calc_mm(self._readnc('stream_function_mar').transpose(),
                                     profile=True)
        elif self.time_avg == 'yearly':
            self.dates = self._ym_dates()
            self.sf = self._calc_ym(self._readnc('stream_function_mar').transpose(),
                                     profile=True)
        else:
            print self.time_avg
            raise ValueError('time_avg must be "monthly" or "yearly"')
            
        if (self.mindt is not None) and (self.maxdt is not None):
            tind = utils.get_dateind(self.dates, self.mindt, self.maxdt)
            self.sf = self.sf[tind,:]
            self.dates = self.dates[tind]

    def _read_dates(self):
        """ Read date information from file """
        nc = Dataset(self.f)
        t = nc.variables['time']
        self.original_dates = num2date(t[:],units=t.units)
        self.hh = np.array([dt.hour for dt in self.original_dates], dtype=np.int)
        self.dd = np.array([dt.day for dt in self.original_dates], dtype=np.int)
        self.mm = np.array([dt.month for dt in self.original_dates], dtype=np.int)
        self.yy = np.array([dt.year for dt in self.original_dates], dtype=np.int)

    def write_to_netcdf(self, ncfile):
        """ Write observation data to netcdf file """
        
        # Open ncfile and create coords
        dataset = Dataset(ncfile, 'w', format='NETCDF4_CLASSIC')
        zdim = dataset.createDimension('depth', self.z.size)
        tdim = dataset.createDimension('time', None)
        
        # Create time coordinate
        time = dataset.createVariable('time',np.float64,(tdim._name,))
        time.units = 'hours since 0001-01-01 00:00:00.0'
        time.calendar = 'gregorian'
        time[:] = date2num(self.dates, time.units, calendar=time.calendar)
        
        # Create depth coordinate 
        z = dataset.createVariable('depth',np.float64,(zdim._name,))
        z.units = 'm'
        z[:] = self.z

        # Create streamfunction variable
        sf = dataset.createVariable('stream_function_mar',np.float64,(tdim._name, zdim._name))
        sf.units = 'Sv'
        sf[:] = self.sf

        # Close file
        print 'SAVING: %s' % ncfile
        dataset.close()
        


class VolumeTransportObs(TransportObs):
    """ 
    Sub-class to hold volume transport observations
    from the RAPID-MOCHA-WBTS array at 26N.
    
    Data source:
    https://www.bodc.ac.uk/data/published_data_library/catalogue/
    10.5285/35784047-9b82-2160-e053-6c86abc0c91b/
    
    
    Data reference:
    Smeed D.; McCarthy G.; Rayner D.; Moat B.I.; Johns W.E.;
    Baringer M.O.; Meinen C.S. (2016). Atlantic meridional 
    overturning circulation observed by the RAPID-MOCHA-WBTS
    (RAPID-Meridional Overturning Circulation and Heatflux 
    Array-Western Boundary Time Series) array at 26N from 
    2004 to 2015. British Oceanographic Data Centre - Natural
    Environment Research Council, UK. doi:10/bkzc.
 
    """
     
    def _read_data(self):
        """ Read data and apply time averaging """
        self._read_dates()
        
        if self.time_avg is None:
            self.dates = self.original_dates
            self.ekman = self._readnc('t_ek10')
            self.umo = self._readnc('t_umo10')
            self.fc = self._readnc('t_gs10')
            self.moc = self._readnc('moc_mar_hc10')           
        elif self.time_avg == 'monthly':
            self.dates = self._mm_dates()
            self.ekman = self._calc_mm(self._readnc('t_ek10'))
            self.umo = self._calc_mm(self._readnc('t_umo10'))
            self.fc = self._calc_mm(self._readnc('t_gs10'))
            self.moc = self._calc_mm(self._readnc('moc_mar_hc10'))
        elif self.time_avg == 'yearly':
            self.dates = self._ym_dates()
            self.ekman = self._calc_ym(self._readnc('t_ek10'))
            self.umo = self._calc_ym(self._readnc('t_umo10'))
            self.fc = self._calc_ym(self._readnc('t_gs10'))
            self.moc = self._calc_ym(self._readnc('moc_mar_hc10'))
        else:
            print self.time_avg
            raise ValueError('time_avg must be "monthly" or "yearly"')
        
        if (self.mindt is not None) and (self.maxdt is not None):
            tind = utils.get_dateind(self.dates, self.mindt, self.maxdt)
            self.ekman = self.ekman[tind]
            self.umo = self.umo[tind]
            self.fc = self.fc[tind]
            self.moc = self.moc[tind]
            self.dates = self.dates[tind]

    def _read_dates(self):
        """ Read date information from file """
        nc = Dataset(self.f)
        t = nc.variables['time']
        self.original_dates = num2date(t[:],units=t.units)
        self.hh = np.array([dt.hour for dt in self.original_dates], dtype=np.int)
        self.dd = np.array([dt.day for dt in self.original_dates], dtype=np.int)
        self.mm = np.array([dt.month for dt in self.original_dates], dtype=np.int)
        self.yy = np.array([dt.year for dt in self.original_dates], dtype=np.int)

    def write_to_netcdf(self, ncfile):
        """ Write observation data to netcdf file """
        
        # Open ncfile and create coords
        dataset = Dataset(ncfile, 'w', format='NETCDF4_CLASSIC')
        tdim = dataset.createDimension('time', None)
        
        # Create time coordinate
        time = dataset.createVariable('time',np.float64,(tdim._name,))
        time.units = 'hours since 0001-01-01 00:00:00.0'
        time.calendar = 'gregorian'
        time[:] = date2num(self.dates, time.units, calendar=time.calendar)

        # Create variables
        ek = dataset.createVariable('t_ek10',np.float64,(tdim._name,))
        ek.units = 'Sv'
        ek[:] = self.ekman

        umo = dataset.createVariable('t_umo10',np.float64,(tdim._name,))
        umo.units = 'Sv'
        umo[:] = self.umo

        fc = dataset.createVariable('t_gs10',np.float64,(tdim._name,))
        fc.units = 'Sv'
        fc[:] = self.fc

        moc = dataset.createVariable('t_moc_mar_hc10',np.float64,(tdim._name,))
        moc.units = 'Sv'
        moc[:] = self.moc
        
        # Close file
        print 'SAVING: %s' % ncfile
        dataset.close()



class HeatTransportObs(TransportObs):
    """ 
    Sub-class to hold meridional heat transport observations
    from the RAPID-MOCHA-WBTS array at 26N.
    
    Data source:
    https://www.rsmas.miami.edu/users/mocha/mocha_results.htm 
    
    Data reference:
    http://journals.ametsoc.org/doi/abs/10.1175/2010JCLI3997.1
    
    """
    
    def _read_data(self):
        """ 
        Read data at original frequency or calculate a time-average
        
        """
        self._read_dates()
        self.z = self._readnc('z')
                
        if self.time_avg is None:
            self.dates = self.original_dates
            self.q_eddy = self._readnc('Q_eddy') / 1e15
            self.q_ek = self._readnc('Q_ek') / 1e15
            self.q_fc = self._readnc('Q_fc') / 1e15
            self.q_gyre = self._readnc('Q_gyre') / 1e15
            self.q_geoint = self._readnc('Q_int') / 1e15
            self.q_mo = self._readnc('Q_mo') / 1e15
            self.q_ot = self._readnc('Q_ot') / 1e15
            self.q_sum = self._readnc('Q_sum') / 1e15
            self.q_wbw = self._readnc('Q_wedge') / 1e15
            self.t_basin = self._readnc('T_basin')
            self.v_basin = self._readnc('V_basin')
            self.v_fc = self._readnc('V_fc')
        elif self.time_avg == 'monthly':
            self.dates = self._mm_dates()
            self.q_eddy = self._calc_mm(self._readnc('Q_eddy')) / 1e15
            self.q_ek = self._calc_mm(self._readnc('Q_ek')) / 1e15
            self.q_fc = self._calc_mm(self._readnc('Q_fc')) / 1e15
            self.q_gyre = self._calc_mm(self._readnc('Q_gyre')) / 1e15
            self.q_geoint = self._calc_mm(self._readnc('Q_int')) / 1e15
            self.q_mo = self._calc_mm(self._readnc('Q_mo')) / 1e15
            self.q_ot = self._calc_mm(self._readnc('Q_ot')) / 1e15
            self.q_sum = self._calc_mm(self._readnc('Q_sum')) / 1e15
            self.q_wbw = self._calc_mm(self._readnc('Q_wedge')) / 1e15
            self.t_basin = self._calc_mm(self._readnc('T_basin'), profile=True)
            self.v_basin = self._calc_mm(self._readnc('V_basin'), profile=True)
            self.v_fc = self._calc_mm(self._readnc('V_fc'), profile=True)   
        elif self.time_avg == 'yearly':
            self.dates = self._ym_dates()
            self.q_eddy = self._calc_ym(self._readnc('Q_eddy')) / 1e15
            self.q_ek = self._calc_ym(self._readnc('Q_ek')) / 1e15
            self.q_fc = self._calc_ym(self._readnc('Q_fc')) / 1e15
            self.q_gyre = self._calc_ym(self._readnc('Q_gyre')) / 1e15
            self.q_geoint = self._calc_ym(self._readnc('Q_int')) / 1e15
            self.q_mo = self._calc_ym(self._readnc('Q_mo')) / 1e15
            self.q_ot = self._calc_ym(self._readnc('Q_ot')) / 1e15
            self.q_sum = self._calc_ym(self._readnc('Q_sum')) / 1e15
            self.q_wbw = self._calc_ym(self._readnc('Q_wedge')) / 1e15
            self.t_basin = self._calc_ym(self._readnc('T_basin'), profile=True)
            self.v_basin = self._calc_ym(self._readnc('V_basin'), profile=True)
            self.v_fc = self._calc_ym(self._readnc('V_fc'), profile=True)
        else:
            print self.time_avg
            raise ValueError('time_avg must be "monthly" or "yearly"')
 
        if (self.mindt is not None) and (self.maxdt is not None):
            tind = utils.get_dateind(self.dates, self.mindt, self.maxdt)
            self.q_eddy = self.q_eddy[tind]
            self.q_ek = self.q_ek[tind]
            self.q_fc = self.q_fc[tind]
            self.q_gyre = self.q_gyre[tind]
            self.q_geoint = self.q_geoint[tind]
            self.q_mo = self.q_mo[tind]
            self.q_ot = self.q_ot[tind]
            self.q_sum = self.q_sum[tind]
            self.q_wbw = self.q_wbw[tind]
            self.t_basin = self.t_basin[tind,:] 
            self.v_basin = self.v_basin[tind,:]
            self.v_fc = self.v_fc[tind,:]
            self.dates = self.dates[tind]

    def _read_dates(self):
        """ Read date information from file """
        dts = []
        self.hh = np.array(self._readnc('hour'), dtype=np.int)
        self.dd = np.array(self._readnc('day'), dtype=np.int)
        self.mm = np.array(self._readnc('month'), dtype=np.int)
        self.yy = np.array(self._readnc('year'), dtype=np.int)

        for ndt in xrange(len(self.hh)):
            dts.append(datetime.datetime(
                self.yy[ndt], self.mm[ndt], self.dd[ndt], self.hh[ndt],0,0))
        
        self.original_dates = np.array(dts)
    
    def write_to_netcdf(self, ncfile):
        """ Write observation data to netcdf file """
        
        # Open ncfile and create coords
        dataset = Dataset(ncfile, 'w', format='NETCDF4_CLASSIC')
        tdim = dataset.createDimension('time', None)
        zdim = dataset.createDimension('depth', self.z.size)

        # Create time coordinate
        time = dataset.createVariable('time',np.float64,(tdim._name,))
        time.units = 'hours since 0001-01-01 00:00:00.0'
        time.calendar = 'gregorian'
        time[:] = date2num(self.dates, time.units, calendar=time.calendar)
        
        # Create depth coordinate 
        z = dataset.createVariable('depth',np.float64,(zdim._name,))
        z.units = 'm'
        z[:] = self.z

        # Create variables
        q_eddy = dataset.createVariable('Q_eddy',np.float64,(tdim._name,))
        q_eddy.units = 'PW'
        q_eddy[:] = self.q_eddy

        q_ek = dataset.createVariable('Q_ek',np.float64,(tdim._name,))
        q_ek.units = 'PW'
        q_ek[:] = self.q_ek

        q_fc = dataset.createVariable('Q_fc',np.float64,(tdim._name,))
        q_fc.units = 'PW'
        q_fc[:] = self.q_fc

        q_gyre = dataset.createVariable('Q_gyre',np.float64,(tdim._name,))
        q_gyre.units = 'PW'
        q_gyre[:] = self.q_gyre

        q_geoint = dataset.createVariable('Q_int',np.float64,(tdim._name,))
        q_geoint.units = 'PW'
        q_geoint[:] = self.q_geoint

        q_mo = dataset.createVariable('Q_mo',np.float64,(tdim._name,))
        q_mo.units = 'PW'
        q_mo[:] = self.q_mo

        q_ot = dataset.createVariable('Q_ot',np.float64,(tdim._name,))
        q_ot.units = 'PW'
        q_ot[:] = self.q_ot

        q_sum = dataset.createVariable('Q_sum',np.float64,(tdim._name,))
        q_sum.units = 'PW'
        q_sum[:] = self.q_sum

        q_wbw = dataset.createVariable('Q_wedge',np.float64,(tdim._name,))
        q_wbw.units = 'PW'
        q_wbw[:] = self.q_wbw

        t_basin = dataset.createVariable('T_basin',np.float64,(tdim._name,zdim._name,))
        t_basin.units = 'degC'
        t_basin[:] = self.t_basin

        v_basin = dataset.createVariable('V_basin',np.float64,(tdim._name,zdim._name,))
        v_basin.units = 'Sv/m'
        v_basin[:] = self.v_basin

        v_fc = dataset.createVariable('V_fc',np.float64,(tdim._name,zdim._name,))
        v_fc.units = 'Sv/m'
        v_fc[:] = self.v_fc

        # Close file
        print 'SAVING: %s' % ncfile
        dataset.close()

                    

class FloridaCurrentObs(TransportObs):
    """ 
    Class to hold Florida current transport estimates derived from
    submarine cable measurements.
        
    Data source:
    http://www.aoml.noaa.gov/phod/floridacurrent/data_access.php
    
    The Florida Current cable and section data are made freely available
    on the Atlantic Oceanographic and Meteorological Laboratory web page 
    (www.aoml.noaa.gov/phod/floridacurrent/) and are funded by the DOC-NOAA
    Climate Program Office - Ocean Observing and Monitoring Division.

    The project scientists would also appreciate it if you informed us of 
    any publications or presentations that you prepare using this data.
    Continued funding of this project depends on us being able to justify 
    to NOAA (and hence the US Congress) the usefulness of this data.
        
    """
    
    def _read_data(self):
        """ Read data and apply time averaging """
        self._read_dates()
        
        if self.time_avg is None:
            self.fc = self._readnc('florida_current_transport')
        elif self.time_avg == 'monthly':
            self.dates = self._mm_dates()
            self.fc = self._calc_mm(self._readnc('florida_current_transport'))
        elif self.time_avg == 'yearly':
            self.dates = self._ym_dates()
            self.fc = self._calc_ym(self._readnc('florida_current_transport'))
        else:
            print self.time_avg
            raise ValueError('time_avg must be "monthly" or "yearly"')
  
        if (self.mindt is not None) and (self.maxdt is not None):
            tind = utils.get_dateind(self.dates, self.mindt, self.maxdt)
            self.fc = self.fc[tind]
            self.dates = self.dates[tind]

    def _read_dates(self):
        """ Read date information from file """
        nc = Dataset(self.f)
        t = nc.variables['time']
        self.original_dates = num2date(t[:],units=t.units)
        self.hh = np.array([dt.hour for dt in self.original_dates], dtype=np.int)
        self.dd = np.array([dt.day for dt in self.original_dates], dtype=np.int)
        self.mm = np.array([dt.month for dt in self.original_dates], dtype=np.int)
        self.yy = np.array([dt.year for dt in self.original_dates], dtype=np.int)


    def write_to_netcdf(self, ncfile):
        """ Write observation data to netcdf file """
        
        # Open ncfile and create coords
        dataset = Dataset(ncfile, 'w', format='NETCDF4_CLASSIC')
        tdim = dataset.createDimension('time', None)

        # Create time coordinate
        time = dataset.createVariable('time',np.float64,(tdim._name,))
        time.units = 'hours since 0001-01-01 00:00:00.0'
        time.calendar = 'gregorian'
        time[:] = date2num(self.dates, time.units, calendar=time.calendar)
 
        # Create variables
        fc = dataset.createVariable('florida_current_transport',np.float64,(tdim._name,))
        fc.units = 'Sv'
        fc[:] = self.fc

        # Close file
        print 'SAVING: %s' % ncfile
        dataset.close()
