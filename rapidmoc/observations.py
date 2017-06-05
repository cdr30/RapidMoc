"""
Module containing code to work with Rapid observational data

"""

from netCDF4 import Dataset, num2date
import datetime
import numpy as np
import abc


class TransportObs(object):
    """ Template class to interface with observed ocean transports """
    
    def __init__(self, f, time_avg=None, mindt=None, maxdt=None):
        """ Create instance holding ocean transport data """
        self.f = f
        self.time_avg = time_avg
        self.mindt = mindt
        self.maxdt = maxdt
        self._read_data()

    @abc.abstractmethod
    def _read_data(self):
        """ Abstract method to read data and apply time averaging """
        pass
    
    @abc.abstractmethod
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
                    ym_data.append(np.mean(data[:,ind],axis=1))
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
                        mm_data.append(np.mean(data[:,ind],axis=1))
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
            self.streamfunction = self._readnc('stream_function_mar').transpose()
        elif self.time_avg == 'monthly':
            self.dates = self._mm_dates()
            self.streamfunction = self._calc_mm(self._readnc('stream_function_mar'),
                                     profile=True)
        elif self.time_avg == 'yearly':
            self.dates = self._ym_dates()
            self.streamfunction = self._calc_ym(self._readnc('stream_function_mar'),
                                     profile=True)
        else:
            print self.time_avg
            raise ValueError('time_avg must be "monthly" or "yearly"')
 
    def _read_dates(self):
        """ Read date information from file """
        nc = Dataset(self.f)
        t = nc.variables['time']
        self.original_dates = num2date(t[:],units=t.units)
        self.hh = np.array([dt.hour for dt in self.original_dates], dtype=np.int)
        self.dd = np.array([dt.day for dt in self.original_dates], dtype=np.int)
        self.mm = np.array([dt.month for dt in self.original_dates], dtype=np.int)
        self.yy = np.array([dt.year for dt in self.original_dates], dtype=np.int)

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
            self.fs = self._readnc('t_gs10')
            self.moc = self._readnc('moc_mar_hc10')           
        elif self.time_avg == 'monthly':
            self.dates = self._mm_dates()
            self.ekman = self._calc_mm(self._readnc('t_ek10'))
            self.umo = self._calc_mm(self._readnc('t_umo10'))
            self.fs = self._calc_mm(self._readnc('t_gs10'))
            self.moc = self._calc_mm(self._readnc('moc_mar_hc10'))
        elif self.time_avg == 'yearly':
            self.dates = self._ym_dates()
            self.ekman = self._calc_ym(self._readnc('t_ek10'))
            self.umo = self._calc_ym(self._readnc('t_umo10'))
            self.fs = self._calc_ym(self._readnc('t_gs10'))
            self.moc = self._calc_ym(self._readnc('moc_mar_hc10'))
        else:
            print self.time_avg
            raise ValueError('time_avg must be "monthly" or "yearly"')
 
    def _read_dates(self):
        """ Read date information from file """
        nc = Dataset(self.f)
        t = nc.variables['time']
        self.original_dates = num2date(t[:],units=t.units)
        self.hh = np.array([dt.hour for dt in self.original_dates], dtype=np.int)
        self.dd = np.array([dt.day for dt in self.original_dates], dtype=np.int)
        self.mm = np.array([dt.month for dt in self.original_dates], dtype=np.int)
        self.yy = np.array([dt.year for dt in self.original_dates], dtype=np.int)


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
        Read data at original frequency or convert to monthly means
        
        """
        self._read_dates()
        
        if self.time_avg is None:
            self.dates = self.original_dates
            self.q_eddy = self._readnc('Q_eddy')
            self.q_ek = self._readnc('Q_ek')
            self.q_fc = self._readnc('Q_fc')
            self.q_gyre = self._readnc('Q_gyre')
            self.q_int = self._readnc('Q_int')
            self.q_mo = self._readnc('Q_mo')
            self.q_ot = self._readnc('Q_ot')
            self.q_sum = self._readnc('Q_sum')
            self.q_wedge = self._readnc('Q_wedge')
        elif self.time_avg == 'monthly':
            self.dates = self._mm_dates()
            self.q_eddy = self._calc_mm(self._readnc('Q_eddy'))
            self.q_ek = self._calc_mm(self._readnc('Q_ek'))
            self.q_fc = self._calc_mm(self._readnc('Q_fc'))
            self.q_gyre = self._calc_mm(self._readnc('Q_gyre'))
            self.q_int = self._calc_mm(self._readnc('Q_int'))
            self.q_mo = self._calc_mm(self._readnc('Q_mo'))
            self.q_ot = self._calc_mm(self._readnc('Q_ot'))
            self.q_sum = self._calc_mm(self._readnc('Q_sum'))
            self.q_wedge = self._calc_mm(self._readnc('Q_wedge'))
        elif self.time_avg == 'yearly':
            self.dates = self._ym_dates()
            self.q_eddy = self._calc_ym(self._readnc('Q_eddy'))
            self.q_ek = self._calc_ym(self._readnc('Q_ek'))
            self.q_fc = self._calc_ym(self._readnc('Q_fc'))
            self.q_gyre = self._calc_ym(self._readnc('Q_gyre'))
            self.q_int = self._calc_ym(self._readnc('Q_int'))
            self.q_mo = self._calc_ym(self._readnc('Q_mo'))
            self.q_ot = self._calc_ym(self._readnc('Q_ot'))
            self.q_sum = self._calc_ym(self._readnc('Q_sum'))
            self.q_wedge = self._calc_ym(self._readnc('Q_wedge'))
        else:
            print self.time_avg
            raise ValueError('time_avg must be "monthly" or "yearly"')
 
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

                    
