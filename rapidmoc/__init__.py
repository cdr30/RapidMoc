"""
RapidMoc is a package to calculate diagnostics of the Atlantic Meridional 
Overturning Circulation (AMOC) using output from an ocean general circulation
model for comparison with observed data from the RAPID-MOCHA-WTBS array
at 26N.

Author:
Chris Roberts      ECMWF     chris.roberts@ecmwf.int

testing command:

# ORCA1 

run_rapidmoc.py config_orca1.ini "/hugetmp/necr/dat/gp69-0/ORCA1_Z75/ann/gp69_1y_197?_grid_T.nc" "/hugetmp/necr/dat/gp69-0/ORCA1_Z75/ann/gp69_1y_197?_grid_T.nc" "/hugetmp/necr/dat/gp69-0/ORCA1_Z75/ann/gp69_1y_197?_grid_U.nc" "/hugetmp/necr/dat/gp69-0/ORCA1_Z75/ann/gp69_1y_197?_grid_V.nc"


run_rapidmoc.py etc/config.ini.orca025 "/hugetmp/necr/dat/gp6b-0/ORCA025_Z75/ann/gp6b_1y_197?_grid_T.nc" "/hugetmp/necr/dat/gp6b-0/ORCA025_Z75/ann/gp6b_1y_197?_grid_T.nc" "/hugetmp/necr/dat/gp6b-0/ORCA025_Z75/ann/gp6b_1y_197?_grid_U.nc" "/hugetmp/necr/dat/gp6b-0/ORCA025_Z75/ann/gp6b_1y_197?_grid_V.nc"

"""
