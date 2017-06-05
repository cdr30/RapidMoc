"""
RapidMoc is a package to calculate diagnostics of the Atlantic Meridional 
Overturning Circulation (AMOC) using output from an ocean general circulation
model for comparison with observed data from the RAPID-MOCHA-WTBS array
at 26N.

Author:
Chris Roberts      ECMWF     chris.roberts@ecmwf.int

testing command:

# ORCA1 

run_rapidmoc.py etc/config.ini.orca1 "/hugetmp/necr/dat/gp69-0/ORCA1_Z75/ann/gp69_1y_197?_grid_T.nc" "/hugetmp/necr/dat/gp69-0/ORCA1_Z75/ann/gp69_1y_197?_grid_T.nc" "/hugetmp/necr/dat/gp69-0/ORCA1_Z75/ann/gp69_1y_197?_grid_U.nc" "/hugetmp/necr/dat/gp69-0/ORCA1_Z75/ann/gp69_1y_197?_grid_V.nc"

TFILE=/hugetmp/necr/ecflow/ecmwf_ifs_lr_spinup_1950_cmor2jas_gp69_0_0/cmorized_data/gp69/0/19500101/1950/thetao_Omon_ECMWF-IFS-LR_spinup-1950_r1i1p1f1_gn_195001-195012.nc
SFILE=/hugetmp/necr/ecflow/ecmwf_ifs_lr_spinup_1950_cmor2jas_gp69_0_0/cmorized_data/gp69/0/19500101/1950/so_Omon_ECMWF-IFS-LR_spinup-1950_r1i1p1f1_gn_195001-195012.nc
UFILE=/hugetmp/necr/ecflow/ecmwf_ifs_lr_spinup_1950_cmor2jas_gp69_0_0/cmorized_data/gp69/0/19500101/1950/tauuo_Omon_ECMWF-IFS-LR_spinup-1950_r1i1p1f1_gn_195001-195012.nc
VFILE=/hugetmp/necr/ecflow/ecmwf_ifs_lr_spinup_1950_cmor2jas_gp69_0_0/cmorized_data/gp69/0/19500101/1950/vo_Omon_ECMWF-IFS-LR_spinup-1950_r1i1p1f1_gn_195001-195012.nc

run_rapidmoc.py etc/config.ini.orca1_cmorized $TFILE $SFILE $UFILE $VFILE
run_rapidmoc.py etc/config.ini.orca025 "/hugetmp/necr/dat/gp6b-0/ORCA025_Z75/ann/gp6b_1y_197?_grid_T.nc" "/hugetmp/necr/dat/gp6b-0/ORCA025_Z75/ann/gp6b_1y_197?_grid_T.nc" "/hugetmp/necr/dat/gp6b-0/ORCA025_Z75/ann/gp6b_1y_197?_grid_U.nc" "/hugetmp/necr/dat/gp6b-0/ORCA025_Z75/ann/gp6b_1y_197?_grid_V.nc"

"""
