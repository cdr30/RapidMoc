#!/bin/ksh

module load nco

# Define help function
display_help() {
    echo "  Retrieve Florida Current data from NOAA website  "
    echo ""
    echo "  USAGE: $(basename $0) MINYR MAXYR"
    echo ""
    echo "  MINYR - First year of data to retrieve"
    echo "  MAXYR - last year of data to retrieve"
}

# Test number of arguments
if [[ $# -lt 2 ]]; then
  display_help
  exit 1
fi

# Get arguments
MINYR=$1
MAXYR=$2


NCFILES=""
# Loop through years and retrieve from web
for YEAR in $(seq $MINYR $MAXYR); do
  wget http://www.aoml.noaa.gov/phod/floridacurrent/FC_cable_transport_${YEAR}.dat
  
  if [[ $? = 0 ]]; then
    python3 dat2nc.py FC_cable_transport_${YEAR}.dat
    NCFILES="$NCFILES FC_cable_transport_${YEAR}.dat.nc"
  else
    echo "ERROR: could not retrieve http://www.aoml.noaa.gov/phod/floridacurrent/FC_cable_transport_${YEAR}.dat"
  fi
done

# Concatenate all files into single netcdf
OUTFILE=FC_cable_transport_${MINYR}-${MAXYR}.nc
echo "SAVING: $OUTFILE"
ncrcat $NCFILES $OUTFILE
rm *.dat.nc

exit 0


