#!/usr/bin/ksh

for f in $(ls test_data_199901-199912_*.png); do
  eog $f expected_output/$f
done
