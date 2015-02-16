#!/usr/bin/bash
# Download millenium simulation data for different snapnums
for i in `seq 0 63`;
do
    python get_data.py datafile_$i.txt --snapnum $i
done
