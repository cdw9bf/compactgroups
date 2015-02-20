#!/usr/bin/bash
# Find compact groups at many snapnums
for i in `seq 0 63`;
do
    echo Working on snapnum $i
    python FindCompactGroups.py datafile_$i.txt --groupfile groups_$i.txt --memberfile members_$i.txt
    echo
    echo
    echo
    echo
done
