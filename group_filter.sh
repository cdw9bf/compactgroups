#!/usr/bin/bash
# Find compact groups at many snapnums
# COMMANDS: % bash group_filter.sh CGDATADIRECTORY

USAGE="Usage:$0 cgDataDir"
if [ "$#" == "0" ]; then
    echo "$USAGE"
    exit 1
fi
dir=$1

tmpmember="tmp_members.txt"
tmpgroup="tmp_groups.txt"
if [[ -f $tmpgroup ]]; then
    rm $tmpgroup
fi

# velocity cut loop
velfile="vel_filter_snapnum.txt"
if [[ -f $velfile ]]; then
    rm $velfile
fi

echo "snapnum,vel_cut,num_gp" >> $velfile
for i in `seq 0 63`; do
    for j in `seq 1 10`; do
        echo Working on snapnum $i
        vel=`echo "$j*100" | bc`
        python FindCompactGroups.py $dir"datafile_$i.txt" --groupfile $tmpgroup --memberfile $tmpmember --velocity_filter $vel
        num_gp=`wc -l < $tmpgroup`
        num_gp=`echo $num_gp-1 | bc`
        echo "$i,$vel,$num_gp" >> $velfile
        echo
        echo
        echo
        echo
        rm $tmpgroup
    done
done

# band width loop
bwfile="bandwidth_snapnum.txt"
if [[ -f $bwfile ]]; then
    rm $bwfile
fi

echo "snapnum,bandwidth,num_gp" >> $bwfile
for i in `seq 0 63`; do
    for j in `seq 1 6`; do
        echo Working on snapnum $i
        bw=`echo "$j*0.05" | bc`
        python FindCompactGroups.py $dir"datafile_$i.txt" --groupfile $tmpgroup --memberfile $tmpmember --bandwidth $bw
        num_gp=`wc -l < $tmpgroup`
        num_gp=`echo $num_gp-1 | bc`
        echo "$i,$bw,$num_gp" >> $bwfile
        echo
        echo
        echo
        echo
        rm $tmpgroup
    done
done

# separation ratio loop
srfile="sep_ratio_snapnum.txt"
if [[ -f $srfile ]]; then
    rm $srfile
fi

echo "snapnum,sep_ratio,num_gp" >> $srfile
for i in `seq 0 63`; do
    for j in `seq 3 8`; do
        echo Working on snapnum $i
        sr=`echo "$j*0.2" | bc`
        python FindCompactGroups.py $dir"datafile_$i.txt" --groupfile $tmpgroup --memberfile $tmpmember --max_sep_ratio $sr
        num_gp=`wc -l < $tmpgroup`
        num_gp=`echo $num_gp-1 | bc`
        echo "$i,$sr,$num_gp" >> $srfile
        echo
        echo
        echo
        echo
        rm $tmpgroup
    done
done

rm $tmpmember

echo "JOBS DONE!"
