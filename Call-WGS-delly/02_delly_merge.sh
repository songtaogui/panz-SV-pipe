#!/bin/bash
#PBS -N delly_merge
#PBS -j oe
#PBS -V
#PBS -l nodes=1:ppn=4
#PBS -q batch

cd $PBS_O_WORKDIR

echo "[`date +"%m-%d %H:%M"`] -----[`du -sh|perl -pe 's/\s+\.//'`]-----> ALL START"

delly merge -o ALL_delly_merge.bcf -c -p -b 1000 -r 0.8 -n 1000000 $WORK/AMP_MEM_SV_TEs/delly/MKdup/delly_raw/Q*_delly/*.bcf $WORK/AMP_MEM_SV_TEs/delly/MKdup/delly_raw_130/*.bcf 1>delly_merge.log 2>&1

qstat -f $PBS_JOBID|grep "used"
echo "[`date +"%m-%d %H:%M"`] -----[`du -sh|perl -pe 's/\s+\.//'`]-----> ALL END"


