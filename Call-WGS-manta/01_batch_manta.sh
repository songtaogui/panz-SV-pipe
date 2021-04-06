ref=$WORK/ref/Zea_mays.AGPv4.dna.toplevel.fa
#include bed should be bgzipped and tab indexed 
inbed=$WORK/AMP_MEM_SV_TEs/exclude_bed/B73_Include.bed.gz
for i in $WORK/AMP_MEM_SV_TEs/MarkDup_MEM/*.MEM.mkdup.bam
do
pre=`basename ${i%%.*}`
name=manta
cpu=2
echo "
#PBS -N gst_${name}_$pre
#PBS -l nodes=1:ppn=$cpu
#PBS -j oe
#PBS -q batch
#PBS -V

cd  \$PBS_O_WORKDIR

echo \"[\`date +\"%m-%d %H:%M\"\`] ---------> ALL START\"

mkdir ${pre}_${name}

cd ${pre}_${name}

configManta.py --bam $i --referenceFasta $ref --callRegions $inbed --runDir ./ --generateEvidenceBam --outputContig

./runWorkflow.py -m local -j $cpu

rm -rf ./workspace

qstat -f \$PBS_JOBID|grep \"used\"

echo \"[\`date +\"%m-%d %H:%M\"\`] ----------> ALL DONE\"

" > $pre.pbs
done
