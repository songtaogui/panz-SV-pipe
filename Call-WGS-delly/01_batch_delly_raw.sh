ref=$WORK/ref/Zea_mays.AGPv4.dna.toplevel.fa
exbed=$WORK/AMP_MEM_SV_TEs/exclude_bed/B73_gap100F50_cent.bed
for i in $WORK/AMP_MEM_SV_TEs/MarkDup_MEM/*.MEM.mkdup.bam
do
pre=`basename ${i%%.MEM.mkdup.bam}`
name=delly
cpu=1
echo "
#PBS -N ${name}_$pre
#PBS -l nodes=1:ppn=$cpu
#PBS -j oe
#PBS -q batch
#PBS -V

cd  \$PBS_O_WORKDIR

echo \"[\`date +\"%m-%d %H:%M\"\`] -----[\`du -sh|perl -pe 's/\s+\.//'\`]-----> ALL START\"

mkdir ${pre}_${name}

cd ${pre}_${name}

delly call -x $exbed -q 20 -r 20 -s 9 -n -g $ref -o ${pre}.bcf $i

echo \"[\`date +\"%m-%d %H:%M\"\`] -----[\`du -sh|perl -pe 's/\s+\.//'\`]-----> ALL DONE\"

" > $pre.pbs
done
