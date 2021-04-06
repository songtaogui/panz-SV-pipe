ref=$WORK/ref/Zea_mays.AGPv4.dna.toplevel.fa
exbed=$WORK/AMP_MEM_SV_TEs/exclude_bed/B73_gap100F50_cent.bed
merge_bcf="$WORK/AMP_MEM_SV_TEs/delly/MKdup/delly_merge/ALL_delly_merge.bcf"
for i in $WORK/AMP_MEM_SV_TEs/MarkDup_MEM/*.MEM.mkdup.bam
do
pre=`basename ${i%%.*}`
name=dgeno
cpu=2
echo "
#PBS -N gst_${name}_$pre
#PBS -l nodes=1:ppn=$cpu
#PBS -j oe
#PBS -q batch
#PBS -V
source \$GST/.zshrc
cd  \$PBS_O_WORKDIR

echo \"[\`date +\"%m-%d %H:%M\"\`] ----------> ALL START\"

mkdir ${pre}_${name}

cd ${pre}_${name}

delly call -x $exbed -g $ref -v $merge_bcf -o ${pre}_geno.bcf $i

bcftools convert -O z -o $WORK/AMP_MEM_SV_TEs/delly/MKdup/delly_geno_merge/${pre}_delly_geno.vcf.gz ${pre}_geno.bcf

bcftool index -t $WORK/AMP_MEM_SV_TEs/delly/MKdup/delly_geno_merge/${pre}_delly_geno.vcf.gz

echo \"[\`date +\"%m-%d %H:%M\"\`] ----------> ALL DONE\"

" > ${pre}_geno.pbs
done
