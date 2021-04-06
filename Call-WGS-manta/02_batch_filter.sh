ref=$WORK/ref/Zea_mays.AGPv4.dna.toplevel.fa
for i in $WORK/AMP_MEM_SV_TEs/manta/MKdup/manta_call/Q*_manta/results/variants/diploidSV.vcf.gz
do
	pre0=${i%%_manta/results/variants/diploidSV.vcf.gz}
	pre=`basename $pre0`
	name=mf
	cpu=1
	echo "
#PBS -N gst_${name}_$pre
#PBS -l nodes=1:ppn=$cpu
#PBS -j oe
#PBS -q batch
#PBS -V
source \$GST/.zshrc
cd  \$PBS_O_WORKDIR

echo \"[\`date +\"%m-%d %H:%M\"\`] ---------> ALL START\"

zcat $i | perl -lane '/^#/ && print && next; print if \$F[6] eq \"PASS\";' | perl -lane '/^#/ && print && next; print unless \$F[7]=~/IMPRECISE/;' | perl -lane '/^#/ && print && next; if( \$F[7]=~/SVLEN=(-*\d+);/ ){ \$len=abs(\$1);print if \$len >=50 && \$len <=1000000; }else{print;}' | gzip > ${pre}_manta_filtered.vcf.gz

echo \"[\`date +\"%m-%d %H:%M\"\`] ---------> DONE filter\"

bayesTyperTools convertAllele -v ${pre}_manta_filtered.vcf.gz -g $ref -o ${pre}_manta_bystypr.vcf.gz -z

echo \"[\`date +\"%m-%d %H:%M\"\`] ---------> ALL DONE\"

	" > $pre.pbs
done
