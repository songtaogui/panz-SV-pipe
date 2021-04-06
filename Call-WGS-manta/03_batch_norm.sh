ref=$WORK/ref/Zea_mays.AGPv4.dna.toplevel.fa
for i in $WORK/AMP_MEM_SV_TEs/manta/MKdup/manta_filter/*_manta_bystypr.vcf.gz
do
	pre0=${i%%_manta_bystypr.vcf.gz}
	pre=`basename $pre0`
	name=nm
	cpu=1
	echo "
#PBS -N ${name}_$pre
#PBS -l nodes=1:ppn=$cpu
#PBS -j oe
#PBS -q batch
#PBS -V

cd  \$PBS_O_WORKDIR

echo \"[\`date +\"%m-%d %H:%M\"\`] ---------> ALL START\"

bcftools norm --threads $cpu -f $ref -o ${pre}_norm.vcf.gz -O z -c x $i

echo \"[\`date +\"%m-%d %H:%M\"\`] ---------> ALL DONE\"

	" > $pre.pbs
done

