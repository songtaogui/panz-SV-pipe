ls $WORK/RMdup_MEM/MarkDup_MEM/*.MEM.mkdup.bam | while read bam
do
pre=$(basename $bam)
pre=${pre%%.*}
name=SVtyper_run
cpu=1
echo "#BSUB -J ${name}_${pre}
#BSUB -n $cpu
#BSUB -R span[hosts=1]
#BSUB -R rusage[mem=1GB]
#BSUB -o ${name}_${pre}.out
#BSUB -e ${name}_${pre}.out
#BSUB -q \"normal\"

svtyper -i PANZ_SV_Caller_merged.vcf -B $bam -l ${bam%%.*}.json -w ${pre}_evidence.bam -o ${pre}_SVTYPER_geno.vcf 1>${pre}_svtyper.log 2>&1

" > ${name}_${pre}.lsf
done
