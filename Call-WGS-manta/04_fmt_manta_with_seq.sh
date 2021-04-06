#!/usr/bin/env bash
usage="
------------------------------------------------------------
combine raw manta with convert allele out, rmdup and rename 
------------------------------------------------------------
Dependency: csvtk bcftools
------------------------------------------------------------
USAGE:
	$0 <raw_manta.vcf> <convert_allele.vcf.gz> <cpus>
------------------------------------------------------------
                                            Songtao Gui
                                       songtaogui@sina.com
"
if [[ $# -ne 3 ]]; then 
	echo "$usage" >&2
	exit 1
fi

raw_vcf=$1
cvt_vcf=$2
cpus=$3
############ Check if all related files exists ###########
num_related_file=1;
for related_file in $raw_vcf $cvt_vcf
do
	if [[ ! -s "$related_file" ]]; then
		printf "[ERROR] --> No file: %s \n" $related_file >&2
		let num_related_file++
	fi
done

if [ $num_related_file -eq 1 ];then
	echo "All related files were found. Proceeding ..." >&2
else
	echo "Check if you miss somthing." >&2
	exit 1;
fi

export pre=$(basename $raw_vcf)
pre=${pre%%.*}
pre=${pre%%_manta*}
# echo "$pre" >&2
# grep "^#" $raw_vcf > ${pre}_manta_addSeq.vcf &&\
# cvt rcd
zcat $cvt_vcf | perl -F"\t" -lane '$,="\t";/^#/ && print && next;$F[2]="$ENV{pre}_$F[2]";print @F;' >${pre}_manta_addSeq.vcf &&\
# cvt id 
csvtk cut -tTH -f 3 -j $cpus -tTH $cvt_vcf -o ${pre}_cvt.id &&\
# non cvt rcd
csvtk grep -tTH -j $cpus -v -f 3 -P ${pre}_cvt.id $raw_vcf | perl -F"\t" -lane '$,="\t";/^#/ && print && next;$F[2]="$ENV{pre}_$F[2]";print @F;' >> ${pre}_manta_addSeq.vcf &&\
# sort
bcftools sort ${pre}_manta_addSeq.vcf -o ${pre}_manta_addSeq_sort.vcf
if [[ $? -ne 0 ]] ; then 
	echo "[ERROR] --> CMD: non-zero exit." >&2
	rm -f ${pre}_manta_addSeq_sort.vcf
	exit 1
fi
if [ -s "${pre}_manta_addSeq_sort.vcf" ];then
	rm -f ${pre}_manta_addSeq.vcf ${pre}_cvt.id
fi
echo "Done for $pre" >&2