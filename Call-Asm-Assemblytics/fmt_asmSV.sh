#!/usr/bin/env bash
usage="
------------------------------------------------------------
A script for formatting asmSV.vcf IDs
------------------------------------------------------------
USAGE:
	$0 <input.vcf (result from run_asmSV.sh)>

Parse input names to get [SP]:
	input names should be in format of "[SP]_asmSV_XXXX.vcf"
------------------------------------------------------------
                                            Songtao Gui
                                       songtaogui@sina.com
"
if [[ $# -ne 1 ]]; then 
	echo "$usage" >&2
	exit 1
fi

invcf=$1

if [ ! -s "$invcf" ];then
	echo "No file: $invcf" >&2
	exit 1 
fi

sp0=$(basename $invcf)
export sp=${sp0%%_asmSV*}
echo "Dealing with $sp ..." >&2
cat $invcf | perl -lane '
/^#CHR/ && s/_asm.Assemblytics_structural_variants.bed//g && print && next;
/^#/ && print && next;
$,="\t";
if($F[7]=~/SVTYPE=(\w+?);/){
	$svtype=$1;
	$h{$svtype}++;
	$id=sprintf("%s_%s%07s",$ENV{sp},$svtype,$h{$svtype});
	$F[2]=$id;
	print @F;
}
' > ${sp}_asmSV_mft.vcf
if [[ $? -ne 0 ]] ; then 
	echo "[ERROR] --> CMD: non-zero exit." >&2
	# rm -f xxx.outputs
	exit 1
fi

echo "Done."

