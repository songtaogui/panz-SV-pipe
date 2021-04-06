#!/usr/bin/env bash
usage="
------------------------------------------------------------
A script for getting the exact seq from Assemblytics results
------------------------------------------------------------
USAGE:
	$0 [options below]

	1. ABS_PATH of vcf file of the Assemblytics results,named:
	   PREFIX_asmSV_mft.vcf
	   (converted by SURVIVOR and reformated by local script)
	2. ABS_PATH of Raw Assemblytics results dir
	3. ABS_PATH of REF_fasta
	4. ABS_PATH of Query_fasta
------------------------------------------------------------
                                            Songtao Gui
                                       songtaogui@sina.com
"
if [[ $# -ne 4 ]]; then 
	echo "$usage" >&2
	exit 1
fi

invcf=$1
asmdir=$2
reffa=$3
queryfa=$4
export pre=$(basename $invcf)
pre=${pre%%_asmSV*}
asmbed=$asmdir/${pre}_asm.Assemblytics_structural_variants.bed

echo "Dealing with $pre ..." >&2

############ Check if all related files exists ###########
num_related_file=1;
for related_file in $invcf $reffa $queryfa $asmbed
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
#################################################################

echo "extract ref cord of $pre ..." >&2

cat $invcf | grep -v "^#" | perl -F"\t" -lane '
$chr1=$F[0];
$start=$F[1]; 
$ID=$F[2];
if(/CHR2=(.*?);END=(\d+);/){
	$chr2=$1;$end=$2;
	next if $chr1 ne $chr2;
	printf "%s\t%s#%s#%s\n",$ID,$chr1,$start,$end;
}
' > ${pre}_id_cord_vcf.tsv.tmp

if [[ $? -ne 0 ]] ; then 
	echo "[ERROR] --> vcf cord: non-zero exit." >&2
	# rm -f xxx.outputs
	exit 1
fi

echo "extract ref & query cord of $pre from Assemblytics BED ..." >&2

cat $asmbed | grep -v "^#" | perl -F"\t" -lane '
$ref_cord=sprintf("%s:%s-%s:%s",$F[0],$F[1],$F[2],$F[5]);
$id=sprintf("%s#%s#%s",$F[0],$F[1],$F[2]);
$query_cord=$F[9];
print "$id\t$ref_cord\t$query_cord"
' > ${pre}_ref_query_cord_bed.tsv.tmp

if [[ $? -ne 0 ]] ; then 
	echo "[ERROR] --> bed cord: non-zero exit." >&2
	# rm -f xxx.outputs
	exit 1
fi

echo "join and get ref & query bed ..." >&2

csvtk join -tTH -f "1;2" ${pre}_id_cord_vcf.tsv.tmp ${pre}_ref_query_cord_bed.tsv.tmp -o ${pre}_id_ref_query.tsv
if [[ $? -ne 0 ]] ; then 
	echo "[ERROR] --> join: non-zero exit." >&2
	# rm -f xxx.outputs
	exit 1
fi
echo "get ref seq ..." >&2
cat ${pre}_id_ref_query.tsv | perl -F"\t" -lane '
BEGIN{
	print STDERR "Generate ref and query bed file for $ENV{pre}";
	open(REF,">$ENV{pre}_ref.bed") or die;
	open(QUERY,">$ENV{pre}_query.bed") or die;
}
($id,undef,$ref_cord,$query_cord)=@F;
if($ref_cord=~s/(.*):(\d+)-(\d+):(.)/$1\t$2\t$3\t$4/){
	$rc=$1;$rs=$2;$re=$3;$rstr=$4;
}
if($query_cord=~s/(.*):(\d+)-(\d+):(.)/$1\t$2\t$3\t$4/){
	$qc=$1;$qs=$2;$qe=$3;$qstr=$4;
}
$,="\t";
print REF $rc,$rs-1,$re,$id,"1",$str;
print QUERY $qc,$qs-1,$qe,$id,"1",$str;
'
if [[ $? -ne 0 ]] ; then 
	echo "[ERROR] --> get ref qeury bed: non-zero exit." >&2
	# rm -f xxx.outputs
	exit 1
fi

seqkit subseq --bed ${pre}_ref.bed $reffa | seqkit fx2tab | perl -lane '$,="\t";print @F[1,2]' > ${pre}_ref_seq.tsv.tmp &&\
seqkit subseq --bed ${pre}_query.bed $queryfa | seqkit fx2tab | perl -lane '$,="\t";print @F[1,2]' > ${pre}_query_seq.tsv.tmp &&\
csvtk join -tTH -f 1 ${pre}_ref_seq.tsv.tmp ${pre}_query_seq.tsv.tmp -o ${pre}_ref_query_seq.tsv 

if [[ $? -ne 0 ]] ; then 
	echo "[ERROR] --> get ref alt seq: non-zero exit." >&2
	# rm -f xxx.outputs
	exit 1
fi

echo "Get new vcf with ref alt seq ..." >&2

cat $invcf | perl -lane '
BEGIN{
	open(IN,"$ENV{pre}_ref_query_seq.tsv") or die;
	$,="\t";
	while(<IN>){
		chomp;
		($id,$ref,$query)=split(/\t/,$_);
		$ref_seq{$id}=$ref;
		$query_seq{$id}=$query;
	}
}
/^#/ && print && next;
$F[3]=$ref_seq{$F[2]} if $ref_seq{$F[2]};
$F[4]=$query_seq{$F[2]} if $query_seq{$F[2]};
print @F;
' > ${pre}_with_seq.vcf
if [[ $? -ne 0 ]] ; then 
	echo "[ERROR] --> with seq vcf: non-zero exit." >&2
	# rm -f xxx.outputs
	exit 1
fi

echo "all done !"

