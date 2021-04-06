#!/usr/bin/env bash
usage="
------------------------------------------------------------
Filter survivor merge vcf files by support samples
------------------------------------------------------------
USAGE:
	$0 <survivor_merge.vcf> <supp cut-off> <Prefix>
	Prefix: for non-sample out vcfs, usually use Method, 
	eg: Manta
------------------------------------------------------------
                                            Songtao Gui
                                       songtaogui@sina.com
"
if [[ $# -ne 3 ]]; then 
	echo "$usage" >&2
	exit 1
fi

invcf=$1
export suppnum=$2
export pre=$3

export outpre=$(basename $invcf)
outpre=${outpre%%.vcf}

cat $invcf | perl -F"\t" -lane '
BEGIN{
	$,="\t";
	open(ALL,">$ENV{outpre}_SUPP$ENV{suppnum}_AllSamples.vcf") or die;
	open(NOSAMPLE,">$ENV{outpre}_SUPP$ENV{suppnum}_NoSample.vcf") or die;
}
	print ALL if /^##/;
	print NOSAMPLE if /^##/;
	print ALL @F if /^#CHR/;
	print NOSAMPLE @F[0..8],"$ENV{pre}" if /^#CHR/;
	#$F[2]="$ENV{pre}".$F[2];
	$F[7]=~s/(SUPP_VEC=\d+;)//;
	if($F[7] =~ /SUPP=(\d+);/){
		$supp=$1;
		print NOSAMPLE @F[0..7],"GT","1/1" if $supp >= $ENV{suppnum};
		print ALL @F if $supp >= $ENV{suppnum};
}' >&2



