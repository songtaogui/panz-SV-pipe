#!/usr/bin/env bash
usage="
------------------------------------------------------------
A script for get ref alt seq from Gridss vcfs
------------------------------------------------------------
USAGE:
	$0 <Gridss_smp.vcf> <Gridss_raw/annot.vcf> <ref.fa> <min SV len>

Gridss_smp.vcf: SV records of simple type, generated from gst_gridss_anno_filt_bed.R
Gridss_raw/annot.vcf" raw Gridss output vcf
------------------------------------------------------------
                                            Songtao Gui
                                       songtaogui@sina.com
"
if [[ $# -ne 4 ]]; then 
	echo "$usage" >&2
	exit 1
fi


#smp_vcf=/public/home/stgui/work/22_SV/Gridss/vcf/Q100_grfmt.vcf
#annot_dir=/public/home/stgui/work/22_SV/Gridss/gridss_seq/anno_vcf
#annot_vcf=${annot_dir}/${pre}_annotated.vcf.gz
#reffa=/public/home/stgui/work/ref/B73_V4/Zea_mays.AGPv4.dna.toplevel.fa
#export minlen=50
smp_vcf=$1
annot_vcf=$2
reffa=$3
export minlen=$4
export pre=""
pre=$(basename $smp_vcf)
pre=${pre%%_*}

############ Check if all related files exists ###########
num_related_file=1;
for related_file in $smp_vcf $annot_vcf $reffa
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
############ Check if all related files exists ###########
echo "
------------------------------------------------
          START Dealing with $pre
------------------------------------------------
SMP_VCF :   $smp_vcf
RAW_VCF :   $annot_vcf
REF_FA  :   $reffa
MINSVLEN:   $minlen
------------------------------------------------
" >&2

echo "
------------------------------------------------
[$(date +%m-%d_%H:%M)] Filter out SV < $minlen and get join seq 
------------------------------------------------
" >&2
if [ -s "${pre}_join_seq_gt${minlen}.tsv" ];then
	echo "SKIP seq info extracting, file exists: ${pre}_join_seq_gt${minlen}.tsv." >&2
else
	# get smp id (no prefix), filter out SVlen < $minlen records
	cat $smp_vcf | perl -F"\t" -lane '
		/^#/ && next;
		$id=$F[2];
		if($F[7]=~/SVLEN=(.*?);/){
			$len=abs($1);
			print $id if $len >= $ENV{minlen};
		}
	' | sed "s/${pre}_//;" > ${pre}_smp_id_gt${minlen}.id_o.tmp &&\
	cat ${pre}_smp_id_gt${minlen}.id_o.tmp | sed 's/o$/h/' > ${pre}_smp_id_gt${minlen}.id_h.tmp &&\
	csvtk grep -tTH -f 3 -P ${pre}_smp_id_gt${minlen}.id_o.tmp $annot_vcf | perl -F"\t" -lane '
		BEGIN{
			print "ID\tSV_O\tLEN_O\tCHR_O\tPOS_O\tREF_O\tALT_O";
			$,="\t";
		}
		$F[2]=~s/[oh]$//;
		$svtype="UNDEF";
		$svtype=$1 if $F[7]=~/SIMPLE_TYPE=(\w+)/;
		$svlen=0;
		$svlen=abs($1) if $F[7]=~/SVLEN=(.*?);/;		
		print $F[2],$svtype,$svlen,@F[0,1,3,4];
	' > ${pre}_seq_o_gt${minlen}.tmp
	csvtk grep -tTH -f 3 -P ${pre}_smp_id_gt${minlen}.id_h.tmp $annot_vcf | perl -F"\t" -lane '
		BEGIN{
			print "ID\tSV\tLEN\tCHR_H\tPOS_H\tREF_H\tALT_H";
			$,="\t";
		}
		$F[2]=~s/[oh]$//;
		$svtype="UNDEF";
		$svtype=$1 if $F[7]=~/SIMPLE_TYPE=(\w+)/;
		$svlen=0;
		$svlen=abs($1) if $F[7]=~/SVLEN=(.*?);/;		
		print $F[2],$svtype,$svlen,@F[0,1,3,4];
	' > ${pre}_seq_h_gt${minlen}.tmp
	csvtk join -tTH -f 1 ${pre}_seq_o_gt${minlen}.tmp ${pre}_seq_h_gt${minlen}.tmp -o ${pre}_join_seq_gt${minlen}.tsv
	if [[ $? -ne 0 ]] ; then 
		echo "[ERROR] --> join_seq: non-zero exit." >&2
		rm -f ${pre}_join_seq_gt${minlen}.tsv
	else
		rm -f ${pre}_*.tmp
	fi
	echo "Done! Result: ${pre}_join_seq_gt${minlen}.tsv" >&2
fi
echo "
------------------------------------------
[$(date +%m-%d_%H:%M)] GET information bed 
------------------------------------------
" >&2
if [ -s "${pre}_ref_alt_seq_gt${minlen}.bed" ];then
	echo "skip running. File exists:${pre}_ref_alt_seq_gt${minlen}.bed " >&2
else
	cat ${pre}_join_seq_gt${minlen}.tsv | sed '1d' | perl -F"\t" -lane '
		BEGIN{ $,="\t"; }
		($ID,$SV_O,$LEN_O,$CHR_O,$POS_O,$REF_O,$ALT_O,$SV_H,$LEN_H,$CHR_H,$POS_H,$REF_H,$ALT_H)=@F;
		# skip inversions because bayestype can not geno it
		$SV_O eq "INV" && next;
		$SV_O ne $SV_H && print STDERR "[WARNING] --> SV_O and SV_H not same for $ENV{pre}:$ID , Skip it" && next;
		#die("[ERROR] --> LEN_O and LEN_H not same for $ENV{pre}:$ID") if $LEN_O ne $LEN_H;
		if( $ALT_O =~ /\[/ && $ALT_H =~ /\]/ ){
			# forward [>>>>>>> seqs <<<<<<] backward 
			$forward_rcd = $ALT_H; 
			$backward_rcd = $ALT_O; 
		} elsif ( $ALT_H =~ /\[/ && $ALT_O =~ /\]/ ){
			$forward_rcd = $ALT_O; 
			$backward_rcd = $ALT_H;
		} else {
			print STDERR "[WARNING] --> Wrong SV orient for $ID , SKIP it";
			next;
		}
		#f ]1:16728055] CCTCCCTGAAT A-----------		
		#b -----------T CCTCCCTGAAT [1:16728290[
		($f_seq,$f_coord) = (split(/\]/,$forward_rcd))[-1,-2];
		($f_chr,$f_posi) = split(/:/,$f_coord);
		($b_seq,$b_coord) = (split(/\[/,$backward_rcd))[0,1];
		($b_chr,$b_posi) = split(/:/,$b_coord);
		# print STDERR "$f_seq,$f_coord,$forward_rcd";
		# print STDERR "$b_seq,$b_coord,$backward_rcd";
		$f_chr ne $b_chr && print STDERR "[WARNING] --> CHR: forward $f_chr and backward $b_chr not same for $ENV{pre}:$ID, SKIP it" && next;
		($f_seq_common,$f_seq_single)=($1,$2) if $f_seq =~ /^(\w*)(\w)$/;
		($b_seq_single,$b_seq_common)=($1,$2) if $b_seq =~ /^(\w)(\w*)$/;
		$diff_len_O=abs($b_posi-$f_posi)+1-$LEN_O;
		$diff_len_H=abs($b_posi-$f_posi)+1-$LEN_H;
		$alt_seq_same=join("",$b_seq,$f_seq_single);
		$info=join("",$forward_rcd,"#"x10,$backward_rcd);
		$svlen_diff_add_1_O=$diff_len_O - length($alt_seq_same)+1;
		$svlen_diff_add_1_H=$diff_len_H - length($alt_seq_same)+1;
		$f_sss=$f_seq_common;
		$f_sss="UNDEF" if !$f_seq_common;
		$b_sss=$b_seq_common;
		$b_sss="UNDEF" if !$b_seq_common;
		$check_posi_svlen=0;
		$check_posi_svlen=1 if $svlen_diff_add_1_O == 1 or $svlen_diff_add_1_H == 1;
		$type="UNDEF";
		if( $check_posi_svlen == 1 and $f_sss eq $b_sss ){
			# perfect match
			$type="perfect";
			$alt_seq_final=$alt_seq_same;
		}else{
			if($f_sss eq $b_sss){
				# raw svlen was wrong but common seqs are same
				$type="re_svlen";
				$alt_seq_final=$alt_seq_same;
			}else{
				# common seqs are diff but one of them is UNDEF
				$type="unsure";
				if($f_sss eq "UNDEF" or $b_sss eq "UNDEF"){
					$alt_seq_final=join("",$b_seq,$f_seq);	
				}else{
					$alt_seq_final=join("",$b_seq,"N"x10,$f_seq);
				}
			}
		}
		if( $f_posi > $b_posi ){
			$bed_str="-";
			$bed_start=$b_posi;
			$bed_end=$f_posi;
			$bed_chr=$f_chr;
		}else{
			$bed_str="+";
			$bed_start=$f_posi;
			$bed_end=$b_posi;
			$bed_chr=$f_chr;
		}
		$bed_start = 1 if $bed_start < 1;
		$bed_id=join("",$ENV{pre},"_",$ID,"o","#ALT#",$alt_seq_final);
		print $bed_chr,$bed_start-1,$bed_end,$bed_id,$type,$bed_str,$info;
	' > ${pre}_ref_alt_seq_gt${minlen}.bed
	if [[ $? -ne 0 ]] ; then 
		echo "[ERROR] --> get ref_alt_seq for $pre: non-zero exit." >&2
		rm -f ${pre}_ref_alt_seq_gt${minlen}.bed
		exit 1
	fi
	echo "Done! Result: ${pre}_ref_alt_seq_gt${minlen}.bed" >&2
fi
echo "
------------------------------------------
[$(date +%m-%d_%H:%M)] GET REF Sequences 
------------------------------------------
" >&2
if [ ! -s "${pre}_ref_alt_seq_gt${minlen}.info" ];then
	seqkit subseq --bed ${pre}_ref_alt_seq_gt${minlen}.bed $reffa | seqkit fx2tab | perl -lane '
	BEGIN{
		print "#Chr\tstart\tend\tstrand\tID\talt\tref";
		$,="\t";
	}
	$F[0]=~s/^(.*)_(\d+)-(\d+):(.)$/$1\t$2\t$3\t$4/g;
	$F[1]=~s/#ALT#/\t/g;
	print @F;' | \
	perl -F"\t" -lane 'BEGIN{$,="\t";} print @F[0,1,2,3,4,6,5];' > ${pre}_ref_alt_seq_gt${minlen}.info
	if [[ $? -ne 0 ]] ; then 
		echo "[ERROR] --> get ref alt info: non-zero exit." >&2
		rm -f ${pre}_ref_alt_seq_gt${minlen}.info
		exit 1
	fi
	echo "Done! Result: ${pre}_ref_alt_seq_gt${minlen}.info" >&2
else
	echo "Skip running. File exists: ${pre}_ref_alt_seq_gt${minlen}.info" >&2
fi

echo "
----------------------------------------------
[$(date +%m-%d_%H:%M)] Generate final VCF file
----------------------------------------------
" >&2
if [ -s "${pre}_Gridss_add_seq_gt${minlen}.vcf" ];then
	echo "Skip running. File exists: ${pre}_Gridss_add_seq_gt${minlen}.vcf" >&2
else
	cat $smp_vcf | perl -F"\t" -lane '
	BEGIN{
		$,="\t";
		open(IN,"$ENV{pre}_ref_alt_seq_gt$ENV{minlen}.info") or die("Cannot open file: $ENV{pre}_ref_alt_seq_gt$ENV{minlen}.info");
		while(<IN>){
			chomp;
			($c,$s,$e,$str,$id,$ref,$alt)=split(/\t/,$_);
			$h_chr{$id}=$c;
			$h_posi{$id}=$s+1;
			$h_end{$id}=$e;
			$h_svlen{$id}=length($alt)-length($ref);
			$h_strand{$id}=$str;
			$h_ref{$id}=$ref;
			$h_alt{$id}=$alt;
		}
	}
	/^#/ && print && next;

	($CHROM,$POS,$ID,$REF,$ALT,$QUAL,$FILTER,$INFO,$FORMAT,$SAMPLE)=@F;
	$CIPOS="";$CIPOS=$1 if $INFO =~ /(CIPOS=\d+,\d+)/ ;
	$CIEND="";$CIEND=$1 if $INFO =~ /(CIEND=\d+,\d+)/ ;
	$SVTYPE="";$SVTYPE=$1 if $INFO=~ /(SVTYPE=.*?);/ ;
	$SVLEN=abs($1) if $F[7]=~/SVLEN=(.*?);/;
	$SVLEN < $ENV{minlen} && next;
	if($h_chr{$ID}){
		$CHROM=$h_chr{$ID};
		$POS=$h_posi{$ID};
		$REF=$h_ref{$ID};
		$ALT=$h_alt{$ID};
		$INFO="CHR2=$h_chr{$ID};END=$h_end{$ID};$SVTYPE;SVLEN=$h_svlen{$ID};$CIPOS;$CIEND";
	}
	print $CHROM,$POS,$ID,$REF,$ALT,$QUAL,$FILTER,$INFO,$FORMAT,$SAMPLE;
	' > ${pre}_Gridss_add_seq_gt${minlen}.vcf
fi
echo "DONE ! Result: ${pre}_Gridss_add_seq_gt${minlen}.vcf" >&2
echo "
----------------------------------------------
[$(date +%m-%d_%H:%M)] ALL DONE !!
----------------------------------------------
" >&2
