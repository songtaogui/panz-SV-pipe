#!/usr/bin/env bash

set -euo pipefail

# usage check_files_exists
function check_files_exists(){
num_related_file=1;
    for related_file in  "$@"
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
}
function gst_log () {
    local info=$1
    echo "[$(date +'%y-%m-%d %H:%M')] $info" >&2
}
function gst_err () {
    local info=$1
    echo "[ERROR] --> $info" >&2
}
function gst_warn () {
    local info=$1
    echo "[WARNING] --> $info" >&2
}

usage="
------------------------------------------------------------
Filter NRINS VCFs and get sequence:
    Keep only both end anchored to same chr records.
    Filter out ref or query region > 1 Mb.
------------------------------------------------------------
USAGE:
    bash $(basename $0) <NRINS.vcf> <Ref_seq.fa> <Query_seq.fa> <prefix>
------------------------------------------------------------
Author: Songtao Gui
E-mail: songtaogui@sina.com
"
if [[ $# -ne 4 ]]; then 
    echo "$usage" >&2
    exit 1
fi

vcf=$1
ref=$2
query=$3
pre=$4

gst_log "Dealing with ${pre} ..."

check_files_exists $vcf $ref $query

#0 10
#1 84296222
#2 BEref_DiffStr # 5to3 # R_L # chr10_100704230_100811999_CML247_scaftig0176907_53808R56683 # RF
#3 N
#4 <INS>
#5 .
#6 PASS
#7 SVTYPE=INS;SVMETHOD=MMref;CHR2=10;END=97041204;SVLEN=2876;CT=5to3;
#8 GT
#9 1/1
NRINS_info=$(
    cat $vcf | perl -F"\t" -lane '
        BEGIN{
            $,="\t";
        }
        /^#/ && next;
        $chr1=$F[0];
        $start=$F[1]-1;
        $ID=$F[2];
        ($type,$orient,$part,$qid,$str)=split(/#/,$ID);
        $qchr="";$qstart=0;$qend=0;
        if($qid=~/^(.*)_(\d+)R(\d+)$/){
            $qchr=$1;
            $qstart=$2-1;
            $qend=$3;
        }else{
            die("[ERROR]-->Wrong query id: $qid for $F[2]");
        }
        if($F[7] =~ /CHR2=(.*?);END=(\d+);/){
            $chr2=$1;
            $end=$2;
        }else{
            die("[ERROR]-->Wrong info format for $F[2]");
        }
        #? same chr
        $chr1 ne $chr2 && next;
        $r_len=$end-$start;
        $q_len=$qend-$qstart;
        #? no more than 1Mb
        $r_len >= 1e6 && next;
        $q_len >= 1e6 && next;
        #? Both end
        $type !~ /^BE/ && next; 
        # svtype: DEL/INS
        $SVTYPE="INS";
        $SVTYPE="DEL" if $r_len > $q_len;
        $STRAND="+";
        $STRAND="-" if $str=~/^R/;
        $SVLEN=$q_len-$r_len;
        #? 0:ID 1:rc 2:rs 3:re 4:str 5:svtype 6:qc 7:qs 8:qe 9:svlen
        print $ID,$chr1,$start,$end,$STRAND,$SVTYPE,$qchr,$qstart,$qend,$SVLEN;
    '
)

gst_log "Getting ref and query sequences ..."

# >>>>>>>>>>>>>>>>>>>>>>>> get sequences >>>>>>>>>>>>>>>>>>>>>>>>
# get beds
NRINS_ref_bed=$( echo "$NRINS_info" | perl -F"\t" -lane '$,="\t";print @F[1,2,3,0,9],"+";')
NRINS_query_bed=$( echo "$NRINS_info" | perl -F"\t" -lane '$,="\t";print @F[6,7,8,0,9,4];')
# get seq
seqkit subseq --bed <(echo "$NRINS_ref_bed") $ref | seqkit fx2tab | perl -lane '$,="\t"; print @F[-2,-1];' > ${pre}_ref_seq.tsv
seqkit subseq --bed <(echo "$NRINS_query_bed") $query | seqkit fx2tab | perl -lane '$,="\t"; print @F[-2,-1];' > ${pre}_query_seq.tsv
# join
csvtk join -tTH -f 1 <(echo "$NRINS_info") ${pre}_ref_seq.tsv ${pre}_query_seq.tsv -o ${pre}_NRINS_info.tsv
# <<<<<<<<<<<<<<<<<<<<<<<< get sequences <<<<<<<<<<<<<<<<<<<<<<<<

# >>>>>>>>>>>>>>>>>>>>>>>> get final vcf >>>>>>>>>>>>>>>>>>>>>>>>
gst_log "Get final vcf ..."

grep "^#" $vcf > ${pre}_NRINS_withSeq.vcf
#? 0:ID 1:rc 2:rs 3:re 4:str 5:svtype 6:qc 7:qs 8:qe 9:svlen 10:refseq 11:queryseq
#0 10
#1 84296222
#2 BEref_DiffStr # 5to3 # R_L # chr10_100704230_100811999_CML247_scaftig0176907_53808R56683 # RF
#3 N
#4 <INS>
#5 .
#6 PASS
#7 SVTYPE=INS;SVMETHOD=MMref;CHR2=10;END=97041204;SVLEN=2876;CT=5to3;
#8 GT
#9 1/1
cat ${pre}_NRINS_info.tsv | perl -F"\t" -lane '
    BEGIN{$,="\t";}
    ($id,$chr,$posi,$end,$str,$svtype,$qc,$qs,$qe,$svlen,$ref,$alt)=@F;
    $ct="5to3";
    $ct="3to5" if $str eq "-";
    print $chr,$posi+1,$id,$ref,$alt,".","PASS","SVTYPE=$svtype;SVMETHOD=PANZ_NRINS;CHR2=$chr;END=$end;SVLEN=$svlen;CT=$ct","GT","1/1";
' >> ${pre}_NRINS_withSeq.vcf
# <<<<<<<<<<<<<<<<<<<<<<<< get final vcf <<<<<<<<<<<<<<<<<<<<<<<<

gst_log "Done. Resulted vcf file: ${pre}_NRINS_withSeq.vcf"

rm -f ${pre}_ref_seq.tsv ${pre}_query_seq.tsv ${pre}_NRINS_info.tsv
