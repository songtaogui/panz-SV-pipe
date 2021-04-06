#!/usr/bin/env bash

set -euo pipefail

# >>>>>>>>>>>>>>>>>>>>>>>> function gst_log_err_warn >>>>>>>>>>>>>>>>>>>>>>>>
function gst_log () {
    local info=$1
    echo "$(tput setaf 6)[$(date +'%y-%m-%d %H:%M')]$(tput sgr0) $info" >&2
}
function gst_err () {
    local info=$1
    echo "$(tput setaf 1)[ERROR] --> $info$(tput sgr0)" >&2
}
function gst_warn () {
    local info=$1
    echo "$(tput setaf 5)[WARNING] --> $info$(tput sgr0)" >&2
}
# <<<<<<<<<<<<<<<<<<<<<<<< function gst_log_err_warn <<<<<<<<<<<<<<<<<<<<<<<<
usage="
------------------------------------------------------------
USAGE:
    bash $(basename $0) <input: 00_Merge_all_gt0.vcf>
------------------------------------------------------------
A script to get represent record for Merged All SV for PANZ
rename the SV, get represent ref-alt seq, and generate these
outputs:
    PANZ_SV_Reformat.vcf
    PANZ_SV_SVtyper.vcf
    PANZ_SV_Delly.vcf
    PANZ_SV_Bayestyper.vcf
    PANZ_SV_ID_each_call.tsv

Record Keep Priority:
    1BSpSV > 2AsmREF > 3AsmAMP > 4NRINS > 5Gridss > 6Manta > 7Delly
Rename rules:
    PZ00001aSVnum[SVTYPE]Ref_Chr1#Posi1#Chr2#Pos2

Note: input vcf should been merged according to the priority order above.
(The script parsed the SUPP_VEC=010101 information without checking the order.)
------------------------------------------------------------
Author: Songtao Gui
E-mail: songtaogui@sina.com
"
if [[ $# -ne 1 ]]; then
    echo "$usage" >&2
    exit 1
fi

vcf=$1
ref=/public/home/stgui/work/ref/B73_V4/Zea_mays.AGPv4.dna.toplevel.fa
snp_vcf=
indel_vcf=
delly_header=/public/home/stgui/work/22_SV/00_Merge_all/Delly_header.txt

if [ ! -s "$vcf" ]; then
    gst_log "No file: $vcf"
    exit 1
fi


# >>>>>>>>>>>>>>>>>>>>>>>>  PANZ_SV_Reformat.vcf >>>>>>>>>>>>>>>>>>>>>>>>
# parse each line function

gst_log "Generating NewIDs, and reformat vcf files ..."

rm -f PANZ_SV_ID_each_call.tsv
cat $vcf | perl -F"\t" -lane '
    BEGIN{
        $,="\t";
        @sample=("BSpSV","asm_REF","asm_AMP","NRINS","gridss","manta","Delly");
        open(OUT,">>PANZ_SV_ID_each_call.tsv") or die("ERROR: Cannot open fie: PANZ_SV_ID_each_call.tsv");
    }
    #? parse samples and print header
    /^##/ && print && next;
    if(/^#CHR/) {
        print OUT "PANZ_SV",@sample;
        print @F[0..8],"PANZSV";
    }else{
        #? parse SUPP_VEC
        $supp_vec=$1 if $F[7]=~/SUPP_VEC=(\d+);/;
        @vec=split("",$supp_vec);
        #? parse info features
        $svtype=$1 if $F[7]=~/SVTYPE=(.*?);/;
        $cipos="0,0";$ciend="0,0";$strands="+-";
        $cipos=$1 if $F[7]=~/CIPOS=(.*?);/;
        $ciend=$1 if $F[7]=~/CIEND=(.*?);/;
        $strands=$1 if $F[7]=~/STRANDS=([\+\-]*)/;
        #? parse each Method info,and save Names into hashes
        %hash_ID=();
        %hash_LN=();
        %hash_ST=();
        %hash_RAL=();
        %hash_AAL=();
        %hash_CO=();
        for ($i = 0; $i <= $#sample; $i++) {
            @FORMAT=split(/:/,$F[$i+9]);
            #? parse ID incase there is : in it
            $tmp_sample=$sample[$i];
            $tmp_ID=join(":",@FORMAT[7..$#FORMAT-3]);
            #? save to each hash
            $hash_ID{$tmp_sample}=$tmp_ID;
            $hash_LN{$tmp_sample}=$FORMAT[2];
            $hash_ST{$tmp_sample}=$FORMAT[4];
            $hash_RAL{$tmp_sample}=$FORMAT[-3];
            $hash_AAL{$tmp_sample}=$FORMAT[-2];
            $hash_CO{$tmp_sample}=$FORMAT[-1];
        }
        #? get the first record in @vec, and print out
        for ($i = 0; $i <= $#vec; $i++) {
            if( $vec[$i] > 0 ) {
                $used_geno=$F[$i+9];
                $used_sample=$sample[$i];
                $used_ID=$hash_ID{$used_sample};
                $used_LN=$hash_LN{$used_sample};
                $used_ST=$hash_ST{$used_sample};
                $used_RAL=$hash_RAL{$used_sample};
                $used_AAL=$hash_AAL{$used_sample};
                $used_CO=$hash_CO{$used_sample};
                $used_posi=(split(/,/,$used_CO))[0];
                if($used_posi =~ /^(.*)_(\d+)-(.*)_(\d+)$/){
                    $used_chr1=$1;
                    $used_start=$2;
                    $used_chr2=$3;
                    $used_end=$4;
                }else{
                    print STDERR "$F[2],[$used_sample]",@sample;#! test
                    die("ERROR: Wrong CO for $used_ID");
                }
                # ? use raw ref alt info if the ref and alt in the record is NA
                # ! print STDERR @F[3,4] if $used_RAL eq "NA"; # ! test
                $used_RAL=$F[3] if $used_RAL eq "NA";
                $used_AAL=$F[4] if $used_AAL eq "NA";
                last;
            }
        }
        #? ChangeID to PZ00001aSVnum[SVTYPE]Ref_Chr1#Posi1#Chr2#Pos2
        #? print out the represent SV info
        # // ++$count_svtype{$svtype};
        # // $SV_TYPE_NUM=$count_svtype{$svtype};
        $SV_TYPE_NUM++;
        $NEW_SV_ID=sprintf("PZ00001aSV%08d%s#%s#%s#%s#%s",$SV_TYPE_NUM,$svtype,$used_chr1,$used_start,$used_chr2,$used_end);
        $NEW_SV_INFO=sprintf("SVLEN=%s;SVTYPE=%s;SVMETHOD=PANZ;CHR2=%s;END=%s;CIPOS=%s;CIEND=%s;STRANDS=%s",$used_LN,$svtype,$used_chr2,$used_end,$cipos,$ciend,$used_ST);
        if($F[4] !~ /\W/ and $F[3] !~ /\W/){
            #? raw record was allready seq-awarable, just changeID and print out the raw record;
            $F[2]=$NEW_SV_ID;
            $F[7]=~s/SUPP.*?;//g;
            print @F[0..7],"GT","1/1";
        }else{
            #? print new record parsed above
            print $used_chr1,$used_start,$NEW_SV_ID,$used_RAL,$used_AAL,$F[5],$F[6],$NEW_SV_INFO,"GT","1/1";
        }
        print OUT $NEW_SV_ID,$hash_ID{$sample[0]},$hash_ID{$sample[1]},$hash_ID{$sample[2]},$hash_ID{$sample[3]},$hash_ID{$sample[4]},$hash_ID{$sample[5]},$hash_ID{$sample[6]};
    }
    ' > PANZ_SV_Reformat.vcf

gst_log "Done! Result files:
$(ls -l PANZ_SV_Reformat.vcf PANZ_SV_ID_each_call.tsv)"

# <<<<<<<<<<<<<<<<<<<<<<<<  PANZ_SV_Reformat.vcf <<<<<<<<<<<<<<<<<<<<<<<<
#? 1 CHROM
#? 2 POS
#? 3 ID
#? 4 REF
#? 5 ALT
#? 6 QUAL
#? 7 FILTER
#? 8 INFO: SUPP=6;SUPP_VEC=1110111;SVLEN=-4693;SVTYPE=DEL;SVMETHOD=SURVIVOR1.0.6;CHR2=1;END=2966815;CIPOS=-9,7;CIEND=-6,71;STRANDS=+-
#? 9 FORMAT: 0GT :1PSV:2LN:3DR :4ST:5QV :6TY :7ID :8RAL:9AAL:10CO
#?           ./.:NaN:0 :0,0:--:NaN:NaN:NaN:NAN:NAN:NAN
#? 10 BSpSV
#? 11 asm_REF
#? 12 asm_AMP
#? 13 NRINS
#? 14 gridss
#? 15 manta
#? 16 Delly

# >>>>>>>>>>>>>>>>>>>>>>>> get Bayestyper VCF >>>>>>>>>>>>>>>>>>>>>>>>
gst_log "Getting Bayestyper VCF ..."

# Get Seq only records
cat PANZ_SV_Reformat.vcf | perl -lane '/^#/ && print && next;$F[3]=uc($F[3]);$F[4]=uc($F[4]);print if $F[4] !~ /\W/ and $F[3] !~ /\W/;' >PANZ_SV_bayestyper.vcf
# Bcftools norm -c to not stop with wrong ref records ( caused by SURVIVOR merge )
# bcftools norm --threads 8 -f $ref -o PANZ_SV_bayestyper.vcf -O v -m +any -c s tmp_PANZ_SV_bayestyper.vcf
# <<<<<<<<<<<<<<<<<<<<<<<< get Bayestyper VCF <<<<<<<<<<<<<<<<<<<<<<<<

# >>>>>>>>>>>>>>>>>>>>>>>> get Delly/SVTYPER VCF >>>>>>>>>>>>>>>>>>>>>>>>
gst_log "Getting Delly/Svtyper VCF ..."
cat $delly_header > PANZ_SV_Delly.vcf
cat PANZ_SV_Reformat.vcf | perl -lane '
    /^#/ && next;
    $,="\t";
    if($F[4] !~ /\W/ and $F[3] !~ /\W/){
        $F[3]=~s/^(\w).*$/$1/;
        $svtype="";
        $svtype=$1 if /SVTYPE=(.*?);/;
        $F[4]="<$svtype>";
    }
    print @F[0..7];
' >> PANZ_SV_Delly.vcf
# <<<<<<<<<<<<<<<<<<<<<<<< get Delly/SVTYPER VCF <<<<<<<<<<<<<<<<<<<<<<<<
