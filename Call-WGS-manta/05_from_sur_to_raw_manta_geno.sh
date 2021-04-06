#!/usr/bin/env bash

# set -o xtrace
# set -o errexit
set -o nounset
set -o pipefail

# >>>>>>>>>>>>>>>>>>>>>>>> Common functions >>>>>>>>>>>>>>>>>>>>>>>>
gst_log () {
    local info=$1
    echo -e "\033[36m[$(date +'%y-%m-%d %H:%M')]\033[0m $info" >&2
}

gst_err () {
    local info=$1
    echo -e "\033[31m\033[7m[ERROR]\033[0m --> $info" >&2
}

gst_warn () {
    local info=$1
    echo -e "\033[35m[WARNING]\033[0m --> $info" >&2
}

check_files_exists(){
    local num_related_file=1
    local related_file=""
    for related_file in  "$@"
    do
        if [[ ! -s "$related_file" ]]; then
            echo -e "\033[31m\033[7m[ERROR]\033[0m --> No file: $related_file $" >&2
            let num_related_file++
        fi
    done
    [ "$num_related_file" -ne 1 ] && exit 1
}

check_abs_path() {
    local var_cc=1
    local check_file=""
    for check_file in "$@";do
        if [[ "${check_file:0:1}" != "/" ]]; then
            echo -e "\033[31m\033[7m[ERROR]\033[0m --> $check_file was not an ABSOLUTE path." >&2
            let var_cc++
        fi
    done
    [ "$var_cc" -ne 1 ] && exit 1
}

check_sftw_path(){
    local num_tp_program=1
    local tp_program=""
    for tp_program in "$@"
    do
        if ! which $tp_program >/dev/null 2>&1 ; then
            echo -e "\033[31m\033[7m[ERROR]\033[0m --> Program not in PATH: $tp_program $" >&2
            let num_tp_program++
        fi
    done
    [ "$num_tp_program" -ne 1 ] && exit 1
}

check_var_empty () {
    local var_cc=1
    local var_name=""
    local var=""
    for var_name in "$@"; do
        var=$(eval echo "$"$var_name)
        case ${var} in
            '')
                echo -e "\033[31m\033[7m[ERROR]\033[0m --> $var_name is empty: '$var' " >&2
                let var_cc++ ;;
            *) ;;
        esac >&2
    done
    [ "$var_cc" -ne 1 ] && exit 1
}

check_var_numeric () {
    local var_cc=1
    local var_name=""
    local var=""
    for var_name in "$@"; do
        var=$(eval echo "$"$var_name)
        # add ${var#prefix} substitution to trim sign
        case ${var#[-+]} in
            '')
                echo -e "\033[31m\033[7m[ERROR]\033[0m --> $var_name is empty: '$var' " >&2
                let var_cc++ ;;
            *.*.*)
                echo -e "\033[31m\033[7m[ERROR]\033[0m --> $var_name has more than one decimal point: '$var' " >&2
                let var_cc++ ;;
            *[!0-9.]*)
                echo -e "\033[31m\033[7m[ERROR]\033[0m --> $var_name has a non-digit somewhere in it: '$var' " >&2
                let var_cc++ ;;
            *) ;;
        esac >&2
    done
    [ "$var_cc" -ne 1 ] && exit 1
}

check_suffix () {
    check_suffix_file=$( basename $1 )
    check_suffix=$2
    # add x incase file has no suffix
    if [[ "${check_suffix_file##*.}"x != "$check_suffix"x ]];then
        echo "[ERROR] --> $check_suffix_file should have suffix: '$check_suffix'." >&2
        exit 1
    fi
}

export -f gst_log gst_warn gst_err check_var_empty check_var_numeric check_sftw_path check_suffix check_files_exists check_abs_path
# <<<<<<<<<<<<<<<<<<<<<<<< Common functions <<<<<<<<<<<<<<<<<<<<<<<<

usage=$(
cat <<EOF
------------------------------------------------------------
Convert Genotype in SURVIVOR merged manta result back into there raw genotype format
------------------------------------------------------------
Dependence: csvtk parallel
------------------------------------------------------------
USAGE:
    bash $(basename $0) <survivor_merged_vcf> <raw_each_sample_vcf_list> <output> <threads>
list format:
    <sample_ID>   <File_Path>
------------------------------------------------------------
Author: Songtao Gui
E-mail: songtaogui@sina.com
EOF
)
if [[ $# -ne 4 ]]; then 
    echo "$usage" >&2
    exit 1
fi

export vcf=$1
export list=$2
export out=$3
export threads=$4

check_files_exists $vcf $list
check_var_numeric threads
# ? set temp dir path
export temp_dir=${vcf##*/}.temp_dir
mkdir -p $temp_dir

# USAGE:$0 sample_ID
fmt_each_sample () {
    local sample=$1
    local ref=$(grep -w "^$sample" $list | cut -f 2)
    check_files_exists $ref
    if [[ ! -s "$temp_dir/$sample.manta" ]];then
        # ? get $sample survivor geno
        if [[ ! -s "$temp_dir/$sample.survivor" ]];then
            grep -v "^##" $vcf | csvtk cut -tT -f "$sample" -C "@" | perl -F":" -lane '
                if(/Manta/){
                    $id=join(":",@F[7..14]);
                }else{
                    $id=$F[0];
                }
                print $id;
            ' > $temp_dir/$sample.survivor
            if [ $? -ne 0 ];then gst_err "get $sample.survivor failed: Non-zero exit"; exit 1;fi
        fi
        # ? convert to raw manta geno
        csvtk join -tTH -k -f 1 $temp_dir/$sample.survivor <(csvtk cut -tTH -f 3,9,10 $ref) | perl -lane '
            $geno="./.:NA:0:0,0,0:0,0:0,0";
            $geno=$F[2] if $F[2];
            $geno=$F[2].":0,0" if $F[2] && $F[1] !~ /SR/;
            $geno=$F[0] if $.==1;print $geno;
        ' > $temp_dir/$sample.manta
        check_files_exists $temp_dir/$sample.manta
        if [ $? -ne 0 ];then gst_err "get $sample.manta failed: Non-zero exit";rm -f $temp_dir/$sample.manta; exit 1;fi
    fi
}
export -f fmt_each_sample

gst_log "Dealing with $vcf ..."
# ? get header and non-geno vcf part of survivor vcf
gst_log "Get header and non-geno part ..."
# header
perl -lne 'if(/^##/){print;}else{last;}' $vcf | grep -v "##FORMAT=" > $temp_dir/01_header.txt &&\
cat <<EOF >> $temp_dir/01_header.txt
##FORMAT=<ID=FT,Number=1,Type=String,Description="Sample filter, 'PASS' indicates that all filters have passed for this sample">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
##FORMAT=<ID=PR,Number=.,Type=Integer,Description="Spanning paired-read support for the ref and alt alleles in the order listed">
EOF
if [ $? -ne 0 ];then gst_err "get header failed: Non-zero exit"; exit 1;fi
# non-geno-part
if [[ ! -s "$temp_dir/02_no_geno_part.tsv" ]];then
    grep -v "^##" $vcf | cut -f 1-9 | perl -F"\t" -lane 'BEGIN{$,="\t";}$F[-1]="GT:FT:GQ:PL:PR:SR" unless $.==1;print @F;' > $temp_dir/02_no_geno_part.tsv
    if [ $? -ne 0 ];then gst_err "get no_geno_part failed: Non-zero exit";rm -f $temp_dir/02_no_geno_part.tsv; exit 1;fi
fi
# sample list
if [[ ! -s "$temp_dir/00_sample_list.txt" ]];then
    grep -m 1 "^#CHR" $vcf | cut -f 10- | csvtk transpose -tTH -o $temp_dir/00_sample_list.txt
fi
# ? fmt each sample
gst_log "Formatting each sample ..."
parallel --bar -j $threads -q -k fmt_each_sample :::: $temp_dir/00_sample_list.txt
if [ $? -ne 0 ];then gst_err "format each sample failed: Non-zero exit"; exit 1;fi

# ? merge all sample
gst_log "Merging all samples ..."
paste $temp_dir/02_no_geno_part.tsv $temp_dir/*.manta > $temp_dir/03_no_header_manta_geno.tsv
if [ $? -ne 0 ];then gst_err "merge all samples manta geno failed: Non-zero exit"; exit 1;fi
gst_log "get final result ..."
cat $temp_dir/01_header.txt $temp_dir/03_no_header_manta_geno.tsv > $out
if [ $? -ne 0 ];then gst_err "get final out failed: Non-zero exit";rm -f $out; exit 1;fi
# clean temp files
rm -rf $temp_dir
gst_log "All Done !"


