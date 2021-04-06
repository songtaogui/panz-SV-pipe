#!/usr/bin/env bash

# set -euo pipefail

function gst_log () {
    local info=$1
    echo -e "\033[36m[$(date +'%y-%m-%d %H:%M')]\033[0m $info" >&2
}
function gst_err () {
    local info=$1
    echo -e "\033[31m\033[7m[ERROR]\033[0m --> $info" >&2
}
function gst_warn () {
    local info=$1
    echo -e "\033[35m[WARNING]\033[0m --> $info" >&2
}
function check_files_exists(){
    local num_related_file=1
    local related_file=""
    for related_file in  "$@"
    do
        if [[ ! -s "$related_file" ]]; then
            echo "$(tput setaf 1)[ERROR] --> No file: $related_file $(tput sgr0)" >&2
            let num_related_file++
        fi
    done
    [ "$num_related_file" -ne 1 ] && exit 1
}
function check_var_numeric () {
    local var_cc=1
    local var_name=""
    local var=""
    for var_name in "$@"; do
        var=$(eval echo "$"$var_name)
        # add ${var#prefix} substitution to trim sign
        case ${var#[-+]} in
            '')
                echo -e "\033[36m[[ERROR]\033[0m --> $var_name is empty: '$var' (tput sgr0)" >&2
                let var_cc++ ;;
            *.*.*)
                echo "\033[36m[[ERROR]\033[0m --> $var_name has more than one decimal point: '$var' (tput sgr0)" >&2
                let var_cc++ ;;
            *[!0-9.]*)
                echo "\033[36m[[ERROR]\033[0m --> $var_name has a non-digit somewhere in it: '$var' (tput sgr0)" >&2
                let var_cc++ ;;
            *) ;;
        esac >&2
    done
    [ "$var_cc" -ne 1 ] && exit 1
}
function check_var_empty () {
    local var_cc=1
    local var_name=""
    local var=""
    for var_name in "$@"; do
        var=$(eval echo "$"$var_name)
        case ${var} in
            '')
                echo "$(tput setaf 1)[ERROR] --> $var_name is empty: '$var' (tput sgr0)" >&2
                let var_cc++ ;;
            *) ;;
        esac >&2
    done
    [ "$var_cc" -ne 1 ] && exit 1
}
# >>>>>>>>>>>>>>>>>>>>>>>> Function to dealing with each genotype >>>>>>>>>>>>>>>>>>>>>>>>
# PZ00001aSV00000629DEL#1#369621#1#369888 svtyper:1/1;Delly:1/1;Manta:1/1;Gridss:1/1;AMPasmSV:1/1
function geno_merge () {
    local geno_pre=$( echo -e "B\tbayestyper\nD\tdelly\nM\tmanta\nS\tsvtyper\nA\tasmsv\nG\tgridss\nN\tnrins" )
    # while read cur_sv cur_geno
    # do
    local cur_sv=$(echo "$1" | cut -f 1 )
    local cur_geno=$(echo "$1" | cut -f 2 )
    check_var_empty cur_sv cur_geno
    cur_geno_tsv=$( csvtk join -tTH -f "2;1" -k <(echo "$geno_pre") <(echo "$cur_geno" | sed 's/;/\n/g;s/:/\t/g') | cut -f 1,3 | perl -lane '$F[1]="./." unless $F[1];$F[2]=1;$F[2]=0 if $F[1] eq "./.";$,="\t";print @F' )
    [ -z "$cur_geno_tsv" ] && gst_err "Wrong geno format for $cur_sv" && exit 1
    # ? Dealing with BDMS
    cur_weighted_BDMS=$( csvtk join -tTH <(echo "$BDMS_weight") <(echo "$cur_geno_tsv") | perl -F"\t" -lane '$,="\t";print $F[0],$F[2],$F[1]*$F[3];' )
    [ -z "$cur_weighted_BDMS" ] && gst_err "Wrong weighted_BDMS for $cur_sv" && exit 1
    # ? total number
    cur_TN_BDMS=$( echo "$cur_weighted_BDMS" | grep -v -P "\./\." | perl -lane 'BEGIN{$sum=0} $sum+=$F[2];END{print $sum}' )
    # ? Ref Alt Hetero Number
    cur_RN_BDMS=$( echo "$cur_weighted_BDMS" | grep -P "0/0" | perl -lane 'BEGIN{$sum=0} $sum+=$F[2];END{print $sum}' )
    cur_HN_BDMS=$( echo "$cur_weighted_BDMS" | grep -P "0/1|1/0" | perl -lane 'BEGIN{$sum=0} $sum+=$F[2];END{print $sum}' )
    cur_AN_BDMS=$( echo "$cur_weighted_BDMS" | grep -P "1/1" | perl -lane 'BEGIN{$sum=0} $sum+=$F[2];END{print $sum}' )
    # ? Ref Alt Hetero Rate
    if [ "$cur_TN_BDMS" -gt "0" ];then
        cur_RR_BDMS=$( echo "scale=2; $cur_RN_BDMS / $cur_TN_BDMS" | bc)
        cur_HR_BDMS=$( echo "scale=2; $cur_HN_BDMS / $cur_TN_BDMS" | bc)
        cur_AR_BDMS=$( echo "scale=2; $cur_AN_BDMS / $cur_TN_BDMS" | bc)
    else
        cur_RR_BDMS=0
        cur_HR_BDMS=0
        cur_AR_BDMS=0
    fi
    # ? chech if valid
    check_var_numeric cur_RR_BDMS cur_AR_BDMS cur_HR_BDMS
    # ? Determine genotype: { Alt_rate > cut-off --> 1/1 }else{ Hetero_rate > cut-off --> 0/1 }else{ Ref_rate > cut-off --> 0/0 }else{ ./. }
    cur_GT_BDMS="./." && cur_ER="0,0,0" && cur_EN="0,0,0"
    [ "$(echo "$cur_RR_BDMS >= $cutoff" | bc)" -eq 1 ] && cur_GT_BDMS="0/0"
    [ "$(echo "$cur_HR_BDMS >= $cutoff" | bc)" -eq 1 ] && cur_GT_BDMS="0/1"
    [ "$(echo "$cur_AR_BDMS >= $cutoff" | bc)" -eq 1 ] && cur_GT_BDMS="1/1"
    cur_ER="$cur_RR_BDMS,$cur_HR_BDMS,$cur_AR_BDMS"
    cur_EN="$cur_RN_BDMS,$cur_HN_BDMS,$cur_AN_BDMS"
    cur_EV=$( echo "$cur_geno_tsv" | head -n 4 | cut -f 2 | perl -pane 's/\n/,/ unless eof' )
    # ? Dealing with AGN
    cur_SP=$( echo "$cur_geno_tsv" | tail -n 3 | cut -f 3 | perl -pane 's/\n/,/ unless eof' )
    cur_SN=$( echo "$cur_SP" | perl -F"," -lane 'BEGIN{use List::Util qw/sum/;} print sum @F' )
    cur_GT_AGN="./."
    [ "$cur_SN" -ge 1 ] && cur_GT_AGN="1/1"
    # ? get final genotype
    cur_GT=$cur_GT_BDMS && cur_GM="BDMS"
    if [[ "$cur_GT_BDMS" == "./." ]];then
        cur_GT=$cur_GT_AGN; cur_GM="AGN"
        [[ "$cur_GT_AGN" == "./." ]] && cur_GT=$cur_GT_AGN && cur_GM="NONE"
    elif [[ "$cur_GT_BDMS" =~ "1" && "$cur_GT_AGN" == "1/1" ]];then
        cur_GT=$cur_GT_BDMS; cur_GM="BOTH"
    elif [[ "$cur_GT_BDMS" == "0/0" && "$cur_GT_AGN" == "1/1" ]];then
        cur_GT="./." && cur_GM="CONF"
    fi
    # [ "$cur_GT_BDMS" == "./." ] && cur_GT=$cur_GT_AGN && cur_GM="AGN"
    # [[ "$cur_GT_BDMS" == "./." && "$cur_GT_AGN" == "./." ]] && cur_GT=$cur_GT_AGN && cur_GM="NONE"
    # # ? BDMS GT and AGN GT conflicts, mark as NA
    # [[ "$cur_GT_BDMS" == "0/0" && "$cur_GT_AGN" == "1/1" ]] && cur_GT="./." && cur_GM="CONF"
    # ? get final outputs: SV\tGT:EV:EN:ER:SP:SN:GM
    echo -e "${cur_sv}\t${cur_GT}:${cur_GM}:${cur_EV}:${cur_EN}:${cur_ER}:${cur_SP}:${cur_SN}"
    # done
}
# ? export functions for parallel
export -f geno_merge check_var_numeric gst_log gst_warn gst_err check_var_empty
# <<<<<<<<<<<<<<<<<<<<<<<< Function to dealing with each genotype <<<<<<<<<<<<<<<<<<<<<<<<

usage="
------------------------------------------------------------
Merging all SV geno evidences to get final genotype
------------------------------------------------------------
USAGE:
    bash $(basename $0) [OPTIONS]

OPTIONS: ([R]:required  [O]:optional)
    -h, --help                show help and exit.
    -v, --vcf        <str>    [R]   Raw SV records vcf file
    -l, --list       <str>    [R]   List of sample ID name map, will use the name and the order of this file 
                                    to generate the final result, so make sure the sample names are unique.
                                        Format (Tab-sep-table):
                                            <ID>    <Sample Name>
                                        Eg.:
                                            Q417    SK
                                            Q441    B73
    --genorate       <0-1>    [O]   Genotype rate cutoff for filtering (Default: 0.1)
    --maf            <0-1>    [O]   Minus allel frequence cutoff for filtering (Default: 0.01)
    -b, --bayestyper <str>    [R]   bayestyper genotype matrix
    -d, --delly      <str>    [R]   delly genotype matrix
    -s, --svtyper    <str>    [R]   svtyper genotype matrix
    -m, --manta      <str>    [R]   manta genotype matrix
    -a, --asmsv      <str>    [R]   AsmSV genotype matrix
    -g, --gridss     <str>    [R]   Gridss genotype matrix
    -n, --nrins      <str>    [R]   NR-INS genotype matrix
    -o, --output     <str>    [O]   Output file name (default: 00_PANZ_merged_geno.matrix)
    -t, --threads    <num>    [O]   set threads (default 2)
    -w, --weight     <str>    [O]   comma seperated integer weight value for setting the genotypes, format:
                                    bayestyper_w,delly_w,manta_w,svtyper_w
                                    (Default: \"2,1,1,1\")
    -c, --cutoff     <0-1>    [O]   Cut-off value to set final genotype (Default 0.5)
------------------------------------------------------------
Author: Songtao Gui
E-mail: songtaogui@sina.com
"
if [[ $# -eq 0 ]]; then
    echo "$usage" >&2
    exit 1
fi

# >>>>>>>>>>>>>>>>>>>>>>>> Parse Options >>>>>>>>>>>>>>>>>>>>>>>>
# Set Default Opt
vcf=
list=
v_bayestyper=
v_delly=
v_svtyper=
v_manta=
v_asmsv=
v_gridss=
v_nrins=
weight="2,1,1,1"
output="00_PANZ_merged_geno.matrix"
export cutoff=0.5
threads=2
# parse args
UNKOWN_ARGS=()
while [[ $# > 0 ]]; do
    case "$1" in
        -h|--help)
            echo "$usage" >&2
            exit 1
        ;;
        -v|--vcf)
            #echo "set argument \"$1\" with value: $2" >&2
            vcf=$2
            shift 2
        ;;
        -l|--list)
            #echo "set argument \"$1\" with value: $2" >&2
            list=$2
            shift 2
        ;;
        --genorate)
            #echo "set argument \"$1\" with value: $2" >&2
            genorate=$2
            shift 2
        ;;
        --maf)
            #echo "set argument \"$1\" with value: $2" >&2
            maf=$2
            shift 2
        ;;
        -b|--bayestyper)
            #echo "set argument \"$1\" with value: $2" >&2
            v_bayestyper=$2
            shift 2
        ;;
        -d|--delly)
            #echo "set argument \"$1\" with value: $2" >&2
            v_delly=$2
            shift 2
        ;;
        -s|--svtyper)
            #echo "set argument \"$1\" with value: $2" >&2
            v_svtyper=$2
            shift 2
        ;;
        -m|--manta)
            #echo "set argument \"$1\" with value: $2" >&2
            v_manta=$2
            shift 2
        ;;
        -a|--asmsv)
            #echo "set argument \"$1\" with value: $2" >&2
            v_asmsv=$2
            shift 2
        ;;
        -g|--gridss)
            #echo "set argument \"$1\" with value: $2" >&2
            v_gridss=$2
            shift 2
        ;;
        -n|--nrins)
            #echo "set argument \"$1\" with value: $2" >&2
            v_nrins=$2
            shift 2
        ;;
        -o|--output)
            #echo "set argument \"$1\" with value: $2" >&2
            output=$2
            shift 2
        ;;
        -w|--weight)
            #echo "set argument \"$1\" with value: $2" >&2
            weight=$2
            shift 2
        ;;
        -c|--cutoff)
            #echo "set argument \"$1\" with value: $2" >&2
            cutoff=$2
            shift 2
        ;;
        -t|--threads)
            #echo "set argument \"$1\" with value: $2" >&2
            threads=$2
            shift 2
        ;;
        *) # unknown flag/switch
            UNKOWN_ARGS+=("$1")
            shift
        ;;
    esac
done
if [ "${#UNKOWN_ARGS[@]}" -gt 0 ];then
    gsterr "[ERROR] --> Wrong options: \"${UNKOWN_ARGS[@]}\"" >&2
    exit 1
fi
unset UNKOWN_ARGS # restore UNKOWN_ARGS params
# Check if required vars are legal
check_var_empty v_bayestyper v_delly v_manta v_asmsv v_gridss v_nrins vcf list weight output
check_var_numeric cutoff threads
check_files_exists $v_bayestyper $v_delly $v_manta $v_asmsv $v_gridss $v_nrins $vcf $list
# <<<<<<<<<<<<<<<<<<<<<<<< Parse Options <<<<<<<<<<<<<<<<<<<<<<<<


# >>>>>>>>>>>>>>>>>>>>>>>> Main: parse total matrix and get merged geno >>>>>>>>>>>>>>>>>>>>>>>>

# >>>>>>>>>>>>>>>>>>>>>>>> Pre-dealing >>>>>>>>>>>>>>>>>>>>>>>>
# ? generate matrix based on list and vcf
mkdir -p ${output}_geno_matrix
# ? set ouputs
m_bayestyper=${output}_geno_matrix/bayestyper.geno_matrix
m_delly=${output}_geno_matrix/delly.geno_matrix
m_manta=${output}_geno_matrix/manta.geno_matrix
m_svtyper=${output}_geno_matrix/svtyper.geno_matrix
m_asmsv=${output}_geno_matrix/asmsv.geno_matrix
m_gridss=${output}_geno_matrix/gridss.geno_matrix
m_nrins=${output}_geno_matrix/nrins.geno_matrix
# ? generate matrix
# get id list and sample head
id_csv=$(cut -f 1 $list | csvtk transpose)
sample_header_tsv=$(cut -f 2 $list | csvtk transpose -tT )
check_var_empty id_csv sample_header_tsv
for cur_method in bayestyper delly manta svtyper asmsv gridss nrins
do
    gst_log "Generating genotype matrix for $cur_method ..."
    cur_vcf_file=$(eval echo "$"v_$cur_method)
    cur_mtx_file=$(eval echo "$"m_$cur_method)
    check_files_exists $cur_vcf_file
    # ? get matrix for cur_method
    gst_log "Getting geno_matrix for $cur_method ..."
    if [ -s "$cur_mtx_file" ];then
        gst_warn "Using pre-exist file: $cur_mtx_file"
    else
        echo -e "ID\t$sample_header_tsv" >$cur_mtx_file
        bcftools query -f "%ID[\t$cur_method:%GT]\n" -s "$id_csv" $cur_vcf_file >>$cur_mtx_file
        if [ $? -ne 0 ];then gst_err "bcftools query failed for $cur_method: Non-zero exit"; rm -f $cur_mtx_file; exit 1;fi
        gst_log "Generated $cur_method geno_matrix: $cur_mtx_file"
    fi
done
# ? dealing with weight
check_weight=$( echo "$weight" | grep -P "^\d+,\d+,\d+,\d+$" )
[ -z "$check_weight" ] && gst_log "Wrong weight format: $weight" && exit 1
export BDMS_weight=$( echo -e "B,D,M,S\n$weight" | csvtk transpose -TH )
echo -e "Using weights:\n$BDMS_weight" >&2
# ? dealing with vcf and header
gst_log "Getting vcf header to ${output}_vcfheader.tmp ..."
head -n 1000 $vcf | grep "^##" | grep -v "FORMAT=<ID=GT," > ${output}_vcfheader.tmp
echo  "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##FORMAT=<ID=EV,Number=1,Type=String,Description=\"Genotype of BDMS (bayestyper,delly,manta,svtyper),format: 0/0,./.,1/1,0/1 for BDMS;\">
##FORMAT=<ID=EN,Number=1,Type=Integer,Description=\"Weighted number of BDMS evidences,format:ref,hetero,alt\">
##FORMAT=<ID=ER,Number=1,Type=Float,Description=\"Weighted rate of the BDSM evidence genotype,format:ref,hetero,alt\">
##FORMAT=<ID=SP,Number=1,Type=String,Description=\"Genotype of AGN (assemblytics,gridss,non_ref_insertions) supports, format: 0,1,1 \">
##FORMAT=<ID=SN,Number=1,Type=Integer,Description=\"Number of AGN evidences\">
##FORMAT=<ID=GM,Number=1,Type=String,Description=\"Genotype Marker: AGN/BDMS/CONF/NONE/BOTH. AGN for AGN supports; BDMS for BDMS supports; CONF for conflicts between AGN and BDMS; BOTH for both AGN and BDMS supports; NONE for No evidences\">" >> ${output}_vcfheader.tmp
# ? dealing with raw vcf records
gst_log "Getting no header pre-vcf file to ${output}_vcfNoheader.tmp ..."
if [ -s "${output}_vcfNoheader.tmp" ];then
    gst_warn "Using pre-exist file: ${output}_vcfNoheader.tmp"
else
    cat $vcf | parallel --pipe -q -k --jobs $threads perl -F"\t" -lane '/^##/ && next;$F[8]="GT:GM:EV:EN:ER:SP:SN" if $F[0] !~ /^#CH/;$,="\t";print @F[0..8];' > ${output}_vcfNoheader.tmp 2>/dev/null
    if [ $? -ne 0 ];then gst_err "parallel get raw vcf: Non-zero exit"; exit 1;fi
fi
# ? merge all matrix
gst_log "Getting merged matrix to ${output}_merged_matrix.tsv ..."
if [ -s "${output}_merged_matrix.tsv" ];then
    gst_log "Using pre-exit file: ${output}_merged_matrix.tsv"
else
    csvtk concat -tT -j $threads $m_bayestyper $m_delly $m_manta $m_svtyper $m_asmsv $m_gridss $m_nrins -o ${output}_merged_matrix.tsv
    if [ $? -ne 0 ];then gst_err "csvtk concat failed: Non-zero exit"; exit 1;fi
fi
# ? sort matrix
gst_log "Getting sorted matrix to ${output}_merged_matrix_sort.tsv ..."
if [ -s "${output}_merged_matrix_sort.tsv" ];then
    gst_warn "Using pre-exist file ${output}_merged_matrix_sort.tsv"
else
    csvtk sort -tT -j $threads ${output}_merged_matrix.tsv -o ${output}_merged_matrix_sort.tsv
    if [ $? -ne 0 ];then gst_err "csvtk sort failed: Non-zero exit"; exit 1;fi
fi
# <<<<<<<<<<<<<<<<<<<<<<<< Pre-dealing <<<<<<<<<<<<<<<<<<<<<<<<

gst_log "Merge genotype for each samples ..."
if [ -s "${output}_all_sample_mgeno.tmp" ];then
    gst_warn "Using pre-exist file ${output}_all_sample_mgeno.tmp"
else
    mkdir -p ${output}_each_sample_tmp
    for cur_sample in $(head -n 1 ${output}_merged_matrix_sort.tsv | cut -f 2- | csvtk transpose -tTH)
    do
        gst_log "Dealing with $cur_sample ..."
        echo -e "ID\t$cur_sample" > ${output}_each_sample_tmp/${cur_sample}.mgeno.tmp
        if [ ! -s "{output}_each_sample_tmp/${cur_sample}.mgeno.tmp" ];then
            [ -s "${output}_each_sample_tmp/${cur_sample}.collapse.tmp" ] || csvtk cut -tT -f "ID,$cur_sample" ${output}_merged_matrix_sort.tsv | sed '1d' | csvtk collapse -tTH -j $threads -f 1 -v 2 -o ${output}_each_sample_tmp/${cur_sample}.collapse.tmp
            if [ $? -ne 0 ];then gst_err "collapse failed for $cur_sample : Non-zero exit"; exit 1;fi
            parallel -k --jobs $threads geno_merge :::: ${output}_each_sample_tmp/${cur_sample}.collapse.tmp >> ${output}_each_sample_tmp/${cur_sample}.mgeno.tmp
            if [ $? -ne 0 ];then gst_err "merge geno for $cur_sample failed: Non-zero exit"; exit 1;fi
        fi
    done
fi
# ? merge all mgeno
gst_log "Merging samples ..."
[ -s "${output}_all_sample_mgeno.tmp" ] || csvtk join -tTH -j $threads -f 1 ${output}_each_sample_tmp/*.mgeno.tmp -o ${output}_all_sample_mgeno.tmp
if [ $? -ne 0 ];then gst_err "csvtk join mgeno failed: Non-zero exit"; exit 1;fi
# ? merge vcf and mgeno and get final vcf results
gst_log "Generating final VCF ..."
[ -s "${output}_all_sample_mgeno.nh.vcf" ] || csvtk join -tTH -C "$" -j $threads -f "3;1" ${output}_vcfNoheader.tmp ${output}_all_sample_mgeno.tmp -o ${output}_all_sample_mgeno.nh.vcf
if [ $? -ne 0 ];then gst_err "csvtk join vcf failed: Non-zero exit"; exit 1;fi
# ? add header
cat ${output}_vcfheader.tmp ${output}_all_sample_mgeno.nh.vcf > ${output}_all_sample_mgeno.vcf
# ? sort bgzip index
gst_log "Sorting VCF ..."
[ -s "${output}_all_sample_mgeno_sort.vcf.gz" ] || bcftools sort ${output}_all_sample_mgeno.vcf -O z -o ${output}_all_sample_mgeno_sort.vcf.gz
if [ $? -ne 0 ];then gst_err "bcftools sort failed: Non-zero exit"; exit 1;fi
gst_log "Indexing VCF ..."
[ -s "${output}_all_sample_mgeno_sort.vcf.gz" ] || bcftools index ${output}_all_sample_mgeno_sort.vcf.gz
if [ $? -ne 0 ];then gst_err "bcftools index failed: Non-zero exit"; exit 1;fi
# ? clean tmp files
# rm -f ${output}_*.tmp ${output}_*.vcf
# rm -rf ${output}_each_sample_tmp
gst_log "All Done !
>>>> Final result:
>>>> ${output}_all_sample_mgeno_sort.vcf.gz
"
# <<<<<<<<<<<<<<<<<<<<<<<< Main: parse total matrix and get merged geno <<<<<<<<<<<<<<<<<<<<<<<<





# # >>>>>>>>>>>>>>>>>>>>>>>> run loop with multi threads >>>>>>>>>>>>>>>>>>>>>>>>
# # >>>>>>>>>>>>>>>>>>>>>>>> gst progress bar >>>>>>>>>>>>>>>>>>>>>>>>
# # gst progress bar (outside_loop)
# # usage: gst_bar_outloop <num_times_of_loop (> 100)>
# function gst_bar_outloop () {
#     progress_bar_tnum=$1
#     echo -n "Progress:( Total records: $progress_bar_tnum )" >&2
#     if [ "$progress_bar_tnum" -gt 100 ];then
#         rm -f .progress_bar_calc_tmp .progress_bar_count_tmp
#         echo "
# 0%   10   20   30   40   50   60   70   80   90   100%
# |----|----|----|----|----|----|----|----|----|----|" >&2
#         # creat tmp files for counting (use tmp files instead of VARs for FIFO-loop)
#         echo "start" > .progress_bar_count_tmp
#         touch .progress_bar_calc_tmp
#     else
#         echo "Loops too short for showing a progress bar" >&2
#     fi
# }

# # gst progress bar (inside_loop)
# # usage(put it in a loop): gst_bar_inloop
# function gst_bar_inloop () {
#     if [ "$progress_bar_tnum" -gt 100 ];then
#         progress_bar_count=$(echo "$(cat .progress_bar_count_tmp | wc -l) * 100 / $progress_bar_tnum" | bc)
#         progress_bar_calc=$(cat .progress_bar_calc_tmp | wc -l)
#         echo "$progress_bar_count" >> .progress_bar_count_tmp
#         if [[ "$progress_bar_count" -eq "$progress_bar_calc" ]];then
#             echo -e "${progress_bar_calc}_1\nprogress_bar_calc_2" >> .progress_bar_calc_tmp
#             echo -n "*" >&2;
#         fi
#     fi
# }
# # <<<<<<<<<<<<<<<<<<<<<<<< gst progress bar <<<<<<<<<<<<<<<<<<<<<<<<
# #progress bar outside loop,
# gst_bar_outloop $(head -n 1 ${output}_merged_matrix_sort.tsv | cut -f 2- | csvtk transpose -tTH | wc -l)
# # set threads
# loop_threads=$threads
# # open fd
# [ -e /tmp/fd1 ] || mkfifo /tmp/fd1
# exec 7<>/tmp/fd1
# rm -rf /tmp/fd1
# for i in $(seq 0 $loop_threads);do
#     echo
# done >&7
# # loop start
# for cur_sample in $(head -n 1 ${output}_merged_matrix_sort.tsv | cut -f 2- | csvtk transpose -tTH)
# do
# read -u7
# {
#     #gst_log "Dealing with $cur_sample ..."
#     gst_bar_inloop
#     echo -e "ID\t$cur_sample" > ${output}_each_sample_tmp/${cur_sample}.mgeno.tmp
#     if [ ! -s "{output}_each_sample_tmp/${cur_sample}.mgeno.tmp" ];then
#         [ -s "${output}_each_sample_tmp/${cur_sample}.collapse.tmp" ] || csvtk cut -tT -f "ID,$cur_sample" ${output}_merged_matrix_sort.tsv | sed '1d' | csvtk collapse -tTH -j $threads -f 1 -v 2 -o ${output}_each_sample_tmp/${cur_sample}.collapse.tmp
#         if [ $? -ne 0 ];then gst_err "collapse failed for $cur_sample : Non-zero exit"; exit 1;fi
#         parallel --jobs 1 geno_merge :::: ${output}_each_sample_tmp/${cur_sample}.collapse.tmp >> ${output}_each_sample_tmp/${cur_sample}.mgeno.tmp
#         if [ $? -ne 0 ];then gst_err "merge geno for $cur_sample failed: Non-zero exit"; exit 1;fi
#     fi
# } &
# done | cat
# wait
# exec 7>&- # close fd
# exec 7<&- # close fd
# # <<<<<<<<<<<<<<<<<<<<<<<< run loop with multi threads <<<<<<<<<<<<<<<<<<<<<<<<