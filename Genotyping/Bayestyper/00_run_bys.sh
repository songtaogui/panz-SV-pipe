#!/usr/bin/env bash

#set -euo pipefail

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

usage="
------------------------------------------------------------
Run bayestyper cluster and genotype
------------------------------------------------------------
USAGE:
    bash $(basename $0) [OPTIONs]

OPT:
    1. simple_ist
    2. vcf_file
    3. ref_fasta
    4. decoy_fasta
    5. threads
------------------------------------------------------------
Author: Songtao Gui
E-mail: songtaogui@sina.com
"
if [[ $# -ne 5 ]]; then 
    echo "$usage" >&2
    exit 1
fi

smp=$1
vcf=$2
ref=$3
decoy=$4
cpu=$5

check_files_exists $smp $vcf $ref $decoy

pre=$(basename $smp)
pre=${pre%%_*}
if [ ! -d "${pre}_geno" ];then
    mkdir ${pre}_geno
fi
cd ${pre}_geno &&\
gst_log "Running Cluster for $pre "
if [ -s "Cluster.done" ];then
    gst_log "Skip running, use pre-exist results"
else
    cmd="bayesTyper cluster -v $vcf -s $smp -g $ref -d $decoy -p $cpu"
    gst_log "CMD: $cmd"
    $cmd
    [ $? -ne 0 ] && gst_err "Non-zero exit for cluster" && exit 1
    echo "Successfully finished" >Cluster.done
fi
gst_log "Run genotype for $pre"
if [ -s "Genotype.done" ];then
    gst_log "Skip running, use pre-exist results"
else
    for unit in bayestyper_unit_*
    do
        geno_cmd="bayesTyper genotype -v $unit/variant_clusters.bin -c bayestyper_cluster_data -s $smp -g $ref -d $decoy -o test_bayestyper -z -p $cpu"
        gst_log "CMD:$geno_cmd"
        $geno_cmd
        [ $? -ne 0 ] && gst_err "Non-zero exit for geno $unit" && exit 1
    done
#    [ $? -ne 0 ] && gst_err "Non-zero exit for geno" && exit 1
    echo "Successfully finished" >Genotype.done
fi

gst_log "All Done!"

