# run modified assemblytics and SURVIVOR

pre=$1
paffile=$2
assemblytics_dir=$SFTW/Assemblytics/
delta=$PWD/delta/${pre}.delta

usage=$(
cat <<EOF
------------------------------------------------------------
run modified assemblytics and SURVIVOR
------------------------------------------------------------
Dependence:
    Assemblytics: modified version that support SVs up to 1Mb
    [SURVIVOR](https://github.com/fritzsedlazeck/SURVIVOR)
	[paf2delta.py](https://gist.github.com/malonge/d347a0061eab77d9ebeb9181b7a9ff31)
------------------------------------------------------------
USAGE:
	bash $(basename $0) <pre> <paf.gz>
	pre: the output prefix, usually be the name of the individual
	paf.gz: path to the WGA result between the query genome and the ref genome, in paf format
NOTE:
	Before running, make sure the "assemblytics_dir" in line 5 of this script was correctly set to the right path of your local Assemblytics installation
------------------------------------------------------------
Author: Songtao Gui
E-mail: songtaogui@sina.com
EOF
)
if [[ $# -ne 2 ]]; then 
	echo "$usage" >&2
	exit 1
fi


check_sftw_path(){
	local num_tp_program=1
	local tp_program=""
	for tp_program in "$@"
	do
		if ! which $tp_program >/dev/null 2>&1 ; then
			echo -e "\033[31m\033[7m[ERROR]\033[0m --> Program not in PATH: $tp_program " >&2
			let num_tp_program++
		fi
	done
	[ "$num_tp_program" -ne 1 ] && exit 1
}

check_sftw_path SURVIVOR Assemblytics paf2delta.py

if [ ! -s $delta ];then
	echo "generating $delta file..." >&2
	if [ -s $paffile ];then
		echo "Found $paffile ..." >&2
	else
		echo "No file : $paffile " >&2
		exit 1
	fi
	paf2delta.py $paffile
	if [[ $? -ne 0 ]] ; then 
		echo "[ERROR] --> paf2: non-zero exit." >&2
		rm -f ${paffile}.delta
		exit 1
	fi
	echo "Done! Moving output ${paffile}.delta to $delta ..." >&2
	mv ${paffile}.delta $delta
else
	echo "Found $delta ..." >&2
fi

if [ ! -s $delta ];then
	echo "No file: $delta !" >&2
	exit 1
fi
mkdir -p ${pre}_asmSV
cd ${pre}_asmSV

echo "# ----> Assemblytics ..." >&2
if [ -s finished.txt ];then
	echo "Skip Assemblytics. Already done." >&2
else
	# running Assemblytics with unique_length=1K (max length = 1M)
	Assemblytics $delta ${pre}_asm 1000 $assemblytics_dir
	if [[ $? -ne 0 ]] ; then echo "[ERROR] --> Assemblytics: non-zero exit." >&2;exit $?;fi
	if [ -s ${pre}_asm.Assemblytics_structural_variants.bed ];then	
		echo "# ----> Assemblytics finished" > finished.txt	
	fi
fi

# running SURVIVOR,filter out SVs < 50 bp
echo "# ----> convert to VCF using SURVIVOR and filter out < 50 bp rcds ..." >&2
if [ -s ${pre}_asm.Assemblytics_structural_variants.bed.gz ];then
	gzip -d ${pre}_asm.Assemblytics_structural_variants.bed.gz
fi

SURVIVOR convertAssemblytics ${pre}_asm.Assemblytics_structural_variants.bed 50 ${pre}_asmSV_un1k_m1M_l50.vcf

if [[ $? -ne 0 ]] ; then echo "[ERROR] --> SURVIVOR: non-zero exit." >&2;exit $?;fi

# # gzip files
# echo "# ----> gzipping files ..."
# bgzip ${pre}_asmSV_un1k_m1M_l50.vcf
# if [[ $? -ne 0 ]] ; then echo "[ERROR] --> bgzip: non-zero exit.";exit $?;fi
# find . -size +1000000c | xargs -I {} gzip {}
# if [[ $? -ne 0 ]] ; then echo "[ERROR] --> find & gzip: non-zero exit.";exit $?;fi

echo "# ----> All Done !" >&2
