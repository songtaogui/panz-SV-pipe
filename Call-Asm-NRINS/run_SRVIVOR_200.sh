#BSUB -J M_NRINS
#BSUB -n 1
#BSUB -R rusage[mem=10GB]
#BSUB -R span[hosts=1]
#BSUB -o M_NRINS.out
#BSUB -e M_NRINS.out
#BSUB -q normal

SURVIVOR merge list.txt 200 1 0 0 0 50 NRINS_survivor200.vcf > logs.txt 2>&1
bash filter_SURVIVOR.sh NRINS_survivor200.vcf 3 NRINS >>logs.txt 2>&1

