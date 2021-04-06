bcftools convert -O b -o all_delly_geno_merge.bcf $WORK/AMP_MEM_SV_TEs/delly/MKdup/delly_geno_merge/all_delly_merged_geno.vcf.gz
delly filter -f germline -o all_delly_geno_merge_filter.bcf -a 0.2 -r 0.75 -q 15 -e 0.8 -u 1.2 all_delly_geno_merge.bcf
