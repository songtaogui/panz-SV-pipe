in=Merge_SV_SNP_INDEL_no_decoy.vcf.gz
out=Merge_SV_no_decoy_Rand10SNPInDels_Seed12345.vcf.gz

pigz -d -p 8 -c $in | parallel --pipe --jobs 8 -k -q perl -lane '
BEGIN{
    srand(12345);
}
$rdnn=int(rand(10));
if(/^#/ or $F[2] =~ /PZ00001a/){
    print;
}else{
    print if $rdnn == 7;
}' | bgzip -@ 8 -c > $out