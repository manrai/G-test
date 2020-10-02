# open terminal session in sb.biodatacatalyst

#bcftools view -R ../../workspace/regions.tsv  \
#  NWD998698.freeze5.v1.vcf.gz \
#  | grep "^[^#;]" > ../../workspace/NWD998698_reduced.vcf

# this works from terminal startup:
bcftools view -r chr11:47332517 ../project-files/JHS/NWD100014.freeze5.v1.vcf.gz | grep "^[^#;]"

# coordinates are GChr38 for TopMed Freeze 5
