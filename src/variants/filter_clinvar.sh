bcftools view -R ../../data/clinvar/regions.tsv  \
  ../../data/clinvar/clinvar_8_31_2020.vcf.gz \
  | grep "^[^#;]" > ../../data/clinvar/clinvar_gnomad.vcf
