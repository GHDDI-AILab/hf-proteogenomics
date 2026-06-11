#!/bin/bash
outdir=../results
header=$outdir/vcf_header.txt
body=$outdir/vcf_body.txt

cat > $header << 'EOF'
##fileformat=VCFv4.2
##INFO=<ID=BETA,Number=1,Type=Float,Description="Effect size">
##INFO=<ID=SE,Number=1,Type=Float,Description="Standard error">
##INFO=<ID=LOG10P,Number=1,Type=Float,Description="-log10(p-value)">
##INFO=<ID=PVAL,Number=1,Type=Float,Description="P-value">
##INFO=<ID=AF,Number=1,Type=Float,Description="Allele frequency">
##INFO=<ID=EXPOSURE,Number=1,Type=String,Description="Exposure">
##INFO=<ID=SAMPLESIZE,Number=1,Type=Integer,Description="Sample size">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
EOF

for i in {1..22} X Y MT; do
  echo "##contig=<ID=$i>" >> $header
done

echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE" >> $header

awk -F',' '
BEGIN {OFS="\t"}
NR==1 {next}
{
  chrom = $1
  pos = $2
  id = "."
  ref = $4
  alt = $5
  af = $6
  beta = $7
  se = $8
  log10p = $9
  pval = $10
  exposure = $11
  samplesize = $12
  info = sprintf("BETA=%s;SE=%s;LOG10P=%s;PVAL=%s;AF=%s;EXPOSURE=%s;SAMPLESIZE=%s",
                 beta, se, log10p, pval, af, exposure, samplesize)
  print chrom, pos, id, ref, alt, ".", ".", info, "GT", "./."
}' \
  $outdir/ukb-ppp.instruments.csv > $body

cat $header $body > $outdir/ukb-ppp.instruments.vcf
bgzip -f $outdir/ukb-ppp.instruments.vcf
tabix -p vcf $outdir/ukb-ppp.instruments.vcf.gz
