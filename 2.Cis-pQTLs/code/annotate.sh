#!/bin/bash
outdir=../results

bcftools norm \
	-c s \
	-f ../data/Ensembl/release-114/GRCh38.fa \
	$outdir/ukb-ppp.instruments.vcf.gz \
	-Oz -o $outdir/ukb-ppp.instruments.norm.vcf.gz && \
tabix -p vcf $outdir/ukb-ppp.instruments.norm.vcf.gz && \
\
bcftools annotate \
	-a ../data/dbsnp/dbsnp.v157.norm.vcf.gz \
	-c ID \
	$outdir/ukb-ppp.instruments.norm.vcf.gz \
	-Oz -o $outdir/instruments.vcf.gz && \
\
echo \
	'chr.exposure,pos.exposure,SNP,other_allele.exposure,effect_allele.exposure,eaf.exposure,beta.exposure,se.exposure,log10p.exposure,pval.exposure,Protein,samplesize.exposure' \
	> $outdir/instruments.new.csv && \
bcftools query \
	-f '%CHROM,%POS,%ID,%REF,%ALT,%INFO/AF,%INFO/BETA,%INFO/SE,%INFO/LOG10P,%INFO/PVAL,%INFO/EXPOSURE,%INFO/SAMPLESIZE\n' \
	$outdir/instruments.vcf.gz >> $outdir/instruments.new.csv
