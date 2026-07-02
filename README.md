# hf-proteogenomics
Proteogenomic analysis of heart failure to identify potential therapeutic targets

## Directory structure
```
.
├── 1.Cohort-study/                 # Cox proportional-hazards survival analysis
│   ├── code/
│   │   └── cox_analysis.R
│   ├── data/
│   │   ├── _prepare_data.R
│   │   ├── diagnoses_of_heart_failure.csv.gz
│   │   ├── pheno_500k_raw.csv.gz
│   │   ├── olink0_data.csv.gz
│   │   ├── mydata.csv.gz
│   │   └── symbol.list
│   └── results/
│       ├── cox_model1_results.csv.gz
│       ├── cox_model2_results.csv.gz
│       └── cox_model3_results.csv.gz
│
├── 2.Cis-pQTLs/                    # Cis-pQTL instrument selection
│   ├── code/
│   │   ├── get_instruments.sh
│   │   ├── cis_pqtl.R
│   │   ├── csv2vcf.sh
│   │   ├── annotate.sh
│   │   └── tidy.R
│   ├── data/
│   │   ├── dbsnp/                   # dbSNP VCF files
│   │   ├── Ensembl/release-114/     # Ensembl annotation
│   │   ├── LD/                      # LD reference panels (bin/, new_1kg_ref_GRCh37/)
│   │   └── UKB-PPP/                 # per-protein UKB-PPP pQTL summary stats
│   └── results/
│       ├── ukb-ppp.instruments.csv
│       ├── ukb-ppp.instruments.vcf.gz(.tbi)
│       ├── ukb-ppp.instruments.norm.vcf.gz(.tbi)
│       ├── instruments.vcf.gz
│       ├── instruments.new.csv
│       └── instruments.csv
│
└── 3.MR-and-coloc/                 # Mendelian randomization & colocalization
    ├── code/
    │   └── mr_and_coloc.R
    ├── data/
    │   └── instruments.csv
    └── results/
        ├── <outcome>.outcome_data.csv
        ├── <outcome>.harmonised_data.csv
        ├── <outcome>.mr_results.csv / .mr_results.fdrtool.csv
        ├── <combination>.meta_analysis[.ivw|.wald].csv
        └── colocalization.{ebi-a-GCST009541,finn-b-I9_HEARTFAIL,ukb-d-I9_HEARTFAIL}.csv
```

