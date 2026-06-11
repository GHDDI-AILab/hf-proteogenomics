#!/bin/bash
./cis_pqtl.R && \
./csv2vcf.sh && \
./annotate.sh && \
./tidy.R
