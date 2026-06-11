#!/usr/bin/env Rscript
library(data.table)

outdir = "../results"
a = file.path(outdir, "instruments.new.csv") |> fread(na.strings = "")

given = a[SNP != ".", unique(Protein)]
a[, log10p.exposure := NULL]

a[SNP != ".", ][order(nchar(chr.exposure), chr.exposure, nchar(pos.exposure), pos.exposure)] |>
  fwrite(file.path(outdir, "instruments.csv"))
