#!/usr/bin/env Rscript
library(data.table)
library(biomaRt)

load_summary_data = function(path) {
  structure(
    data.table::fread(path),
    class = c("SummaryData", "data.table", "data.frame")
  )
}

get_gene_info = function(genename, GRCh = NULL) {
  GRCh = if (length(GRCh) == 1 && GRCh == 37) GRCh else NULL
  biomaRt::getBM(
    attributes = c(
      "hgnc_symbol",
      "ensembl_gene_id",
      "chromosome_name",
      "start_position",
      "end_position",
      "strand"
    ),
    filters = "hgnc_symbol",
    values = genename,
    mart = biomaRt::useEnsembl(
      biomart = "genes",
      dataset = "hsapiens_gene_ensembl",
      GRCh = GRCh
    )
  )
}

ld_based_clump = function(
  data,
  clump_kb = 10000,
  clump_r2 = 0.1,
  clump_p = 1,
  pop = "EUR",
  bfile = "../data/LD/new_1kg_ref_GRCh37/EUR",
  plink_bin = "../data/LD/bin/plink"
) {
  nowtime = format(Sys.time(), "%Y%m%d.%H%M%S")
  fn = sprintf("/tmp/ld_clump_%s.%06X.txt", nowtime, floor(runif(1) * 16^6))
  write.table(data.frame(SNP=data[["rsid"]], P=data[["pval"]]), file=fn, row.names=FALSE, col.names=TRUE, quote=FALSE)
  
  command = paste(
    plink_bin,
    " --bfile ", bfile,
    " --clump ", fn,
    " --clump-p1 ", clump_p, 
    " --clump-r2 ", clump_r2, 
    " --clump-kb ", clump_kb, 
    " --out ", fn
  )
  system(command)
  res = read.table(paste0(fn, ".clumped"), header = TRUE)
  unlink(paste(fn, "*", sep=""))
  y <- subset(data, !data[["rsid"]] %in% res[["SNP"]])
  if(nrow(y) > 0) {
    message("Removing ", length(y[["rsid"]]), " of ",
	    nrow(data), " variants due to LD with other variants or absence from LD reference panel")
  }
  return(subset(data, data[["rsid"]] %in% res[["SNP"]]))
}

cis_pqtl = function(
  summary_data,
  genename,
  p_threshold = 5e-8,
  cis_window = 1e6
) {
  if (class(summary_data)[1] != "SummaryData") {
    stop("summary_data must be of class 'SummaryData'")
  }
  if (length(genename) != 1) {
    stop("genename must be a single gene symbol")
  }
  gene_info = info[hgnc_symbol == genename & !grepl("^H", chromosome_name), ]
  chr = gene_info$chromosome_name
  start = max(1, gene_info$start_position - cis_window)
  end = gene_info$end_position + cis_window

  summary_data[
    CHROM == chr & GENPOS >= start & GENPOS <= end & LOG10P > -log10(p_threshold),
    list(
      CHROM, GENPOS, rsid = sub(":imp:v.$", "", ID), ALLELE0, ALLELE1, A1FREQ,
      BETA, SE, LOG10P, pval = 10^(-LOG10P), exposure = genename, SAMPLESIZE = N
    )
  ][
    pval < .Machine$double.xmin, pval := .Machine$double.xmin
  ]
}

if (sys.nframe() == 0L) {
  datadir = "../data"
  outdir = "../results"
  dir.create(outdir,  recursive = T)
  
  get_infile = function(genename) {
    gene_info = info[hgnc_symbol == genename & !grepl("^H", chromosome_name), ]
    Sys.glob(sprintf(
      "%s/UKB-PPP/%s_*/discovery_chr%s_%s*.gz",
      datadir, genename, gene_info$chromosome_name, genename
    ))
  }

  get_independent_cis_pqtl = function(genename) {
    print(genename)
    print(get_infile(genename))
    get_infile(genename) |>
      load_summary_data() |>
      cis_pqtl(genename, p_threshold = 1e-5, cis_window = 1e6) |>
      ld_based_clump(clump_kb = 10000, clump_r2 = 0.1)
  }

  genes = file.path(datadir, "UKB-PPP") |> dir() |> grep("OID", x=_, value = T) |> sub("_.+$", "", x=_) |> unique() |> sort()
  info = get_gene_info(sub("CERT", "CERT1", genes)) |> setDT()
  info[hgnc_symbol == "CERT1", hgnc_symbol := "CERT"]
  d.t = lapply(genes, get_independent_cis_pqtl) |> rbindlist()
  d.t[order(CHROM, GENPOS)] |> fwrite(file.path(outdir, "ukb-ppp.instruments.csv"))

}
