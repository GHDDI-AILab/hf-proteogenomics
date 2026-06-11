## ----setup--------------------------------------------------------------------
library(data.table)
library(magrittr)
library(fdrtool)

library(MendelianRandomization)
library(TwoSampleMR)
library(metafor)
library(coloc)


## ----instruments--------------------------------------------------------------
# Get instruments
datadir = "../data"
exposure_data = file.path(datadir, "instruments.csv") |> fread(na.strings = "")


## ----number_of----------------------------------------------------------------
# Check the number of variants in the exposure data
nrow(exposure_data)

# Check the number of Proteins in the exposure data
exposure_data[!is.na(Protein), Protein |> unique() |> length()]


## ----id.exposure--------------------------------------------------------------
exposure_data[, id.exposure := Protein][, exposure := Protein]


## ----outcomes-----------------------------------------------------------------
# Get outcomes from IEU OpenGWAS database
outcome_id = c(
  "ebi-a-GCST009541",  # a GWAS study for heart failure
  "ukb-d-I9_HEARTFAIL",
  "finn-b-I9_HEARTFAIL"
)

## Extract the outcome data from OpenGWAS API
outcome_data = sapply(
  outcome_id,
  function(i) {
    cat(sprintf("\nExtracting outcome data for: %s\n", i))
    d.t = extract_outcome_data(
      snps = exposure_data$SNP,
      outcomes = i,
      proxies = FALSE
    ) |>
      setDT()
    fwrite(d.t, sprintf("%s.outcome_data.csv", i))
    return(d.t)
  },
  simplify = FALSE
)

# Add sample sizes manually
outcome_data[["finn-b-I9_HEARTFAIL"]
  ][, samplesize.outcome := as.integer(samplesize.outcome)
    ][, samplesize.outcome := (13087 + 195091)]

# Check the number of variants in the outcome data
data.table(
  outcome_id = outcome_id,
  N = sapply(outcome_data, nrow),
  N_Var = sapply(outcome_data, function(d.t) d.t[!is.na(SNP), SNP |> unique() |> length()])
)


## ----harmonise----------------------------------------------------------------
# Merge
harmonised = sapply(
  outcome_id,
  function(i) {
    harmonise_data(
      as.data.frame(exposure_data),
      as.data.frame(outcome_data[[i]])
    ) |>
      setDT()
  },
  simplify = FALSE
)

# Check the number of rows, variants, and proteins after harmonization
data.table(
  outcome_id = outcome_id,
  N_Rows = sapply(
    harmonised, nrow
  ),
  N_Var = sapply(
    harmonised, function(d.t) d.t[!is.na(SNP), length(unique(SNP))]
  ),
  N_Prot = sapply(
    harmonised, function(d.t) d.t[!is.na(Protein), length(unique(Protein))]
  ),
  N_Var.mr_keep = sapply(
    harmonised, function(d.t) d.t[!is.na(SNP) & mr_keep == T, length(unique(SNP))]
  ),
  N_Prot.mr_keep = sapply(
    harmonised, function(d.t) d.t[!is.na(Protein) & mr_keep == T, length(unique(Protein))]
  )
)

lapply(outcome_id, function(i)
  harmonised[[i]] |>
    fwrite(sprintf("%s.harmonised_data.csv", i))
)


## ----steiger------------------------------------------------------------------
# Keep the qualified variants in the harmonised data
harmonised = sapply(harmonised, function(d.t) d.t[mr_keep == TRUE, ], simplify = FALSE)

# Steiger test
steiger = sapply(harmonised, directionality_test, simplify = FALSE)

# Get the correct direction
steiger_1 = sapply(
  steiger,
  function(x) subset(x, correct_causal_direction == TRUE),
  simplify = FALSE
)

harmonised_steiger = sapply(
  outcome_id,
  function(i) subset(harmonised[[i]], harmonised[[i]]$exposure %in% steiger_1[[i]]$exposure),
  simplify = FALSE
)

# Check the number of rows, variants, and proteins
data.table(
  outcome_id = outcome_id,
  N_Rows = sapply(harmonised_steiger, nrow),
  N_Variants = sapply(
    harmonised_steiger, function(d.t) d.t[!is.na(SNP), SNP |> unique() |> length()]
  ),
  N_Proteins = sapply(
    harmonised_steiger, function(d.t) d.t[!is.na(Protein), Protein |> unique() |> length()]
  )
)


## ----MR-----------------------------------------------------------------------
# Mendelian Randomization analysis

#' Create an MRInput object for MendelianRandomization.
harmonised_to_MRInput = function(data) {
  MendelianRandomization::mr_input(
    bx = data$beta.exposure,
    bxse = data$se.exposure,
    by = data$beta.outcome,
    byse = data$se.outcome,
    snps = data$SNP
  )
}

#' Collapse a character vector to a length-1 ID.
id_unique = function(x) {
  unique(x) |> paste(collapse = ";")
}

#' Perform MR for one MR method and one exposure-outcome pair.
mr_by_method = function(data, method) {
  if (method == "mr_wald_ratio") {
    if (nrow(data) == 1) TwoSampleMR::mr(data, method_list = method) else data.frame()
  } else if (nrow(data) <= 1) {
    data.frame()
  } else if (method == "mr_ivw_fe") {
    message(sprintf("Analysing '%s' on '%s'", data$id.exposure[1], data$id.outcome[1]))
    res = MendelianRandomization::mr_ivw(
      harmonised_to_MRInput(data),
      model = "fixed"
    )
    data.frame(
      id.exposure = data$id.exposure |> id_unique(),
      id.outcome = data$id.outcome |> id_unique(),
      outcome = data$outcome |> id_unique(),
      exposure = data$exposure |> id_unique(),
      method = "Inverse variance weighted (fixed effects)",
      nsnp = length(data$SNP),
      b = res@Estimate,
      se = res@StdError,
      pval = res@Pvalue
    )
  } else if (method %in% c("mr_ivw_mre", "mr_ivw_re")) {
    message(sprintf("Analysing '%s' on '%s'", data$id.exposure[1], data$id.outcome[1]))
    res = MendelianRandomization::mr_ivw(
      harmonised_to_MRInput(data),
      model = "random"
    )
    data.frame(
      id.exposure = data$id.exposure |> id_unique(),
      id.outcome = data$id.outcome |> id_unique(),
      outcome = data$outcome |> id_unique(),
      exposure = data$exposure |> id_unique(),
      method = "Inverse variance weighted (random effects)",
      nsnp = length(data$SNP),
      b = res@Estimate,
      se = res@StdError,
      pval = res@Pvalue
    )
  } else {
    TwoSampleMR::mr(data, method_list = method)
  }
}

# Perform MR
mr_by_methods = function(dat, method_list) {
  rbindlist(lapply(method_list, mr_by_method, data = dat))
}

mr_calculation = function(d.t, method_list) {
  setDT(d.t)
  #d.t[, mr_by_methods(.SD, method_list), by = list(id.exposure, id.outcome)]
  groups = d.t[, .N, by = list(id.exposure, id.outcome)]
  rbindlist(lapply(
    1:nrow(groups),
    function(i)
      d.t[id.exposure == groups[i, id.exposure] & id.outcome == groups[i, id.outcome]] |>
      mr_by_methods(method_list)
  ))
}

# Perform MR for multiple datasets with multiple methods
mr_methods = c(
    "mr_wald_ratio",
    "mr_ivw_fe",
    "mr_ivw_mre",
    "mr_egger_regression",
    "mr_penalised_weighted_median"
)

res = sapply(
  harmonised_steiger,
  mr_calculation,
  method_list = mr_methods,
  simplify = FALSE
)

# Check the number of rows, and exposures
data.table(
  outcome_id = outcome_id,
  N_Rows = sapply(res, nrow),
  N_Exposures = sapply(res, function(d.t) d.t[!is.na(id.exposure), id.exposure |> unique() |> length()])
)

lapply(outcome_id, function(i) fwrite(res[[i]], sprintf("%s.mr_results.csv", i)))


## ----heterogeneity------------------------------------------------------------
# Heterogeneity statistics
heterogeneity1 = sapply(
  harmonised_steiger,
  function(x) mr_heterogeneity(x) |> setDT(),
  simplify = FALSE
)

# Check the number of rows
data.table(
  outcome_id = outcome_id,
  N_Rows = sapply(heterogeneity1, nrow)
)

# Check the number of records with heterogeneity
ivw = "Inverse variance weighted"
egger = "MR Egger"

ht = heterogeneity1[[outcome_id[1]]]
data.table(
  mr_method = c(ivw, egger),
  N = c(
    ht[method == ivw, .N],
    ht[method == egger, .N]
  ),
  heterogeneity_count = c(
    ht[method == ivw & Q_pval < 0.05, .N],
    ht[method == egger & Q_pval < 0.05, .N]
  ),
  heterogeneity_percentage = c(
    ht[method == ivw & Q_pval < 0.05, .N] / ht[method == ivw, .N] * 100,
    ht[method == egger & Q_pval < 0.05, .N] / ht[method == egger, .N] * 100
  )
)

# Horizontal pleiotropy
pleiotropy1 = sapply(
  harmonised_steiger,
  function(x) mr_pleiotropy_test(x) |> setDT(),
  simplify = FALSE
)

# Check the number of rows and the number of records with horizontal pleiotropy
data.table(
  outcome_id = outcome_id,
  N_Rows = sapply(pleiotropy1, nrow),
  N_Pleiotropies = sapply(pleiotropy1, function(x) x[pval < 0.05, .N])
)


## ----restrict-----------------------------------------------------------------
# Restrict MR results to a single method per exposure x disease combination.
ivw = "Inverse variance weighted"
egger = "MR Egger"

# Identify exposure IDs with significant heterogeneity in IVW method
hetero.ids = sapply(
  heterogeneity1,
  function(x)
    x[method == ivw & Q_pval < 0.05, id.exposure],
  simplify = FALSE
)

# Identify exposure IDs with significant pleiotropy (MR-Egger intercept test)
pleio.ids = sapply(
  pleiotropy1,
  function(x)
    x[pval < 0.05, id.exposure],
  simplify = FALSE
)

# Add logical column indicating if each exposure has significant heterogeneity
v = lapply(outcome_id, function(i)
  res[[i]][, has_heterogeneity := ifelse(id.exposure %chin% hetero.ids[[i]], T, F)]
)

# Add logical column indicating if each exposure has significant pleiotropy
v = lapply(outcome_id, function(i)
  res[[i]][, has_pleiotropy := ifelse(id.exposure %chin% pleio.ids[[i]], T, F)]
)

# Select Wald ratio results: 
# single SNP analyses using Wald ratio method
wald_ = sapply(
  res,
  function(x) x[
    nsnp == 1
    & method == "Wald ratio",
  ],
  simplify = FALSE
)

# Select IVW fixed effects results: 
# multiple SNPs, no pleiotropy, no heterogeneity
fe_ = sapply(
  res,
  function(x) x[
    nsnp > 1
    & has_pleiotropy == FALSE
    & has_heterogeneity == FALSE
    & method == paste(ivw, "(fixed effects)"),
  ],
  simplify = FALSE
)

# Select IVW random effects results: 
# multiple SNPs, no pleiotropy, but heterogeneity present
re_ = sapply(
  res,
  function(x) x[
    nsnp > 1
    & has_pleiotropy == FALSE
    & has_heterogeneity == TRUE
    & method == paste(ivw, "(random effects)"),
  ],
  simplify = FALSE
)

# Select MR Egger results: 
# multiple SNPs with significant pleiotropy
egger_ = sapply(
  res,
  function(x) x[
    nsnp > 1
    & has_pleiotropy == TRUE
    & method == "MR Egger",
  ],
  simplify = FALSE
)

# Combine all selected results into a single table
res_subset = sapply(
  outcome_id,
  function(i) rbindlist(list(wald_[[i]], fe_[[i]], re_[[i]], egger_[[i]])),
  simplify = FALSE
)

# Check the number of rows and exposures
data.table(
  outcome_id = outcome_id,
  N_Rows = sapply(res_subset, nrow),
  N_Exposures = sapply(res_subset, function(x) x[, id.exposure |> unique() |> length()])
)

# Check the numbers of results with pval < 0.05
data.table(
  outcome_id = outcome_id,
  N = sapply(res_subset, function(x) x[pval < 0.05, .N])
)


## ----qval---------------------------------------------------------------------
# Estimate false discovery rates
fdr_df_list = sapply(
  res_subset,
  function(x) x[, fdrtool(pval, statistic = "pvalue")],
  simplify = FALSE
)
v = lapply(
  outcome_id,
  function(i) res_subset[[i]][, FDR := fdr_df_list[[i]][["qval"]]]
)

# View
fdr_df_list[[outcome_id[1]]] |> head(5)

# Check the results with FDR < 0.05
res_subset[[outcome_id[1]]][FDR < 0.05, ][order(FDR, pval)]
res_subset[[outcome_id[2]]][FDR < 0.05, ][order(FDR, pval)]
res_subset[[outcome_id[3]]][FDR < 0.05, ][order(FDR, pval)]


# Output
lapply(outcome_id, function(i)
  res_subset[[i]][order(FDR, pval, -nsnp)] |>
    fwrite(sprintf("%s.mr_results.fdrtool.csv", i))
)


## ----meta---------------------------------------------------------------------
# Subset and tidy the MR result table
tidy = function(d.t, method_patt) {
  d.t[grepl(method_patt, method), list(exposure, method, nsnp, b, se)]
}

# Run meta analysis for one exposure-outcome contrast
run_meta = function(beta, se, meth = "FE") {
  rma(yi = beta, sei = se, method = meth) %$%
    list(b = b, se = se, pval = pval, I2 = I2, H2 = H2)
}

# Run multiple meta analyses and fdr calculation
run_meta_and_fdr = function(datasets, mr_method, meta_method = "FE") {
  res = 
    lapply(datasets, tidy, method_patt = mr_method) %>%
    rbindlist() %>%
    .[, run_meta(b, se, meth = meta_method), by = list(exposure)]
    # .[, FDR := p.adjust(pval, method = "fdr")]
  fdr_df = res[, fdrtool(pval, statistic = "pvalue")]
  res[, FDR := fdr_df[["qval"]]]
}

## ebi and finn
ids.ebi.finn =
  c("ebi-a-GCST009541", "finn-b-I9_HEARTFAIL")
res.wald.ebi.finn =
  run_meta_and_fdr(res_subset[ids.ebi.finn], mr_method = "Wald ratio", meta_method = "FE")
res.ivw.ebi.finn =
  run_meta_and_fdr(res_subset[ids.ebi.finn], mr_method = "Inverse variance weighted", meta_method = "FE")
res.ebi.finn =
  run_meta_and_fdr(res_subset[ids.ebi.finn], mr_method = ".", meta_method = "FE")


## ebi and ukb
ids.ebi.ukb =
  c("ebi-a-GCST009541", "ukb-d-I9_HEARTFAIL")
res.wald.ebi.ukb =
  run_meta_and_fdr(res_subset[ids.ebi.ukb], mr_method = "Wald ratio", meta_method = "FE")
res.ivw.ebi.ukb =
  run_meta_and_fdr(res_subset[ids.ebi.ukb], mr_method = "Inverse variance weighted", meta_method = "FE")
res.ebi.ukb =
  run_meta_and_fdr(res_subset[ids.ebi.ukb], mr_method = ".", meta_method = "FE")


## finn and ukb
ids.finn.ukb =
  c("finn-b-I9_HEARTFAIL", "ukb-d-I9_HEARTFAIL")
res.wald.finn.ukb =
  run_meta_and_fdr(res_subset[ids.finn.ukb], mr_method = "Wald ratio", meta_method = "FE")
res.ivw.finn.ukb =
  run_meta_and_fdr(res_subset[ids.finn.ukb], mr_method = "Inverse variance weighted", meta_method = "FE")
res.finn.ukb =
  run_meta_and_fdr(res_subset[ids.finn.ukb], mr_method = ".", meta_method = "FE")


## ebi, finn, and ukb
ids.ebi.finn.ukb =
  c("ebi-a-GCST009541", "finn-b-I9_HEARTFAIL", "ukb-d-I9_HEARTFAIL")
res.wald.ebi.finn.ukb =
  run_meta_and_fdr(res_subset[ids.ebi.finn.ukb], mr_method = "Wald ratio", meta_method = "FE")
res.ivw.ebi.finn.ukb =
  run_meta_and_fdr(res_subset[ids.ebi.finn.ukb], mr_method = "Inverse variance weighted", meta_method = "FE")
res.ebi.finn.ukb =
  run_meta_and_fdr(res_subset[ids.ebi.finn.ukb], mr_method = ".", meta_method = "FE")


## ----save meta results--------------------------------------------------------
## ebi and finn
paste(ids.ebi.finn, collapse = ".") |>
  paste("meta_analysis.wald.csv", sep = ".") |>
  fwrite(res.wald.ebi.finn[order(FDR, pval)], file = _)

paste(ids.ebi.finn, collapse = ".") |>
  paste("meta_analysis.ivw.csv", sep = ".") |>
  fwrite(res.ivw.ebi.finn[order(FDR, pval)], file = _)

paste(ids.ebi.finn, collapse = ".") |>
  paste("meta_analysis.csv", sep = ".") |>
  fwrite(res.ebi.finn[order(FDR, pval)], file = _)

## ebi and ukb
paste(ids.ebi.ukb, collapse = ".") |>
  paste("meta_analysis.wald.csv", sep = ".") |>
  fwrite(res.wald.ebi.ukb[order(FDR, pval)], file = _)

paste(ids.ebi.ukb, collapse = ".") |>
  paste("meta_analysis.ivw.csv", sep = ".") |>
  fwrite(res.ivw.ebi.ukb[order(FDR, pval)], file = _)

paste(ids.ebi.ukb, collapse = ".") |>
  paste("meta_analysis.csv", sep = ".") |>
  fwrite(res.ebi.ukb[order(FDR, pval)], file = _)

## finn and ukb
paste(ids.finn.ukb, collapse = ".") |>
  paste("meta_analysis.wald.csv", sep = ".") |>
  fwrite(res.wald.finn.ukb[order(FDR, pval)], file = _)

paste(ids.finn.ukb, collapse = ".") |>
  paste("meta_analysis.ivw.csv", sep = ".") |>
  fwrite(res.ivw.finn.ukb[order(FDR, pval)], file = _)

paste(ids.finn.ukb, collapse = ".") |>
  paste("meta_analysis.csv", sep = ".") |>
  fwrite(res.finn.ukb[order(FDR, pval)], file = _)

## ebi, finn, and ukb
paste(ids.ebi.finn.ukb, collapse = ".") |>
  paste("meta_analysis.wald.csv", sep = ".") |>
  fwrite(res.wald.ebi.finn.ukb[order(FDR, pval)], file = _)

paste(ids.ebi.finn.ukb, collapse = ".") |>
  paste("meta_analysis.ivw.csv", sep = ".") |>
  fwrite(res.ivw.ebi.finn.ukb[order(FDR, pval)], file = _)

paste(ids.ebi.finn.ukb, collapse = ".") |>
  paste("meta_analysis.csv", sep = ".") |>
  fwrite(res.ebi.finn.ukb[order(FDR, pval)], file = _)


# Collect the significant results
mr.prots = c(
  res_subset[[outcome_id[1]]][FDR < 0.05, exposure],
  res_subset[[outcome_id[2]]][FDR < 0.05, exposure],
  res_subset[[outcome_id[3]]][FDR < 0.05, exposure]
) |>
  sort() |> unique()

meta.prots = c(
  res.wald.ebi.finn[FDR < 0.05, exposure],
  res.ivw.ebi.finn[FDR < 0.05, exposure],
  res.ebi.finn[FDR < 0.05, exposure],

  res.wald.ebi.ukb[FDR < 0.05, exposure],
  res.ivw.ebi.ukb[FDR < 0.05, exposure],
  res.ebi.ukb[FDR < 0.05, exposure],

  res.wald.finn.ukb[FDR < 0.05, exposure],
  res.ivw.finn.ukb[FDR < 0.05, exposure],
  res.finn.ukb[FDR < 0.05, exposure],

  res.wald.ebi.finn.ukb[FDR < 0.05, exposure],
  res.ivw.ebi.finn.ukb[FDR < 0.05, exposure],
  res.ebi.finn.ukb[FDR < 0.05, exposure]
) |>
  sort() |> unique()

# View
data.frame(
  Prot = c(
    paste(mr.prots, collapse = "; "),
    paste(meta.prots, collapse = "; ")
  ),
  Item = c("MR", "Meta")
)


## ----colocalization-----------------------------------------------------------
# Colocalization analysis
get_coloc_summary = function(dat, cc_ratio, p12 = 1e-5) {
  tryCatch({

    # 确保 N 是标量
    N1 = unique(dat$samplesize.exposure)
    N2 = unique(dat$samplesize.outcome)

    # 计算 MAF
    maf = pmin(dat$eaf.exposure, 1 - dat$eaf.exposure)

    res = coloc.abf(
      dataset1 = list(
        beta = dat$beta.exposure,
        varbeta = dat$se.exposure^2,
        N = N1,
        type = "quant",
        MAF = maf,
        snp = dat$SNP
      ),
      dataset2 = list(
        beta = dat$beta.outcome,
        varbeta = dat$se.outcome^2,
        N = N2,
        type = "cc",
        s = cc_ratio,
        snp = dat$SNP
      ),
      p12 = p12
    )

    out = data.table(t(as.data.frame(res$summary)))
    out[, exposure := unique(dat$exposure)]

    return(out)

  }, error = function(e) {
    message("coloc error: ", e$message)
    return(NULL)
  })
}

perform_coloc = function(exposures, cc.ratio, data) {
  rbindlist(lapply(
    exposures,
    function(gene) data[exposure == gene & mr_keep == TRUE, ] |> get_coloc_summary(cc.ratio)
  ))
}

# Combine the results from MR and from meta-analysis
prots = c(mr.prots, meta.prots) |> sort() |> unique()

# ebi
coloc.ebi = perform_coloc(
  prots, cc.ratio = 47309 / 977323, data = harmonised[[outcome_id[1]]]
)[order(-PP.H4.abf)]

coloc.ebi %T>% fwrite("colocalization.ebi-a-GCST009541.csv")
coloc.ebi[PP.H4.abf > 0.5, ]

# ukb
coloc.ukb = perform_coloc(
  prots, cc.ratio = 1405 / (1405 + 359789), data = harmonised[[outcome_id[2]]]
)[order(-PP.H4.abf)]

coloc.ukb %T>% fwrite("colocalization.ukb-d-I9_HEARTFAIL.csv")
coloc.ukb[PP.H4.abf > 0.5, ]

# finn
coloc.finn = perform_coloc(
  prots, cc.ratio = 13087 / (13087 + 195091), data = harmonised[[outcome_id[3]]]
)[order(-PP.H4.abf)]

coloc.finn %T>% fwrite("colocalization.finn-b-I9_HEARTFAIL.csv")
coloc.finn[PP.H4.abf > 0.5, ]


## ----Integrate Cox, MR, and Coloc---------------------------------------------
# Get the candidates with colocalization evidences
candidates = rbindlist(list(
  coloc.ebi [PP.H4.abf > 0.5, list(exposure, PP.H4.abf)],
  coloc.ukb [PP.H4.abf > 0.5, list(exposure, PP.H4.abf)],
  coloc.finn[PP.H4.abf > 0.5, list(exposure, PP.H4.abf)]
))[order(-PP.H4.abf)][, head(.SD, 1), by = list(exposure)]

candidates


# Cox analysis results
cox_model_3 = fread("../../1.Cohort-study/results/cox_model3_results.csv.gz")
cox_model_3[name %chin% candidates$exposure, ]
cox_model_3[name %chin% candidates$exposure & FDR < 0.05, ]

# Mendelian randomization results
res_subset[["ebi-a-GCST009541"]][candidates, on = "exposure"
  ][order(FDR), list(exposure, PP.H4.abf, nsnp, b, se, pval, FDR)][FDR < 0.05]

res_subset[["ukb-d-I9_HEARTFAIL"]][candidates, on = "exposure"
  ][order(FDR), list(exposure, PP.H4.abf, nsnp, b, se, pval, FDR)][FDR < 0.05]

res_subset[["finn-b-I9_HEARTFAIL"]][candidates, on = "exposure"
  ][order(FDR), list(exposure, PP.H4.abf, nsnp, b, se, pval, FDR)][FDR < 0.05]

