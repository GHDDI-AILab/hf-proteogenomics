## ---- packages ---------------------------------------------------------------
library(data.table)
library(parallel)
library(survival)
library(fdrtool)

## ---- Cox models -------------------------------------------------------------
cox_model = function(formula, mydata) {
	idx = 1  # the index of the given gene
	fit = coxph(formula, data = mydata)
	ci_lower = exp(confint(fit))[idx, 1]
	ci_upper = exp(confint(fit))[idx, 2]

	data.frame(
		"name" = names(coef(fit))[idx],
		"HR" = exp(coef(fit))[idx] |> round(3),
		"CI" = sprintf("%.3f-%.3f", ci_lower, ci_upper),
		"Pval" = coef(summary(fit))[idx, 5]
	)
}

cox_model_1 = function(gene, mydata) {
	cox_model(as.formula(paste("Surv(time, bin)", "~", gene)), mydata)
}

cox_model_2 = function(gene, mydata) {
	cox_model(as.formula(paste("Surv(time, bin)", "~", gene, " + age + sex + assessment_center + age_completed_education + employment + income + tdi")), mydata)
}

cox_model_3 = function(gene, mydata) {
	cox_model(as.formula(paste("Surv(time, bin)", "~", gene, " + age + sex + assessment_center + age_completed_education + employment + income + tdi + smoking + alcohol + pa + sleep + bmi")), mydata)
}

## ---- load data --------------------------------------------------------------
datadir = "../data"
symbols = file.path(datadir, "symbol.list") |> fread(header = T) |> unlist()
mydata = file.path(datadir, "mydata.csv.gz") |> fread(na.strings = "")

# Type conversion: integer to factor
mydata[, sex := as.factor(sex)]
mydata[, assessment_center := as.factor(assessment_center)]
mydata[, employment := as.factor(employment)]
mydata[, income := as.ordered(income)]
mydata[, alcohol := as.factor(alcohol)]
mydata[, smoking := as.factor(smoking)]
mydata[, pa := as.ordered(pa)]

## ---- Perform Cox and calculate FDR ------------------------------------------
run = function(model) {
	cox = mclapply(symbols, model, mydata = mydata, mc.cores = 32) |> rbindlist()
	fdr_df = cox[, fdrtool(Pval, statistic = "pvalue")]
	cox[, LocalFDR := fdr_df[["lfdr"]]]
	cox[, FDR := fdr_df[["qval"]]]
	cox[order(FDR, Pval)]
}

# Perform survival analysis with Model 1, Model 2, and Model 3
cox1 = run(cox_model_1)
cox2 = run(cox_model_2)
cox3 = run(cox_model_3)

## ---- Cox results ------------------------------------------------------------
outdir = "../results"
cox1 |> fwrite(file.path(outdir, "cox_model1_results.csv.gz"))
cox2 |> fwrite(file.path(outdir, "cox_model2_results.csv.gz"))
cox3 |> fwrite(file.path(outdir, "cox_model3_results.csv.gz"))
