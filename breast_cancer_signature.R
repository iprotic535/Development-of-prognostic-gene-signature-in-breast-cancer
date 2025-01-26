library(survival)
library(dplyr)

# Load gene expression data from different cohorts
FRA_expr = read.delim("FRA_expr.txt", row.names = "gene")
GER_expr = read.delim("GER_expr.txt", row.names = "gene")
SWE_expr = read.delim("SWE_expr.txt", row.names = "gene")
TWN_expr = read.delim("TWN_expr.txt", row.names = "gene")
USA1_expr = read.delim("USA1_expr.txt", row.names = "gene")
USA2_expr = read.delim("USA2_expr.txt", row.names = "gene")

# Load phenotype data corresponding to each cohort
FRA_pheno = read.delim("FRA_pheno.txt")
GER_pheno = read.delim("GER_pheno.txt")
SWE_pheno = read.delim("SWE_pheno.txt")
TWN_pheno = read.delim("TWN_pheno.txt")
USA1_pheno = read.delim("USA1_pheno.txt")
USA2_pheno = read.delim("USA2_pheno.txt")

# Identify prognostic genes in FRA cohort using Cox regression analysis
z_FRA = match(colnames(FRA_expr), FRA_pheno$sample)
pheno = FRA_pheno[z_FRA,]
hr = rep(NA, nrow(FRA_expr))
p = rep(NA, nrow(FRA_expr))
for (i in 1:nrow(FRA_expr)) {
  g = summary(coxph(Surv(pheno$survival, pheno$death) ~ t(FRA_expr[i,])))
  hr[i] = g$coefficients[2]
  p[i] = g$coefficients[5]
}
# Apply Bonferroni correction for multiple comparisons
result = data.frame(rownames(FRA_expr), hr, p)
sig = result[p.adjust(result$p, method="bonferroni")<0.05,]
signature = data.frame(sig[,1], sign(log2(sig$hr)))
colnames(signature) = c("gene", "weight")
write.table(signature, "signatureFRA.txt", row.names=F, quote=F, sep="	")

# Identify prognostic genes in TWN cohort
z_TWN = match(colnames(TWN_expr), TWN_pheno$sample)
pheno = pheno[z_TWN,]
hr = rep(NA, nrow(TWN_expr))
p = rep(NA, nrow(TWN_expr))
for (i in 1:nrow(TWN_expr)) {
  g = summary(coxph(Surv(TWN_pheno$survival, TWN_pheno$death) ~ t(TWN_expr[i,])))
  hr[i] = g$coefficients[2]
  p[i] = g$coefficients[5]
}
# Apply Benjamini-Hochberg correction for multiple testing
result = data.frame(rownames(TWN_expr), hr, p)
sig = result[p.adjust(result$p, method="BH")<0.05,]
signature = data.frame(sig[,1], sign(log2(sig$hr)))
colnames(signature) = c("gene", "weight")
write.table(signature, "signatureTWN.txt", row.names=F, quote=F, sep="	")

# Merge significant genes from USA2 and GER discovery cohorts
sig_USA2 = read.delim("signatureUSA2.txt")
sig_GER = read.delim("signatureGER.txt")
sig_USA2_GER = right_join(sig_USA2, sig_GER, by = "gene")
sig = sig_USA2_GER %>% select("gene", "weight.x")
write.table(sig, "signature_USA2_GER.txt", row.names=F, quote=F, sep="	")

# Validate the gene signature in the FRA cohort
z = match(colnames(FRA_expr), FRA_pheno$sample)
FRA_pheno = FRA_pheno[z,]
signature = read.delim("signature_USA2_GER.txt")
score = risk.score(as.matrix(FRA_expr), signature)
cl = ifelse(score <= 0, 0, 1)
# Cox regression for validation
g = summary(coxph(Surv(FRA_pheno$survival, FRA_pheno$death) ~ cl))
print(g)
# Kaplan-Meier survival analysis
survdiff(Surv(FRA_pheno$survival, FRA_pheno$death) ~ cl)
plot(survfit(Surv(FRA_pheno$survival, FRA_pheno$death) ~ cl), col=c("blue","red"), 
     main="Validation Cohort France", font.main=1, xlab="Survival time (years)", 
     ylab="Survival proportion", axes=TRUE)
legend("topright", lty=1, col=c("red", "blue"), c("High risk score", "Low risk score"), bty="n")
