# Assign covid hgi top variants to subtypes depending whether they affect 
# (1) susceptibility to infection or
# (2) severity of the disease.
# Matti Pirinen (matti.pirinen@helsinki.fi)
# 3-Dec-2021

###
### Part 1. Pre-processing the data
###

#Input data are GWAS results from two GWAS:
#    B2 = hospitalized covid cases vs. controls
#    C2 = all reported infected vs. controls

# Next: Compute expected correlation in GWAS effect estimates between B2 and C2 GWAS.
#        B2 and C2 share majority of controls and we assume that B2 cases are among C2 cases
#        Because of these shared samples, the effect estimators are correlated.

# Read in the study information. Sample sizes of cases and controls are given for
# inf = C2 GWAS (infected)
# sev = B2 GWAS (severely ill)
# nonsev = imaginary set of "non-severe infected", i.e., all infected that are not included in B2
# NOTE: nonsev is a crude approximation since we DON'T know whether those infected
#       who were not in B2 have been also hospitalized. 

stu = read.table("covid_hgi_v6_studies_in_subtype_analyses.txt", as.is = T, header = T)
stu[2:ncol(stu)] = lapply(stu[2:ncol(stu)], as.numeric)
head(stu, 3)
colSums(stu[,2:7]) # matches with sample sizes in C2 and B2 GWAS

#Compute effective sample sizes
neff.inf = stu$inf_ncase * stu$inf_ncontrol / (stu$inf_ncase + stu$inf_ncontrol)
neff.sev = stu$sev_ncase * stu$sev_ncontrol / (stu$sev_ncase + stu$sev_ncontrol)
neff.nonsev = stu$nonsev_ncase * stu$nonsev_ncontrol / (stu$nonsev_ncase + stu$nonsev_ncontrol)
neff.inf[is.nan(neff.inf)] = 0
neff.sev[is.nan(neff.sev)] = 0
neff.nonsev[is.nan(neff.nonsev)] = 0

# Correlation between C2 and B2 GWAS effect estimators due to overlapping samples.
# Formula based on 
# Bhattacharjee et al. (2012) Am J Hum Genet. 2012 May 4;90(5):821-35. 
# doi: 10.1016/j.ajhg.2012.03.015. PMID: 22560090; PMCID: PMC3376551.
r = sqrt( neff.inf * neff.sev ) * 
  (apply(stu[,c("inf_ncase","sev_ncase")], 1, min) / (stu[,"inf_ncase"]*stu[,"sev_ncase"]) +
   apply(stu[,c("inf_ncontrol","sev_ncontrol")], 1, min) / (stu[,"inf_ncontrol"]*stu[,"sev_ncontrol"]))
r[is.nan(r)] = 0

# Weighted correlation estimate across studies
sum(sqrt(neff.inf * neff.sev) * r) / sqrt(sum(neff.inf) * sum(neff.sev))
# = 0.4539485

# Which proportion of neff in C2 would come from B2 individuals?
sum(neff.sev)/sum(neff.inf)
# = 0.2083978

#This proportion 'w' can be used to define the expected effect size observed in C2
# analysis for a variant that affected severity of disease (but not susceptibility to infection)
# Namely, beta_C2 = w*beta_B2


#We can compute the value for `w` also by taking into account the cohort wise correlations.
# For that,
# let's estimate beta and SE for an imaginary GWAS of non-severe cases vs. controls. 
# Such GWAS would result if we removed B2 cases from C2 cases.
# Hence we call it also C2_B2 GWAS (and also as "nonsevere" GWAS).
# Note that we haven't ever done such a GWAS but here we simply approximate results
# by assuming that C2 was a meta-analysis between B2 and this imaginary C2_B2 GWAS.

# Correlation between effect estimators of sev and nonsev GWAS across studies
#  (These two do not share any cases but do share controls.)
r = sqrt( neff.nonsev * neff.sev ) * 
  (apply(stu[,c("sev_ncontrol","nonsev_ncontrol")], 1, min) / 
     (stu[,"sev_ncontrol"]*stu[,"nonsev_ncontrol"])) # assumes max overlap in controls
r[is.nan(r)] = 0

# Weighted correlation across studies
r = sum(sqrt(neff.sev*neff.nonsev) * r) / sqrt(sum(neff.sev)*sum(neff.nonsev))
r
#r = 0.02157236

#Let's read in GWAS meta-analysis results from such B2 and C2 meta-analyses
# that use only cohorts that have contributed data to both B2 and C2
dat = read.table("covid_hgi_v6_B2_C2_common.tsv", as.is = T, header = T)

#Estimate SE for nonsev analysis using observed SEs and neff estimates: 
se.nonsev.1 = dat$C2_sebeta * sqrt(sum(neff.inf) / sum(neff.nonsev)) #one estimate
se.nonsev.2 = dat$B2_sebeta * sqrt(sum(neff.sev) / sum(neff.nonsev)) #another estimate
se.nonsev = apply(cbind(se.nonsev.1, se.nonsev.2), 1 , max) #robust estimate for SE

# Assume that C2 scan was combination of B2 (SEV) and C2_B2 (NSEV).
# The 'w' below is the weight given to B2 when it is combined with C2_B2 to make C2:
w = (se.nonsev^2 - dat$B2_sebeta * se.nonsev * r)/(dat$B2_sebeta^2 + se.nonsev^2 - 2*dat$B2_sebeta*se.nonsev*r)
# Weigths of these studies would be approximately (assuming same MAF and beta=0):
c(SEV = median(w), NSEV = 1-median(w))
#SEV       NSEV 
#0.1994855 0.8005145

# We can also compute beta and P-value for imaginary C2_B2 scan:
b.nonsev = (dat$C2_beta - w*dat$B2_beta) / (1 - w)
p.nonsev = pchisq( (b.nonsev / se.nonsev)^2, df = 1, lower = F)

#Write to file the results
dat = cbind(dat[c(1:5,14,6:8,15:17)], b.nonsev, se.nonsev, p.nonsev)
names(dat)[c(6,13:15)] = c("RSID", "C2_B2_beta", "C2_B2_sebeta", "C2_B2_pval")

write.table(dat, file = "covid_hgi_v6_top_shared_c2_b2.txt", quote = FALSE, row.names = FALSE)




###
###  PART 2. SUBTYPE ANALYSES 
###

# Includes:
#  Bayesian model comparison between infection susceptibility (INF) and 
#  disease severity model (SEV).
#  Each model has prior prob of 0.50 for each variant.
#  INF model is correlated (0.98) effect model.
#  Severity model says that C2 = 0.20*B2 effect
#
# Explanation of the Bayesian approach: 
#  https://www.mv.helsinki.fi/home/mjxpirin/GWAS_course/material/GWAS9.html

# Read in functions:
log.dmvnorm <- function(x, mu = rep(0, length(x)), S = diag(1, length(x)) ){
  #returns log of density of MV-Normal(mean = mu, var = S) at x 
  K = length(mu)
  stopifnot(all(dim(S) == K))
  stopifnot(length(x) == K)
  chol.S = chol(S) #Cholesky decomposition
  log.det = 2*sum(log(diag(chol.S))) #log of det(S)
  inv.chol.S = solve(t(chol.S)) #inverse of cholesky^T
  return(-K/2*log(2*pi) - 0.5*(log.det + crossprod(inv.chol.S %*% as.numeric(x-mu))))
}

abf.mv <- function(b.est, Sigmas, prior = rep(1,length(Sigmas))){
  #Returns posterior probabilities of the models listed in Sigmas by their 
  #   total variance matrix (= sum of prior + likelihood variance matrices)
  #Returns also approximate Bayes factors (ABFs) w.r.t the first model in Sigmas.
  
  M = length(Sigmas) #number of models
  K = length(b.est) #number of studies
  prior = prior/sum(prior)
  log.abf = sapply(Sigmas, function(x){log.dmvnorm(b.est, S = x)})
  abf = exp(log.abf - log.abf[1]) #abf w.r.t the first model
  posterior = prior*abf
  posterior = posterior/sum(posterior)
  return(list(posterior = posterior, abf = abf))
}

#Input data has effect estimates for COVID risk variants.
# Estimates are given for 3 phenotypes: 
# C2 (infection vs population, "susceptibility scan"), 
# B2 (hospitalized vs population "severity scan")
# C2_B2 (statistically inferred results for C2 \ B2)
x = read.table("covid_hgi_v6_top_shared_c2_b2.txt", as.is = T, header = T)
# We only use estimates for B2 and C2 directly (i.e. we don't use atatistically inferred C2-B2)

#column names
se.b2 = "B2_sebeta"
se.c2 = "C2_sebeta"
beta.b2 = "B2_beta"
beta.c2 = "C2_beta"


r = 0.4539485 # correlation for the likelihood; was computed from overlaps in case-control counts
w = 0.20 #weight of B2 in beta observed in C2 (computed from effective sample size ratio)
r.inf = 0.999 # How correlated the real effects in two studies are?
r.sev = 0.987 # These are visualized below.
tau2 = 0.1^2 #prior variance of effect size -- same for all non-null models
S.inf = tau2 * matrix(c(1, r.inf, r.inf, 1), 2, 2) #correlated effects but not the same
S.sev = tau2 * matrix(c(1, r.sev*w, r.sev*w, w^2), 2, 2) #Effect only in B2 but not in C2-B2
priors = c(0.5, 0.5) # INF, SEV

pr = c() #model posterior probabilties
for(ii in 1:nrow(x)){
  #Sigma is the variance of likelihood function 
  Sigma = diag(x[ii,c(se.b2,se.c2)]) %*% matrix(c(1,r,r,1),2,2) %*% diag(x[ii,c(se.b2,se.c2)])
  Var.matrices = list(Sigma + S.inf, Sigma + S.sev)
  pr = rbind(pr, abf.mv(b.est = x[ii,c(beta.b2, beta.c2)], 
                        Sigmas = Var.matrices, prior = priors)$posterior)
}
colnames(pr) = c("inf","sev")

y = cbind(x, pr)
names(y)[16:17] = c("prob_inf","prob_sev")
y = y[y[,"RSID"]!="rs4767023",] #remove one SNP that is in LD with another -- 23 remain
nrow(y)
write.table(y, file ="covid_hgi_v6_topvariants_prob.txt", 
            quote = FALSE, row.names = FALSE)


#Show the prior distribution of effect size under INF and SEV models
# using simulations from the priors

png("covid_b2_c2_effect_priors.png", width = 600, height = 600)
library(mvtnorm)
npoints = 10000
x.inf = rmvnorm(npoints, sigma = S.inf)
x.sev = rmvnorm(npoints, sigma = S.sev)
cols = rep(c("orange","red"), rep(npoints,2))
plot(rbind(x.inf, x.sev), col = cols, pch = 3, cex = 0.5,
     xlab = "B2", ylab ="C2", 
     sub = paste0("r.inf=",r.inf," r.sev=",r.sev))
abline(0,1, col = "orange", lty = 2)
abline(0,w, col = "red", lty = 2)
dev.off()

summary(abs(x.inf[,1]-x.inf[,2]))
quantile(abs(x.inf[,1]-x.inf[,2]),probs = 0.95)
summary(abs(w*x.sev[,1]-x.sev[,2]))
quantile(abs(w*x.sev[,1]-x.sev[,2]),probs = 0.95)



##
## Make plots of effect sizes
## 

y = read.table("covid_hgi_v6_topvariants_prob.txt", as.is = T, header = T)
n = nrow(y)

#flip effects so that B2 has a positive beta
ii = (y[,"B2_beta"] < 0)
y[ii,grep("_beta",names(y))] = -y[ii,grep("_beta",names(y))]
#Flip also alleles
tmp = y[ii,"REF"]
y[ii,"REF"] = y[ii,"ALT"]
y[ii,"ALT"] = tmp

library(ggplot2)
library(ggrepel)
library(gridExtra)
set.seed(42)

pr.thr = 0.99
ii.sev = (y[,"prob_sev"] > pr.thr)
ii.inf = (y[,"prob_inf"] > pr.thr)

group = rep("other",n)
group[ii.sev] = "SEV"
group[ii.inf] = "INF"

ii = (ii.sev | ii.inf)
texts = rep("", n)
texts[ii] = y$RSID[ii]

ran = range(-0.1, 0.65, y[,"B2_beta"], y[,"C2_beta"])

p1 <- ggplot(y, aes(B2_beta, C2_beta, label = texts, color = group)) +
  geom_point() + geom_text_repel(max.overlaps = 50, show.legend = FALSE) + xlim(ran) + ylim(ran) + 
  geom_abline(linetype = "dashed") + geom_abline(slope = 0.20, linetype = "dashed") +
  theme(text = element_text(size = 10)) 

#95% confidence intervals
y.ci.min = y[,"C2_beta"]-1.96*y[,"C2_sebeta"]
y.ci.max = y[,"C2_beta"]+1.96*y[,"C2_sebeta"]
x.ci.min = y[,"B2_beta"]-1.96*y[,"B2_sebeta"]
x.ci.max = y[,"B2_beta"]+1.96*y[,"B2_sebeta"]

p2 <- ggplot(y, aes(B2_beta, C2_beta, label = texts, color = group)) +
  geom_point() + geom_errorbar(aes(ymin = y.ci.min, ymax = y.ci.max), width=0) + 
  geom_errorbar(aes(xmin = x.ci.min, xmax = x.ci.max), width = 0) +
  xlim(ran) + ylim(ran) + theme(text = element_text(size = 10)) +
  geom_abline(linetype = "dashed") + geom_abline(slope = 0.20, linetype = "dashed") 

pdf("covid_hgi_v6_b2_c2_effects_plot.pdf", width = 10, height = 4)
grid.arrange(p1, p2, nrow = 1)
dev.off()



