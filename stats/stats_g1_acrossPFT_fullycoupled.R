#!/usr/bin/Rscript
# Stats on g1
# - comparison between PFTs for the three methods separately.
#
# Author: Belinda E. Medlyn
# Date: April 2016
#

path = "/Users/bmedlyn/research/g1_leaf_ecosystem_scale"
user <- Sys.info()[["user"]]
path <- gsub("bmedlyn", user, path)

setwd(path)

# packages
library(doBy)
library(multcomp)
library(lme4)
library(LMERConvenienceFunctions)

# read data
isotope <- read.csv("data/processed/g1_isotope_screened.csv",
                    stringsAsFactors=FALSE)
leaf <- read.csv("data/processed/g1_leaf_gas_exchange.csv",
                    stringsAsFactors=FALSE)
fluxnet <- read.csv("data/processed/g1_fluxnet_screened.csv",
                    stringsAsFactors=FALSE)

# Fix levels of PFT.
leaf$PFT[leaf$PFT %in% c("ESAV","DSAV")] <- "SAV"
leaf$PFT[leaf$PFT == "CRO"] <- "C3C"
leaf$PFT[leaf$PFT == "TropRF"] <- "TRF"
fluxnet$PFT[fluxnet$PFT == "WSA"] <- "SAV"
fluxnet$PFT[fluxnet$PFT == "GRA"] <- "C3G"
fluxnet$PFT[fluxnet$PFT %in% c("CSH","OSH")] <- "SHB"
isotope$PFT[isotope$PFT == "TropRF"] <- "TRF"
fluxnet$PFT[fluxnet$PFT == "TropRF"] <- "TRF"

list_PFTs <- c("ENF", "EBF", "DBF", "TRF", "SAV", "SHB", "C3G", "C4G", "C3C", "C4C")
isotope <- subset(isotope, PFT %in% list_PFTs)
leaf <- subset(leaf, PFT %in% list_PFTs)
fluxnet <- subset(fluxnet, PFT %in% list_PFTs)

# for isotopes, fit model by PFT and site_ID
# need to take logs because g1 are clearly not normally distributed
hist(isotope$g1)
hist(log(isotope$g1))
qqnorm(isotope$g1)
qqnorm(log(isotope$g1))

with(isotope,boxplot(log(g1)~PFT,main="Log Isotope"))
iso_smry <- summaryBy(g1~PFT+site_ID,data=isotope, FUN=mean)
with(iso_smry,boxplot(log(g1.mean)~PFT,main="Log Isotope Site Means"))
iso_n <- with(iso_smry,table(PFT))
isotope$Method <- as.factor(isotope$PFT)
iso_log <- lmer(log(g1) ~ Method + (1|site_ID), data=isotope)
logLik(iso_log)
tuk_log <- glht(iso_log, linfct = mcp(Method="Tukey"))
siglet_iso_log <- cld(tuk_log)$mcletters$Letters
siglet_iso_log
# result - all are same except ENF - low; SAV - high; TRF - even higher


### FLUXNET
# use site as random factor
hist(fluxnet$g1)
hist(log(fluxnet$g1))
qqnorm(fluxnet$g1)
qqnorm(log(fluxnet$g1))

fluxnet$Method <- as.factor(fluxnet$PFT)
with(fluxnet, boxplot(g1~PFT,main="Fluxnet",ylab="g1"))
flux_smry <- summaryBy(g1~PFT+site,data=fluxnet, FUN=mean)
with(flux_smry,boxplot(g1.mean~PFT,main="Fluxnet Site Means",ylab="g1 mean"))
with(flux_smry,boxplot(log(g1.mean)~PFT,main="Log Fluxnet Site Means",ylab="log g1 mean"))
flux_n <- with(flux_smry,table(PFT))
flux <- lmer(g1~Method+(1|site),data=fluxnet)
logLik(flux)
tuk_flux <- glht(flux, linfct = mcp(Method="Tukey"))
siglet_flux <- cld(tuk_flux)$mcletters$Letters
siglet_flux
# result
#C4 crops less than C3 crops
# C3G greater than DBF and ENF but not EBF or TropRF
# SHB greater than all other than SAV
flux_log <- lmer(log(g1)~Method+(1|site),data=fluxnet)
logLik(flux_log)
tuk_flux_log <- glht(flux_log, linfct = mcp(Method="Tukey"))
siglet_flux_log <- cld(tuk_flux_log)$mcletters$Letters
siglet_flux_log

### LEAF GSX
hist(leaf$g1)
hist(log(leaf$g1))
qqnorm(leaf$g1)
qqnorm(log(leaf$g1))

leaf$Method <- as.factor(leaf$PFT)
# Use Location as random factor
with(leaf, boxplot(g1~PFT,main="Leaf GSX",ylim=c(0,15)))
leaf_smry <- summaryBy(g1~PFT+Location,data=leaf, FUN=mean)
with(leaf_smry,boxplot(g1.mean~PFT,main="Leaf GSX Site Means",
                       ylab="g1 mean",ylim=c(0,15)))
with(leaf_smry,boxplot(log(g1.mean)~PFT,main="Log Leaf GSX Site Means",
                       ylab="log g1 mean"))
leaf_n <- with(leaf_smry,table(PFT))
gsx <- lmer(g1~Method+(1|Location),data=leaf)
gsx_log <- lmer(log(g1)~Method+(1|Location),data=leaf)
summary(gsx)
tuk_gsx <- glht(gsx, linfct = mcp(Method="Tukey"))
siglet_gsx <- cld(tuk_gsx)$mcletters$Letters
siglet_gsx
# result
# C4G less than C3G
# ENF less than DBF
tuk_gsx_log <- glht(gsx_log, linfct = mcp(Method="Tukey"))
siglet_gsx_log <- cld(tuk_gsx_log)$mcletters$Letters
siglet_gsx_log
logLik(gsx)
logLik(gsx_log)

## How to merge and output these?
siglet_gsx_log
siglet_iso_log
siglet_flux_log

leaf_n
iso_n
flux_n
