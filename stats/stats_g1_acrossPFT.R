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
#leaf$PFT[leaf$PFT %in% c("ESAV","DSAV")] <- "SAV"
#leaf$PFT[leaf$PFT == "CRO"] <- "C3C"
#fluxnet$PFT[fluxnet$PFT == "WSA"] <- "SAV"
#fluxnet$PFT[fluxnet$PFT == "GRA"] <- "C3G"
#fluxnet$PFT[fluxnet$PFT %in% c("CSH","OSH")] <- "SHB"

list_PFTs <- c("ENF", "EBF", "DBF", "TRF", "SAV", "SHB", "C3G", "C4G", "C3C", "C4C")
isotope <- subset(isotope, PFT %in% list_PFTs)
leaf <- subset(leaf, PFT %in% list_PFTs)
fluxnet <- subset(fluxnet, PFT %in% list_PFTs)

# for isotopes, fit model by PFT and site_ID
# need to take logs because g1 are not normally distributed
hist(log(isotope$g1))
with(isotope,boxplot(log(g1)~PFT))
isotope$Method <- as.factor(isotope$PFT)
iso <- lmer(log(g1) ~ Method + (1|site_ID), data=isotope)
summary(iso)
tuk <- glht(iso, linfct = mcp(Method="Tukey"))
siglet_iso <- cld(tuk)$mcletters$Letters
siglet_iso
# result - all are same except ENF - low; SAV - high; TRF - even higher


### FLUXNET
# use site as random factor
fluxnet$Method <- as.factor(fluxnet$PFT)
with(fluxnet, boxplot(g1~PFT))
flux <- lmer(g1~Method+(1|site),data=fluxnet)
summary(flux)
tuk_flux <- glht(flux, linfct = mcp(Method="Tukey"))
siglet_flux <- cld(tuk_flux)$mcletters$Letters
siglet_flux
# result
#C4 crops less than C3 crops
# C3G greater than DBF and ENF but not EBF or TropRF
# SHB greater than all other than SAV

### LEAF GSX
leaf$Method <- as.factor(leaf$PFT)
# Use Location as random factor
gsx <- lmer(g1~Method+(1|Location),data=leaf)
summary(gsx)
tuk_gsx <- glht(gsx, linfct = mcp(Method="Tukey"))
siglet_gsx <- cld(tuk_gsx)$mcletters$Letters
siglet_gsx
with(leaf,boxplot(g1~drop.levels(PFT)))
# result
# C4G less than C3G
# ENF less than DBF

## How to merge and output these?
siglet_gsx_df <- data.frame(PFT=names(siglet_gsx), gsx=as.vector(siglet_gsx))
siglet_iso_df <- data.frame(PFT=names(siglet_iso), iso=as.vector(siglet_iso))
siglet_flux_df <- data.frame(PFT=names(siglet_flux), flux=as.vector(siglet_flux))

smry <- merge(siglet_gsx_df,siglet_iso_df,all=TRUE)
smry <- merge(smry,siglet_flux_df,all=TRUE)
write.csv(smry, "data/processed/g1_signifletters_across_PFT.csv",
          row.names=FALSE)
