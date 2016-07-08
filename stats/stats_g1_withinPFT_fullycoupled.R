#!/usr/bin/Rscript
# Stats on g1
# - comparison between method (leaf,isotope,fluxnet) on g1, for every PFT
#   separately.
#
# Author: Remko A. Duursma
# Date: 1/03/2016
#

# packages
library(doBy)
library(multcomp)

# for a PFT, get sign. letters.
get_signif_pft <- function(x) {
  dat <- alldf[alldf$PFT == x,]
  dat$Method <- as.factor(dat$Method)  # needed for glht

  # If only one method return emtpy strings
  if(nlevels(dat$Method) == 1) {
    return(data.frame(Method=c("fluxnet","isotope","leaf"), PFT=x, siglet=""))
  }
  # anova
  fit <- lm(log(g1) ~ Method, data=dat)

  # Tukey multiple comparison
  tuk <- glht(fit, linfct = mcp(Method="Tukey"))

  # significance letters
  siglet <- cld(tuk)$mcletters$Letters
  siglet <- data.frame(Method=names(siglet), siglet=unname(siglet),
                       stringsAsFactors = FALSE)

  # Make sure we always have length 3 letters (in case of missing methods)
  d <- data.frame(PFT=x, Method=c("fluxnet","isotope","leaf"),
                  stringsAsFactors=FALSE)
  d <- merge(d, siglet, all=TRUE)
  d[is.na(d)] <- ""

  return(d)
}

isotope <- read.csv("data/processed/g1_isotope_screened.csv",
                    stringsAsFactors=FALSE)
leaf <- read.csv("data/processed/g1_leaf_gas_exchange.csv",
                 stringsAsFactors=FALSE)
#fluxnet <- read.csv("data/processed/g1_fluxnet_screened_matching_sites.csv",
#                    stringsAsFactors=FALSE)
fluxnet <- read.csv("data/processed/g1_fluxnet_screened.csv",
                    stringsAsFactors=FALSE)


# Fix levels of PFT.
#leaf$PFT[leaf$PFT %in% c("ESAV","DSAV")] <- "SAV"
#fluxnet$PFT[fluxnet$PFT == "WSA"] <- "SAV"
#fluxnet$PFT[fluxnet$PFT == "GRA"] <- "C3G"
#fluxnet$PFT[fluxnet$PFT %in% c("CSH","OSH")] <- "SHB"

# - for isotope, average by species? I have not done that yet, but probably
#   makes little to no difference as the vast majority (78%) of species have
#   one observation.

# - For fluxnet, average g1 by site (because years are not independent).
fluxnet2 <- summaryBy(. ~site, FUN=mean, data=fluxnet, keep.names=TRUE, id=~PFT)

# - Make single dataset with 'method' (isotope, leaf, fluxnet), and g1.
df1 <- isotope[,c("PFT","g1")]
df1$Method <- "isotope"
df2 <- leaf[,c("PFT","g1")]
df2$Method <- "leaf"
df3 <- fluxnet2[,c("PFT","g1")]
df3$Method <- "fluxnet"

alldf <- rbind(df1, df2, df3)
alldf$PFT <- as.factor(alldf$PFT)

# Inspect sample size by method and PFT
# with(alldf, table(PFT, Method))

# Get all signif letters by PFT
allsig <- lapply(levels(alldf$PFT),get_signif_pft)
allsig <- do.call(rbind, allsig)

# export
write.csv(allsig, "data/processed/g1_siglet_withinPFT_coupled.csv",
          row.names=FALSE)
