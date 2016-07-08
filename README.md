# Discrepancies between leaf and ecosystem measures of water-use efficiency

[Belinda E. Medlyn](https://bmedlyn.wordpress.com/),
[Martin G. De Kauwe](https://mdekauwe.github.io/),
[Yan-Shih Lin](https://sites.google.com/site/yanshihlin/),
Juergen Knauer,
[Remko A. Duursma](http://www.remkoduursma.com/),
Christopher A. Williams,
Almut Arneth,
Rob Clement, Peter Isaac,
Jean-Marc Limousin,
Maj-Lena Linderson,
Patrick Meir,
Nicolas Martin-StPaul,
Eva van Gorsel,
Lisa Wingate.

## Overview
Repository containing all the code associated with our paper.

## Datasets
* Leaf gas exchange data is from [Lin et al. (2015) Optimal stomatal behaviour around the world. Nature Climate Change, 5, 459–464.](http://www.nature.com/nclimate/journal/v5/n5/full/nclimate2550.html)
* Stable isotope dataset: [Dryad Digital Repository](http://dx.doi.org/10.5061/dryad.3jh61)
* Eddy covariance dataset: [FLUXNET](http://www.fluxdata.org/DataInfo/default.aspx)
* Elevation data: [GTOPO30](http://www.geonames.org/export/ws-overview.html)
* NACP Synergetic (SYNMAP) C<sub>4</sub> relative fraction: [Jung, M., Henkel, K., Herold, M., and Churkina, G.: Exploiting synergies of
global land cover products for carbon cycle modeling. Remote Sens. Environ.,
101, 534–553, doi:10.1016/j.rse.2006.01.020, 2006.](https://www.bgc-jena.mpg.de/bgi/uploads/Publ/Publications/Jung_et_al_2006.pdf).
* Global CO<sub>2</sub> data: [Ed Dlugokencky and Pieter Tans, NOAA/ESRL](www.esrl.noaa.gov/gmd/ccgg/trends/)

## Instructions

To generate all g1 fits from leaf, isotope and fluxnet data and manuscript figures, simple type:

```
$ make
```

## Dependancies

You will need a few python packages, namely `numpy`, `pandas`, `matplotlib` `scipy` and `lmfit`. The fitting code also exploits the python MPI library `multiprocessing` to speed things up. The R code depends on the `doby`, `multcomp`, `lme4` and `LMERConvenienceFunctions` packages.

## Contacts
- Belinda Medlyn: B.Medlyn at westernsydney.edu.au
- Martin De Kauwe: mdekauwe at gmail.com
