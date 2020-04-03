<img src="https://andyphilips.github.io/dynamac/img/logo.png" alt="drawing" width="200"/>

## Description
`dynamac` performs dynamic simulation and testing for single-equation ARDL models in R and Stata. Link to the webpage can be found [here](https://andyphilips.github.io/dynamac/).

## Installation (R)
[![Build Status](https://travis-ci.com/andyphilips/dynamac.svg?branch=master)](https://travis-ci.com/andyphilips/dynamac) [![Build status](https://ci.appveyor.com/api/projects/status/o8h5gdh5cuah359y?svg=true)](https://ci.appveyor.com/project/andyphilips/dynamac) ![](https://www.r-pkg.org/badges/version/dynamac) [![codecov](https://codecov.io/gh/andyphilips/dynamac/branch/master/graph/badge.svg)](https://codecov.io/gh/andyphilips/dynamac) [![metacran downloads](https://cranlogs.r-pkg.org/badges/dynamac)](https://cran.r-project.org/package=dynamac) [![doi](https://img.shields.io/badge/DOI:-10.32614%2FRJ--2018--076-blueviolet.svg)](https://doi.org/10.32614/RJ-2018-076)


`dynamac` in R can easily be installed from its CRAN repository:
```
install.packages("dynamac")
library(dynamac)
```

Alternatively, you can use the `devtools` package to load directly from GitHub:
```
library(devtools)
devtools::install_github("andyphilips/dynamac")
library(dynamac)
```
You should now be able to use the package.

## Installation (Stata)
To install `dynamac` in Stata type one of the following into the console:
```
net sj 18-4 st0545
```
```
findit dynamac
```
```
cap ado uninstall dynamac
cap ado uninstall dynardl
cap ado uninstall pssbounds
net install dynamac, from(https://github.com/andyphilips/dynamac/raw/gh-pages/Stata/src/)
```

Alternatively, you download the [raw files](https://github.com/andyphilips/dynamac/raw/gh-pages/Stata/src/) and either call directly to the .ado files or place them in your "ado/plus/" folder.

## I've found a bug!
Great! Please post it on the [issues page](https://github.com/andyphilips/dynamac/issues), or email the authors.
