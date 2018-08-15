# dynamac
Dynamic simulation and testing for single-equation ARDL models in R and Stata. Link to the webpage can be found [here](https://andyphilips.github.io/dynamac/).

## Installation (R)
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
To install `dynamac` in Stata, you can download this [zip file](https://andyphilips.github.io/dynamac/Stata/dynamac.zip) and either call directly to the .ado files or place them in your "ado/plus/" folder. A SSC version will be added soon.
