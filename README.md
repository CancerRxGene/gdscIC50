# gdscIC50

The Genomics of Drugs Sensitivity in Cancer (GDSC) is one of the largest public resources of information for drug sensitivity in cancer cells and molecular markers of drug response. 

**gdscIC50** is an R package that fits dose response curves for data from the GDSC project. This package uses the non-linear mixed effects model as described in [(Vis, D.J. et al. Pharmacogenomics 2016, 17(7):691-700)](https://www.ncbi.nlm.nih.gov/pubmed/27180993). Once curves are fitted, the generated IC50 values can be used as input to analysis using [gdsctools](https://github.com/CancerRxGene/gdsctools). These two software libraries are used to generate the dose reponse data and statistical analysis presented on the website http://www.cancerrxgene.org.

You can install the package in R using the devtools library. If you haven’t already installed devtools run the following:
 ```
> install.packages("devtools")
```
Then install the gdscIC50 package. Building the vignette is recommended.
```
> devtools::install_github("cancerrxgene/gdscIC50", build_vignettes=TRUE)
```
Once installed, the vignette gives a guide to the curve fitting process.
```
> library(gdscIC50)
> vignette("gdscIC50")
```
