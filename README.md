# SpectreMAP
An extension of the 'Spectre' package (https://github.com/ImmuneDynamics/Spectre) for spatial analysis of high-dimensional imaging data.

### Current version
`v0.2.2`

### About
This package expands on the functionality of 'Spectre' to enable spatial analysis of high-dimensional imaging data, such as that generated by Imaging Mass Cytometry (IMC).

### Instructions
Spatial analysis workflows are listed on the Spectre Home Page: https://wiki.centenary.org.au/display/SPECTRE. 

### Installing Spectre
In R, install and load the 'devtools' library.

```     
if(!require('devtools')) {install.packages('devtools')}
library('devtools')
```

Subsequently, use the 'install_github' function to install and load the Spectre and SpectreMAP packages. By default this will load the 'master' branch, which is the same as the latest stable release version. To install a specific release version, see https://cran.r-project.org/web/packages/githubinstall/vignettes/githubinstall.html.

```
    install_github("immunedynamics/spectre")
    install_github("immunedynamics/SpectreMAP")
```

Next install the additional packages.

```
### Install additional packages
    ## Install BiocManager to download packages from Bioconductor
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
 
    ## Download additional BioConductor packages
    if(!require('flowCore')) {BiocManager::install('flowCore')}
    if(!require('Biobase')) {BiocManager::install('Biobase')}
    if(!require('flowViz')) {BiocManager::install('flowViz')}
    if(!require('FlowSOM')) {BiocManager::install('FlowSOM')}
 
    ## Velox for fast extraction
    install_github("hunzikp/velox")
```

Load packages to ensure installation has been successful.

```
### Load packages
 
    library(Spectre)
    library(SpectreMAP)
 
    package.check(type = 'spatial')
    package.load(type = 'spatial')
```
