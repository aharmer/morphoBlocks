
Overview
--------

The 'morphoBlocks' package provides a workflow for constructing a multiple-part morphospace with regularised consensus principal component analysis (RCPCA) using either traditional landmarks or pseudolandmarks.

<br />  

<img src="https://user-images.githubusercontent.com/10540385/80056592-b22bb500-8578-11ea-8081-cfc300538eb3.png" width="600" style="display: block; margin: auto;" />


<br />  


Updates
-------

This is the original release of the package.


Installation
------------

#### morphoBlocks

You can install *morphoBlocks* directly within R using the *install\_github()* function from the [devtools](https://www.rstudio.com/products/rpackages/devtools/) package:

``` r
install_github("aharmer/morphoBlocks", build_vignettes = TRUE)
```

Depending on your setup, you may also need to install *Rtools* first. If you need *Rtools* you will get an error message during *morphoBlocks* installation. Just install *Rtools* then reinstall *bonespace*.

Alternatively, you can manually download and install the source package from the [morphoBlocks](https://morphoBlocks) webpage.

If you manually install *morphoBlocks*, the following packages are also required: *adephylo, data.table, geomorph, Morpho, phytools, RGCCA, Rvcg*.


How to use morphoBlocks
---------------------

