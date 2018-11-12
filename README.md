# Body size, reef area, and temperature predict global reef-fish species richness across spatial scales

This repository contains code and data needed to reproduce the article:

**Barneche DR, Rezende EL, Parravicini V, Maire E, Edgar GJ, Stuart-Smith RD, Arias-González JE, Ferreira CEL, Friedlander AM, Green AL, Luiz OJ, Rodríguez-Zaragoza FA, Vigliola L, Kulbicki M, Floeter SR**, Body size, reef area, and temperature predict global reef-fish species richness across spatial scales. *Global Ecology and Biogeography*. 10.1111/geb.12851

## Instructions

All analyses were done in `R`. To compile the figures, analyses and tables we use the [remake](https://github.com/richfitz/remake) package for R. The `remake` package depends on `storr`, so install it first like this:

```r
devtools::install_github('richfitz/storr', dependencies=TRUE)
```

Run `install.packages('devtools')` to install `devtools` first if needed. Then you can install remake using the `devtools` package:

```r
devtools::install_github('richfitz/remake', dependencies=TRUE)
```

Next you need to open an R session with working directory set to the root of the project.

We use a number of packages, missing packages can be easily installed by remake:

```r
remake::install_missing_packages()
```

And then install the package `fontcm`, via `extrafont`. This installs the font `CM Roman` we use in our figures (for more information on see [these instructions](https://cran.r-project.org/web/packages/fontcm/README.html):

```r
extrafont::font_install('fontcm')
```

Then, to generate all figures, analyses, and tables simply do:

```r
remake::make()
```

All output will be automatically placed in a directory called `output` (it is going to be automatically created for you).

Also notice that the Bayesian analysis in this paper might take up to 2 hours to run on a regular computer. If you find remake confusing and prefer to run plain R, you can use remake to build a script `build.R` that produces a given output, e.g.

```r
remake::make_script(filename = 'build.R')
```

### The paper can be reproduced using the following software and associated packages:
```
R version 3.4.3 (2017-11-30)
Platform: x86_64-apple-darwin15.6.0 (64-bit)
Running under: macOS High Sierra 10.13.6

Matrix products: default
BLAS: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRblas.0.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRlapack.dylib

locale:
[1] en_AU.UTF-8/en_AU.UTF-8/en_AU.UTF-8/C/en_AU.UTF-8/en_AU.UTF-8

attached base packages:
[1] tools     parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] Hmisc_4.1-1        Formula_1.2-2      survival_2.41-3    mapproj_1.2-5      maps_3.2.0         fontcm_1.1        
 [7] extrafont_0.17     LoLinR_0.0.0.9000  vegan_2.4-6        lattice_0.20-35    permute_0.9-4      ape_5.0           
[13] brms_2.6.0         Rcpp_0.12.19       iNEXT_2.0.12       plyr_1.8.4         rstan_2.18.2       StanHeaders_2.18.0
[19] ggplot2_3.0.0     

loaded via a namespace (and not attached):
 [1] nlme_3.1-131         matrixStats_0.53.0   xts_0.11-1           RColorBrewer_1.1-2   threejs_0.3.1       
 [6] backports_1.1.2      R6_2.3.0             DT_0.4               rpart_4.1-11         lazyeval_0.2.1      
[11] mgcv_1.8-22          colorspace_1.3-2     nnet_7.3-12          withr_2.1.2          tidyselect_0.2.4    
[16] gridExtra_2.3        processx_3.1.0       Brobdingnag_1.2-4    compiler_3.4.3       extrafontdb_1.0     
[21] htmlTable_1.11.2     shinyjs_1.0          colourpicker_1.0     checkmate_1.8.5      scales_1.0.0        
[26] dygraphs_1.1.1.6     lmtest_0.9-36        mvtnorm_1.0-7        ggridges_0.5.1       callr_2.0.4         
[31] stringr_1.3.1        digest_0.6.18        foreign_0.8-69       base64enc_0.1-3      pkgconfig_2.0.2     
[36] htmltools_0.3.6      htmlwidgets_1.3      rlang_0.2.2          rstudioapi_0.8       shiny_1.1.0         
[41] bindr_0.1.1          zoo_1.8-4            crosstalk_1.0.0      gtools_3.8.1         acepack_1.4.1       
[46] dplyr_0.7.6          inline_0.3.15        magrittr_1.5         loo_2.0.0            bayesplot_1.6.0     
[51] Matrix_1.2-14        munsell_0.5.0        abind_1.4-5          stringi_1.2.4        yaml_2.2.0          
[56] MASS_7.3-47          storr_1.1.3          pkgbuild_1.0.0       grid_3.4.3           promises_1.0.1      
[61] crayon_1.3.4         miniUI_0.1.1         splines_3.4.3        knitr_1.20           pillar_1.1.0        
[66] igraph_1.1.2         markdown_0.8         shinystan_2.5.0      codetools_0.2-15     reshape2_1.4.3      
[71] stats4_3.4.3         rstantools_1.5.1     glue_1.3.0           latticeExtra_0.6-28  data.table_1.10.4-3 
[76] httpuv_1.4.5         Rttf2pt1_1.3.5       gtable_0.2.0         purrr_0.2.5          assertthat_0.2.0    
[81] mime_0.6             xtable_1.8-3         coda_0.19-1          later_0.7.5          rsconnect_0.8.8     
[86] tibble_1.4.2         shinythemes_1.1.1    bindrcpp_0.2.2       cluster_2.0.6        remake_0.3.0        
[91] bridgesampling_0.4-0
```

## How to download this project for people not familiar with GitHub:  
* on the project main page on GitHub, click on the green button `clone or download` and then click on `Download ZIP`  

## Bug reporting
* Please [report any issues or bugs](https://github.com/dbarneche/speciespackinggeb/issues).
