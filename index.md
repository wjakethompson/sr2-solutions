--- 
knit: "bookdown::render_book"
title: "Statistical Rethinking: Solutions (2nd Edition)"
author: "Jake Thompson"
description: "This project contains solutions to exercises and homework for the second edition of Richard McElreath's Statistical Rethinking: A Bayesian Course Using R and Stan. Solutions are offered using the tidyverse and brms packages. Thus, this project should be viewed as a companion to both the original book and Solomon Kurz's translation of the text to tidyverse and brms syntax. Here, we extend the translation to also include solutions to the exercises and homework assignments associated with the course."
date: "2022-02-21"
url: 'https\://sr2-solutions.wjakethompson.com'
github-repo: wjakethompson/sr2-solutions
twitter-handle: wjakethompson
site: bookdown::bookdown_site
documentclass: book
bibliography: [bib/refs.bib, bib/pkgs.bib]
biblio-style: apa
csl: csl/apa.csl
link-citations: yes
---

# Welcome {-}

This project is a companion to Richard McElreath's *Statistical Rethinking (2nd Edition)* [@sr2], and the 2022 version of the accompanying course. For each of the 10 weeks, of the course (materials provided [here](https://github.com/rmcelreath/stat_rethinking_2022)), I work through the exercises in each chapter covered that week and the assigned homework problems. Like [Solomon Kurz](https://twitter.com/SolomonKurz), I take a {[tidyverse](https://tidyverse.tidyverse.org/)} [@R-tidyverse; @tidyverse2019] and {[brms](https://paul-buerkner.github.io/brms/)} [@R-brms; @brms2017; @brms2018] approach to solving these problems. For a translation of the actual book text to {tidyverse} and {brms} style code, please check out their project, [*Statistical rethinking* with brms, ggplot2, and the tidyverse: Second edition](https://bookdown.org/content/4857).

You can purchase *Statistical Rethinking: A Bayesian Course in R and Stan* [@sr2] from [CRC Press](https://www.routledge.com/Statistical-Rethinking-A-Bayesian-Course-with-Examples-in-R-and-STAN/McElreath/p/book/9780367139919?utm_source=crcpress.com&utm_medium=referral).

## Disclaimer {-}

This project is a work in progress. If you'd like to follow along, you can find the GitHub repository [here](https://github.com/wjakethompson/sr2-solutions). These solutions have not been checked by anybody, so there will undoubtedly be errors. If you find any, please contribute or let me know!

There are several ways to contribute. For simple edits or suggestions, you can use the "Edit this page" link in the sidebar at the right of the screen. You can also create a fork of the [repository](https://github.com/wjakethompson/sr2-solutions) and submit a pull request or [open an issue](https://github.com/wjakethompson/sr2-solutions/issues).

Please note that this project is released with a [Contributor Code of Conduct](https://github.com/wjakethompson/sr2-solutions/blob/main/CODE_OF_CONDUCT.md). By participating in this project you agree to abide by its terms.

## License {-}

<a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/"><img alt="Creative Commons Licence" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-sa/4.0/88x31.png" /></a>

This project is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/">Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License</a>.

## Colophon {-}

This project was built with:


```r
deps <- renv::dependencies() %>% 
  distinct(Package) %>% 
  arrange() %>% 
  pull()
#> Finding R package dependencies ... [8/18] [9/18] [10/18] [11/18] [12/18] [13/18] [14/18] [15/18] [16/18] [17/18] [18/18] Done!

sessioninfo::session_info(deps)
#> ─ Session info ───────────────────────────────────────────────────────────────
#>  setting  value
#>  version  R version 4.1.2 (2021-11-01)
#>  os       Ubuntu 20.04.3 LTS
#>  system   x86_64, linux-gnu
#>  ui       X11
#>  language (EN)
#>  collate  C.UTF-8
#>  ctype    C.UTF-8
#>  tz       UTC
#>  date     2022-02-21
#>  pandoc   2.14.2 @ /usr/bin/ (via rmarkdown)
#> 
#> ─ Packages ───────────────────────────────────────────────────────────────────
#>  ! package        * version    date (UTC) lib source
#>  P abind            1.4-5      2016-07-21 [?] CRAN (R 4.1.2)
#>  P arrayhelpers     1.1-0      2020-02-04 [?] CRAN (R 4.1.2)
#>    askpass          1.1        2019-01-13 [1] CRAN (R 4.1.2)
#>  P assertthat       0.2.1      2019-03-21 [?] CRAN (R 4.1.2)
#>  P backports        1.4.1      2021-12-13 [?] CRAN (R 4.1.2)
#>  P base64enc        0.1-3      2015-07-28 [?] CRAN (R 4.1.2)
#>  P bayesplot        1.8.1      2021-06-14 [?] CRAN (R 4.1.2)
#>    BH               1.78.0-0   2021-12-15 [1] CRAN (R 4.1.2)
#>    bit              4.0.4      2020-08-04 [1] CRAN (R 4.1.2)
#>    bit64            4.0.5      2020-08-30 [1] CRAN (R 4.1.2)
#>    bitops           1.0-7      2021-04-24 [1] CRAN (R 4.1.2)
#>    blob             1.2.2      2021-07-23 [1] CRAN (R 4.1.2)
#>  P bookdown         0.24       2021-09-02 [?] CRAN (R 4.1.2)
#>    boot             1.3-28     2021-05-03 [2] CRAN (R 4.1.2)
#>    brew             1.0-6      2011-04-13 [1] CRAN (R 4.1.2)
#>  P bridgesampling   1.1-2      2021-04-16 [?] CRAN (R 4.1.2)
#>    brio             1.1.3      2021-11-30 [1] CRAN (R 4.1.2)
#>  P brms           * 2.16.3     2021-11-22 [?] CRAN (R 4.1.2)
#>  P Brobdingnag      1.2-6      2018-08-13 [?] CRAN (R 4.1.2)
#>  P broom            0.7.12     2022-01-28 [?] CRAN (R 4.1.2)
#>  P bslib            0.3.1      2021-10-06 [?] CRAN (R 4.1.2)
#>  P cachem           1.0.6      2021-08-19 [?] CRAN (R 4.1.2)
#>  P callr            3.7.0      2021-04-20 [?] CRAN (R 4.1.2)
#>  P cellranger       1.1.0      2016-07-27 [?] CRAN (R 4.1.2)
#>  P checkmate        2.0.0      2020-02-06 [?] CRAN (R 4.1.2)
#>    class            7.3-19     2021-05-03 [2] CRAN (R 4.1.2)
#>    classInt         0.4-3      2020-04-07 [1] CRAN (R 4.1.2)
#>  P cli              3.1.1.9000 2022-02-21 [?] Github (r-lib/cli@a1d4569)
#>    clipr            0.7.1      2020-10-08 [1] CRAN (R 4.1.2)
#>  P cmdstanr       * 0.4.0.9001 2022-02-21 [?] Github (stan-dev/cmdstanr@a2a97d9)
#>  P coda             0.19-4     2020-09-30 [?] CRAN (R 4.1.2)
#>    codetools        0.2-18     2020-11-04 [2] CRAN (R 4.1.2)
#>  P colorspace       2.0-2      2021-06-24 [?] CRAN (R 4.1.2)
#>  P colourpicker     1.1.1      2021-10-04 [?] CRAN (R 4.1.2)
#>    commonmark       1.7        2018-12-01 [1] CRAN (R 4.1.2)
#>    cpp11            0.4.2      2021-11-30 [1] CRAN (R 4.1.2)
#>  P crayon           1.4.2      2021-10-29 [?] CRAN (R 4.1.2)
#>    credentials      1.3.2      2021-11-29 [1] CRAN (R 4.1.2)
#>  P crosstalk        1.2.0      2021-11-04 [?] CRAN (R 4.1.2)
#>  P curl             4.3.2      2021-06-23 [?] CRAN (R 4.1.2)
#>  P dagitty        * 0.3-1      2021-01-21 [?] CRAN (R 4.1.2)
#>    data.table       1.14.2     2021-09-27 [1] CRAN (R 4.1.2)
#>  P DBI              1.1.2      2021-12-20 [?] CRAN (R 4.1.2)
#>  P dbplyr           2.1.1      2021-04-06 [?] CRAN (R 4.1.2)
#>    desc             1.4.0      2021-09-28 [1] CRAN (R 4.1.2)
#>    devtools         2.4.3      2021-11-30 [1] CRAN (R 4.1.2)
#>    diffobj          0.3.5      2021-10-05 [1] CRAN (R 4.1.2)
#>  P digest           0.6.29     2021-12-01 [?] CRAN (R 4.1.2)
#>  P distributional   0.3.0      2022-01-05 [?] CRAN (R 4.1.2)
#>  P downlit          0.4.0.9000 2022-02-21 [?] Github (r-lib/downlit@4011b2f)
#>  P dplyr          * 1.0.7      2021-06-18 [?] CRAN (R 4.1.2)
#>  P DT               0.20       2021-11-15 [?] CRAN (R 4.1.2)
#>    dtplyr           1.2.1      2022-01-19 [1] CRAN (R 4.1.2)
#>  P dygraphs         1.1.1.6    2018-07-11 [?] CRAN (R 4.1.2)
#>    e1071            1.7-9      2021-09-16 [1] CRAN (R 4.1.2)
#>  P ellipsis         0.3.2      2021-04-29 [?] CRAN (R 4.1.2)
#>    emo              0.0.0.9000 2022-02-21 [1] Github (hadley/emo@3f03b11)
#>  P evaluate         0.14       2019-05-28 [?] CRAN (R 4.1.2)
#>  P extrafont        0.17       2014-12-08 [?] CRAN (R 4.1.2)
#>  P extrafontdb      1.0        2012-06-11 [?] CRAN (R 4.1.2)
#>  P fansi            1.0.2      2022-01-14 [?] CRAN (R 4.1.2)
#>  P farver           2.1.0      2021-02-28 [?] CRAN (R 4.1.2)
#>  P fastmap          1.1.0      2021-01-25 [?] CRAN (R 4.1.2)
#>    fontawesome      0.2.2      2021-07-02 [1] CRAN (R 4.1.2)
#>  P forcats        * 0.5.1      2021-01-27 [?] CRAN (R 4.1.2)
#>  P fs               1.5.2      2021-12-08 [?] CRAN (R 4.1.2)
#>    future           1.23.0     2021-10-31 [1] CRAN (R 4.1.2)
#>    gargle           1.2.0      2021-07-02 [1] CRAN (R 4.1.2)
#>  P gdtools          0.2.3      2021-01-06 [?] CRAN (R 4.1.2)
#>  P generics         0.1.1      2021-10-25 [?] CRAN (R 4.1.2)
#>  P geomtextpath   * 0.1.0.9000 2022-02-21 [?] Github (allancameron/geomtextpath@0612ba8)
#>    gert             1.5.0      2022-01-03 [1] CRAN (R 4.1.2)
#>    gganimate        1.0.7      2020-10-15 [1] CRAN (R 4.1.2)
#>    ggdag            0.2.4      2021-10-10 [1] CRAN (R 4.1.2)
#>  P ggdist         * 3.0.1      2021-11-30 [?] CRAN (R 4.1.2)
#>    ggforce          0.3.3      2021-03-05 [1] CRAN (R 4.1.2)
#>    gghighlight      0.3.2      2021-06-05 [1] CRAN (R 4.1.2)
#>  P ggplot2        * 3.3.5      2021-06-25 [?] CRAN (R 4.1.2)
#>    ggraph           2.0.5      2021-02-23 [1] CRAN (R 4.1.2)
#>  P ggrepel        * 0.9.1      2021-01-15 [?] CRAN (R 4.1.2)
#>  P ggridges       * 0.5.3      2021-01-08 [?] CRAN (R 4.1.2)
#>  P ggtext         * 0.1.1      2020-12-17 [?] CRAN (R 4.1.2)
#>    gh               1.3.0      2021-04-30 [1] CRAN (R 4.1.2)
#>    gifski           1.4.3-1    2021-05-02 [1] CRAN (R 4.1.2)
#>    gitcreds         0.1.1      2020-12-04 [1] CRAN (R 4.1.2)
#>    globals          0.14.0     2020-11-22 [1] CRAN (R 4.1.2)
#>  P glue           * 1.6.1      2022-01-22 [?] CRAN (R 4.1.2)
#>    googledrive      2.0.0      2021-07-08 [1] CRAN (R 4.1.2)
#>    googlesheets4    1.0.0      2021-07-21 [1] CRAN (R 4.1.2)
#>    graphlayouts     0.8.0      2022-01-03 [1] CRAN (R 4.1.2)
#>  P gridExtra        2.3        2017-09-09 [?] CRAN (R 4.1.2)
#>  P gridtext         0.1.4      2020-12-10 [?] CRAN (R 4.1.2)
#>    gt               0.3.1      2021-08-07 [1] CRAN (R 4.1.2)
#>  P gtable           0.3.0      2019-03-25 [?] CRAN (R 4.1.2)
#>  P gtools           3.9.2      2021-06-06 [?] CRAN (R 4.1.2)
#>  P haven            2.4.3      2021-08-04 [?] CRAN (R 4.1.2)
#>    HDInterval       0.2.2      2020-05-23 [1] CRAN (R 4.1.2)
#>  P here           * 1.0.1      2020-12-13 [?] CRAN (R 4.1.2)
#>    highr            0.9        2021-04-16 [1] CRAN (R 4.1.2)
#>  P hms              1.1.1      2021-09-26 [?] CRAN (R 4.1.2)
#>  P hrbrthemes       0.8.0      2020-03-06 [?] CRAN (R 4.1.2)
#>  P htmltools        0.5.2      2021-08-25 [?] CRAN (R 4.1.2)
#>  P htmlwidgets      1.5.4      2021-09-08 [?] CRAN (R 4.1.2)
#>  P httpuv           1.6.5      2022-01-05 [?] CRAN (R 4.1.2)
#>  P httr             1.4.2      2020-07-20 [?] CRAN (R 4.1.2)
#>    ids              1.0.1      2017-05-31 [1] CRAN (R 4.1.2)
#>  P igraph           1.2.11     2022-01-04 [?] CRAN (R 4.1.2)
#>    ini              0.3.1      2018-05-20 [1] CRAN (R 4.1.2)
#>  P inline           0.3.19     2021-05-31 [?] CRAN (R 4.1.2)
#>    isoband          0.2.5      2021-07-13 [1] CRAN (R 4.1.2)
#>    jpeg             0.1-9      2021-07-24 [1] CRAN (R 4.1.2)
#>  P jquerylib        0.1.4      2021-04-26 [?] CRAN (R 4.1.2)
#>  P jsonlite         1.7.3      2022-01-17 [?] CRAN (R 4.1.2)
#>    kableExtra       1.3.4      2021-02-20 [1] CRAN (R 4.1.2)
#>    KernSmooth       2.23-20    2021-05-03 [2] CRAN (R 4.1.2)
#>  P knitr            1.37       2021-12-16 [?] CRAN (R 4.1.2)
#>    labeling         0.4.2      2020-10-20 [1] CRAN (R 4.1.2)
#>  P later            1.3.0      2021-08-18 [?] CRAN (R 4.1.2)
#>    lattice          0.20-45    2021-09-22 [2] CRAN (R 4.1.2)
#>    lazyeval         0.2.2      2019-03-15 [1] CRAN (R 4.1.2)
#>  P lifecycle        1.0.1.9000 2022-02-21 [?] Github (r-lib/lifecycle@56eafa4)
#>    listenv          0.8.0      2019-12-05 [1] CRAN (R 4.1.2)
#>  P loo            * 2.4.1      2020-12-09 [?] CRAN (R 4.1.2)
#>    lpSolve          5.6.15     2020-01-24 [1] CRAN (R 4.1.2)
#>  P lubridate        1.8.0      2021-10-07 [?] CRAN (R 4.1.2)
#>  P magrittr         2.0.2      2022-01-26 [?] CRAN (R 4.1.2)
#>  P markdown         1.1        2019-08-07 [?] CRAN (R 4.1.2)
#>    MASS             7.3-54     2021-05-03 [2] CRAN (R 4.1.2)
#>    Matrix           1.3-4      2021-06-01 [2] CRAN (R 4.1.2)
#>  P matrixStats      0.61.0     2021-09-17 [?] CRAN (R 4.1.2)
#>  P memoise          2.0.1      2021-11-26 [?] CRAN (R 4.1.2)
#>    mgcv             1.8-38     2021-10-06 [2] CRAN (R 4.1.2)
#>  P mime             0.12       2021-09-28 [?] CRAN (R 4.1.2)
#>  P miniUI           0.1.1.1    2018-05-18 [?] CRAN (R 4.1.2)
#>  P modelr           0.1.8      2020-05-19 [?] CRAN (R 4.1.2)
#>  P munsell          0.5.0      2018-06-12 [?] CRAN (R 4.1.2)
#>  P mvtnorm          1.1-3      2021-10-08 [?] CRAN (R 4.1.2)
#>    nleqslv          3.3.2      2018-05-17 [1] CRAN (R 4.1.2)
#>    nlme             3.1-153    2021-09-07 [2] CRAN (R 4.1.2)
#>    numDeriv         2016.8-1.1 2019-06-06 [1] CRAN (R 4.1.2)
#>    officedown       0.2.3      2021-11-16 [1] CRAN (R 4.1.2)
#>    officer          0.4.1      2021-11-14 [1] CRAN (R 4.1.2)
#>    openssl          1.4.6      2021-12-19 [1] CRAN (R 4.1.2)
#>    packrat          0.7.0      2021-08-20 [1] CRAN (R 4.1.2)
#>    parallelly       1.30.0     2021-12-17 [1] CRAN (R 4.1.2)
#>    patchwork        1.1.1      2020-12-17 [1] CRAN (R 4.1.2)
#>  P pillar           1.6.5      2022-01-25 [?] CRAN (R 4.1.2)
#>  P pkgbuild         1.3.1      2021-12-20 [?] CRAN (R 4.1.2)
#>  P pkgconfig        2.0.3      2019-09-22 [?] CRAN (R 4.1.2)
#>    pkgload          1.2.4      2021-11-30 [1] CRAN (R 4.1.2)
#>  P plyr             1.8.6      2020-03-03 [?] CRAN (R 4.1.2)
#>    png              0.1-7      2013-12-03 [1] CRAN (R 4.1.2)
#>    polyclip         1.10-0     2019-03-14 [1] CRAN (R 4.1.2)
#>  P posterior        1.2.0      2022-01-05 [?] CRAN (R 4.1.2)
#>    praise           1.0.0      2015-08-11 [1] CRAN (R 4.1.2)
#>  P prettyunits      1.1.1      2020-01-24 [?] CRAN (R 4.1.2)
#>  P processx         3.5.2      2021-04-30 [?] CRAN (R 4.1.2)
#>    progress         1.2.2      2019-05-16 [1] CRAN (R 4.1.2)
#>  P promises         1.2.0.1    2021-02-11 [?] CRAN (R 4.1.2)
#>    prompt           1.0.1      2022-02-21 [1] Github (gaborcsardi/prompt@7ef0f2e)
#>    proxy            0.4-26     2021-06-07 [1] CRAN (R 4.1.2)
#>  P ps               1.6.0      2021-02-28 [?] CRAN (R 4.1.2)
#>  P purrr          * 0.3.4      2020-04-17 [?] CRAN (R 4.1.2)
#>  R R                <NA>       <NA>       [?] <NA>
#>  P R6               2.5.1      2021-08-19 [?] CRAN (R 4.1.2)
#>  P ragg             1.2.1      2021-12-06 [?] CRAN (R 4.1.2)
#>    rappdirs         0.3.3      2021-01-31 [1] CRAN (R 4.1.2)
#>  P ratlas           0.0.0.9000 2022-02-21 [?] Github (atlas-aai/ratlas@ebb795b)
#>    rcmdcheck        1.4.0      2021-09-27 [1] CRAN (R 4.1.2)
#>    RColorBrewer     1.1-2      2014-12-07 [1] CRAN (R 4.1.2)
#>  P Rcpp           * 1.0.8      2022-01-13 [?] CRAN (R 4.1.2)
#>    RcppArmadillo    0.10.8.1.0 2022-01-24 [1] CRAN (R 4.1.2)
#>    RcppEigen        0.3.3.9.1  2020-12-17 [1] CRAN (R 4.1.2)
#>  P RcppParallel     5.1.5      2022-01-05 [?] CRAN (R 4.1.2)
#>    RCurl            1.98-1.5   2021-09-17 [1] CRAN (R 4.1.2)
#>  P readr          * 2.1.1      2021-11-30 [?] CRAN (R 4.1.2)
#>  P readxl           1.3.1      2019-03-13 [?] CRAN (R 4.1.2)
#>    rematch          1.0.1      2016-04-21 [1] CRAN (R 4.1.2)
#>    rematch2         2.1.2      2020-05-01 [1] CRAN (R 4.1.2)
#>    remotes          2.4.2      2021-11-30 [1] CRAN (R 4.1.2)
#>    renv             0.15.2     2022-01-24 [1] CRAN (R 4.1.2)
#>  P reprex           2.0.1      2021-08-05 [?] CRAN (R 4.1.2)
#>  P reshape2         1.4.4      2020-04-09 [?] CRAN (R 4.1.2)
#>  P rethinking     * 2.21       2022-02-21 [?] Github (rmcelreath/rethinking@783d111)
#>    rJava            1.0-6      2021-12-10 [1] CRAN (R 4.1.2)
#>  P rlang            1.0.0      2022-01-26 [?] CRAN (R 4.1.2)
#>  P rmarkdown        2.11       2021-09-14 [?] CRAN (R 4.1.2)
#>    roxygen2         7.1.2      2021-09-08 [1] CRAN (R 4.1.2)
#>  P rprojroot        2.0.2      2020-11-15 [?] CRAN (R 4.1.2)
#>  P rsconnect        0.8.25     2021-11-19 [?] CRAN (R 4.1.2)
#>  P rstan          * 2.21.3     2021-12-19 [?] CRAN (R 4.1.2)
#>  P rstantools       2.1.1      2020-07-06 [?] CRAN (R 4.1.2)
#>    rstudioapi       0.13       2020-11-12 [1] CRAN (R 4.1.2)
#>  P Rttf2pt1         1.3.9      2021-07-22 [?] CRAN (R 4.1.2)
#>    rversions        2.1.1      2021-05-31 [1] CRAN (R 4.1.2)
#>  P rvest            1.0.2      2021-10-16 [?] CRAN (R 4.1.2)
#>    rvg              0.2.5      2020-06-30 [1] CRAN (R 4.1.2)
#>    s2               1.0.7      2021-09-28 [1] CRAN (R 4.1.2)
#>  P sass             0.4.0.9000 2022-02-21 [?] Github (rstudio/sass@f7a9540)
#>  P scales           1.1.1      2020-05-11 [?] CRAN (R 4.1.2)
#>    selectr          0.4-2      2019-11-20 [1] CRAN (R 4.1.2)
#>    servr            0.24       2021-11-16 [1] CRAN (R 4.1.2)
#>  P sessioninfo      1.2.2      2021-12-06 [?] CRAN (R 4.1.2)
#>    sf               1.0-6      2022-02-04 [1] CRAN (R 4.1.2)
#>  P shape            1.4.6      2021-05-19 [?] CRAN (R 4.1.2)
#>  P shiny            1.7.1      2021-10-02 [?] CRAN (R 4.1.2)
#>  P shinyjs          2.1.0      2021-12-23 [?] CRAN (R 4.1.2)
#>  P shinystan        2.5.0      2018-05-01 [?] CRAN (R 4.1.2)
#>  P shinythemes      1.2.0      2021-01-25 [?] CRAN (R 4.1.2)
#>    sourcetools      0.1.7      2018-04-25 [1] CRAN (R 4.1.2)
#>  P StanHeaders    * 2.21.0-7   2020-12-17 [?] CRAN (R 4.1.2)
#>    staplr           3.1.1      2021-01-11 [1] CRAN (R 4.1.2)
#>  P stringi          1.7.6      2021-11-29 [?] CRAN (R 4.1.2)
#>  P stringr        * 1.4.0      2019-02-10 [?] CRAN (R 4.1.2)
#>    svglite          2.0.0      2021-02-20 [1] CRAN (R 4.1.2)
#>  P svUnit           1.0.6      2021-04-19 [?] CRAN (R 4.1.2)
#>    sys              3.4        2020-07-23 [1] CRAN (R 4.1.2)
#>  P systemfonts      1.0.3      2021-10-13 [?] CRAN (R 4.1.2)
#>  P tensorA          0.36.2     2020-11-19 [?] CRAN (R 4.1.2)
#>    testthat         3.1.2      2022-01-20 [1] CRAN (R 4.1.2)
#>  P textshaping      0.3.6      2021-10-13 [?] CRAN (R 4.1.2)
#>  P threejs          0.3.3      2020-01-21 [?] CRAN (R 4.1.2)
#>  P tibble         * 3.1.6      2021-11-07 [?] CRAN (R 4.1.2)
#>  P tidybayes      * 3.0.2      2022-01-05 [?] CRAN (R 4.1.2)
#>    tidygraph        1.2.0      2020-05-12 [1] CRAN (R 4.1.2)
#>  P tidyr          * 1.1.4      2021-09-27 [?] CRAN (R 4.1.2)
#>  P tidyselect       1.1.1      2021-04-30 [?] CRAN (R 4.1.2)
#>  P tidyverse      * 1.3.1      2021-04-15 [?] CRAN (R 4.1.2)
#>    tinytex          0.36       2021-12-19 [1] CRAN (R 4.1.2)
#>    transformr       0.1.3      2020-07-05 [1] CRAN (R 4.1.2)
#>    tweenr           1.0.2      2021-03-23 [1] CRAN (R 4.1.2)
#>  P tzdb             0.2.0      2021-10-27 [?] CRAN (R 4.1.2)
#>    units            0.7-2      2021-06-08 [1] CRAN (R 4.1.2)
#>    usethis          2.1.5      2021-12-09 [1] CRAN (R 4.1.2)
#>  P utf8             1.2.2      2021-07-24 [?] CRAN (R 4.1.2)
#>    uuid             1.0-3      2021-11-01 [1] CRAN (R 4.1.2)
#>  P V8               4.0.0      2021-12-23 [?] CRAN (R 4.1.2)
#>  P vctrs            0.3.8      2021-04-29 [?] CRAN (R 4.1.2)
#>    viridis          0.6.2      2021-10-13 [1] CRAN (R 4.1.2)
#>    viridisLite      0.4.0      2021-04-13 [1] CRAN (R 4.1.2)
#>    vroom            1.5.7      2021-11-30 [1] CRAN (R 4.1.2)
#>    waldo            0.3.1      2021-09-14 [1] CRAN (R 4.1.2)
#>    webshot          0.5.2      2019-11-22 [1] CRAN (R 4.1.2)
#>    whisker          0.4        2019-08-28 [1] CRAN (R 4.1.2)
#>  P withr            2.4.3      2021-11-30 [?] CRAN (R 4.1.2)
#>  P wjake          * 0.1.0      2022-02-21 [?] Github (wjakethompson/wjake@d3f00de)
#>    wk               0.6.0      2022-01-03 [1] CRAN (R 4.1.2)
#>    xaringan         0.22       2021-06-23 [1] CRAN (R 4.1.2)
#>  P xfun             0.29       2021-12-14 [?] CRAN (R 4.1.2)
#>    XML              3.99-0.8   2021-09-17 [1] CRAN (R 4.1.2)
#>  P xml2             1.3.3      2021-11-30 [?] CRAN (R 4.1.2)
#>    xopen            1.0.0      2018-09-17 [1] CRAN (R 4.1.2)
#>  P xtable           1.8-4      2019-04-21 [?] CRAN (R 4.1.2)
#>  P xts              0.12.1     2020-09-09 [?] CRAN (R 4.1.2)
#>  P yaml             2.2.2      2022-01-25 [?] CRAN (R 4.1.2)
#>    zip              2.2.0      2021-05-31 [1] CRAN (R 4.1.2)
#>  P zoo              1.8-9      2021-03-09 [?] CRAN (R 4.1.2)
#> 
#>  [1] /home/runner/.cache/R/renv/library/sr2-solutions-9e1c7ff5/R-4.1/x86_64-pc-linux-gnu
#>  [2] /opt/R/4.1.2/lib/R/library
#> 
#>  P ── Loaded and on-disk path mismatch.
#>  R ── Package was removed from disk.
#> 
#> ──────────────────────────────────────────────────────────────────────────────
```
