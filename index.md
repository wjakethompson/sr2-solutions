--- 
knit: "bookdown::render_book"
title: "Statistical Rethinking: Solutions (2nd Edition)"
author: "Jake Thompson"
description: "This project contains solutions to exercises and homework for the second edition of Richard McElreath's Statistical Rethinking: A Bayesian Course Using R and Stan. Solutions are offered using the tidyverse and brms packages. Thus, this project should be viewed as a companion to both the original book and Solomon Kurz's translation of the text to tidyverse and brms syntax. Here, we extend the translation to also include solutions to the exercises and homework assignments associated with the course."
date: "2022-01-14"
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
desc <- desc::description$new()
deps <- desc$get_deps() %>% 
  filter(package != "R") %>% 
  pull(package)

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
#>  date     2022-01-14
#>  pandoc   2.14.2 @ /usr/bin/ (via rmarkdown)
#> 
#> ─ Packages ───────────────────────────────────────────────────────────────────
#>  package              * version    date (UTC) lib source
#>  abind                  1.4-5      2016-07-21 [1] CRAN (R 4.1.2)
#>  arrayhelpers           1.1-0      2020-02-04 [1] CRAN (R 4.1.2)
#>  askpass                1.1        2019-01-13 [1] CRAN (R 4.1.2)
#>  assertthat             0.2.1      2019-03-21 [1] CRAN (R 4.1.2)
#>  backports              1.4.1      2021-12-13 [1] CRAN (R 4.1.2)
#>  base64enc              0.1-3      2015-07-28 [1] CRAN (R 4.1.2)
#>  bayesplot              1.8.1      2021-06-14 [1] CRAN (R 4.1.2)
#>  BH                     1.78.0-0   2021-12-15 [1] CRAN (R 4.1.2)
#>  bit                    4.0.4      2020-08-04 [1] CRAN (R 4.1.2)
#>  bit64                  4.0.5      2020-08-30 [1] CRAN (R 4.1.2)
#>  blob                   1.2.2      2021-07-23 [1] CRAN (R 4.1.2)
#>  bookdown               0.24       2021-09-02 [1] CRAN (R 4.1.2)
#>  boot                   1.3-28     2021-05-03 [2] CRAN (R 4.1.2)
#>  bridgesampling         1.1-2      2021-04-16 [1] CRAN (R 4.1.2)
#>  brio                   1.1.3      2021-11-30 [1] CRAN (R 4.1.2)
#>  brms                 * 2.16.3     2021-11-22 [1] CRAN (R 4.1.2)
#>  Brobdingnag            1.2-6      2018-08-13 [1] CRAN (R 4.1.2)
#>  broom                  0.7.11     2022-01-03 [1] CRAN (R 4.1.2)
#>  bslib                  0.3.1      2021-10-06 [1] CRAN (R 4.1.2)
#>  cachem                 1.0.6      2021-08-19 [1] CRAN (R 4.1.2)
#>  callr                  3.7.0      2021-04-20 [1] CRAN (R 4.1.2)
#>  cellranger             1.1.0      2016-07-27 [1] CRAN (R 4.1.2)
#>  checkmate              2.0.0      2020-02-06 [1] CRAN (R 4.1.2)
#>  cli                    3.1.0      2021-10-27 [1] CRAN (R 4.1.2)
#>  clipr                  0.7.1      2020-10-08 [1] CRAN (R 4.1.2)
#>  cmdstanr             * 0.4.0.9001 2022-01-14 [1] Github (stan-dev/cmdstanr@a2a97d9)
#>  coda                   0.19-4     2020-09-30 [1] CRAN (R 4.1.2)
#>  codetools              0.2-18     2020-11-04 [2] CRAN (R 4.1.2)
#>  colorblindr            0.1.0      2022-01-06 [1] Github (clauswilke/colorblindr@e6730be)
#>  colorspace             2.0-2      2021-06-24 [1] CRAN (R 4.1.2)
#>  colourpicker           1.1.1      2021-10-04 [1] CRAN (R 4.1.2)
#>  commonmark             1.7        2018-12-01 [1] CRAN (R 4.1.2)
#>  cowplot                1.1.1      2020-12-30 [1] CRAN (R 4.1.2)
#>  cpp11                  0.4.2      2021-11-30 [1] CRAN (R 4.1.2)
#>  crayon                 1.4.2      2021-10-29 [1] CRAN (R 4.1.2)
#>  crosstalk              1.2.0      2021-11-04 [1] CRAN (R 4.1.2)
#>  curl                   4.3.2      2021-06-23 [1] CRAN (R 4.1.2)
#>  dagitty                0.3-1      2021-01-21 [1] CRAN (R 4.1.2)
#>  data.table             1.14.2     2021-09-27 [1] CRAN (R 4.1.2)
#>  DBI                    1.1.2      2021-12-20 [1] CRAN (R 4.1.2)
#>  dbplyr                 2.1.1      2021-04-06 [1] CRAN (R 4.1.2)
#>  desc                   1.4.0      2021-09-28 [1] CRAN (R 4.1.2)
#>  digest                 0.6.29     2021-12-01 [1] CRAN (R 4.1.2)
#>  distributional         0.3.0      2022-01-05 [1] CRAN (R 4.1.2)
#>  downlit                0.4.0      2021-10-29 [1] CRAN (R 4.1.2)
#>  dplyr                * 1.0.7      2021-06-18 [1] CRAN (R 4.1.2)
#>  DT                     0.20       2021-11-15 [1] CRAN (R 4.1.2)
#>  dtplyr                 1.2.0      2021-12-05 [1] CRAN (R 4.1.2)
#>  dygraphs               1.1.1.6    2018-07-11 [1] CRAN (R 4.1.2)
#>  ellipsis               0.3.2      2021-04-29 [1] CRAN (R 4.1.2)
#>  emo                    0.0.0.9000 2022-01-06 [1] Github (hadley/emo@3f03b11)
#>  evaluate               0.14       2019-05-28 [1] CRAN (R 4.1.2)
#>  extrafont              0.17       2014-12-08 [1] CRAN (R 4.1.2)
#>  extrafontdb            1.0        2012-06-11 [1] CRAN (R 4.1.2)
#>  fansi                  1.0.0      2022-01-10 [1] CRAN (R 4.1.2)
#>  farver                 2.1.0      2021-02-28 [1] CRAN (R 4.1.2)
#>  fastmap                1.1.0      2021-01-25 [1] CRAN (R 4.1.2)
#>  fontawesome            0.2.2      2021-07-02 [1] CRAN (R 4.1.2)
#>  forcats              * 0.5.1      2021-01-27 [1] CRAN (R 4.1.2)
#>  fs                     1.5.2      2021-12-08 [1] CRAN (R 4.1.2)
#>  future                 1.23.0     2021-10-31 [1] CRAN (R 4.1.2)
#>  gargle                 1.2.0      2021-07-02 [1] CRAN (R 4.1.2)
#>  gdtools                0.2.3      2021-01-06 [1] CRAN (R 4.1.2)
#>  generics               0.1.1      2021-10-25 [1] CRAN (R 4.1.2)
#>  geomtextpath           0.1.0      2022-01-14 [1] Github (AllanCameron/geomtextpath@6299079)
#>  gganimate              1.0.7      2020-10-15 [1] CRAN (R 4.1.2)
#>  ggdag                  0.2.4      2021-10-10 [1] CRAN (R 4.1.2)
#>  ggdist                 3.0.1      2021-11-30 [1] CRAN (R 4.1.2)
#>  ggforce                0.3.3      2021-03-05 [1] CRAN (R 4.1.2)
#>  gghighlight            0.3.2      2021-06-05 [1] CRAN (R 4.1.2)
#>  ggplot2              * 3.3.5      2021-06-25 [1] CRAN (R 4.1.2)
#>  ggraph                 2.0.5      2021-02-23 [1] CRAN (R 4.1.2)
#>  ggrepel                0.9.1      2021-01-15 [1] CRAN (R 4.1.2)
#>  ggridges               0.5.3      2021-01-08 [1] CRAN (R 4.1.2)
#>  gifski                 1.4.3-1    2021-05-02 [1] CRAN (R 4.1.2)
#>  globals                0.14.0     2020-11-22 [1] CRAN (R 4.1.2)
#>  glue                 * 1.6.0      2021-12-17 [1] CRAN (R 4.1.2)
#>  googledrive            2.0.0      2021-07-08 [1] CRAN (R 4.1.2)
#>  googlesheets4          1.0.0      2021-07-21 [1] CRAN (R 4.1.2)
#>  graphlayouts           0.8.0      2022-01-03 [1] CRAN (R 4.1.2)
#>  gridExtra              2.3        2017-09-09 [1] CRAN (R 4.1.2)
#>  gtable                 0.3.0      2019-03-25 [1] CRAN (R 4.1.2)
#>  gtools                 3.9.2      2021-06-06 [1] CRAN (R 4.1.2)
#>  haven                  2.4.3      2021-08-04 [1] CRAN (R 4.1.2)
#>  HDInterval             0.2.2      2020-05-23 [1] CRAN (R 4.1.2)
#>  here                 * 1.0.1      2020-12-13 [1] CRAN (R 4.1.2)
#>  highr                  0.9        2021-04-16 [1] CRAN (R 4.1.2)
#>  hms                    1.1.1      2021-09-26 [1] CRAN (R 4.1.2)
#>  hrbrthemes           * 0.8.0      2020-03-06 [1] CRAN (R 4.1.2)
#>  htmltools              0.5.2      2021-08-25 [1] CRAN (R 4.1.2)
#>  htmlwidgets            1.5.4      2021-09-08 [1] CRAN (R 4.1.2)
#>  httpuv                 1.6.5      2022-01-05 [1] CRAN (R 4.1.2)
#>  httr                   1.4.2      2020-07-20 [1] CRAN (R 4.1.2)
#>  ids                    1.0.1      2017-05-31 [1] CRAN (R 4.1.2)
#>  igraph                 1.2.11     2022-01-04 [1] CRAN (R 4.1.2)
#>  inline                 0.3.19     2021-05-31 [1] CRAN (R 4.1.2)
#>  isoband                0.2.5      2021-07-13 [1] CRAN (R 4.1.2)
#>  jquerylib              0.1.4      2021-04-26 [1] CRAN (R 4.1.2)
#>  jsonlite               1.7.2      2020-12-09 [1] CRAN (R 4.1.2)
#>  kableExtra             1.3.4      2021-02-20 [1] CRAN (R 4.1.2)
#>  knitr                  1.37       2021-12-16 [1] CRAN (R 4.1.2)
#>  labeling               0.4.2      2020-10-20 [1] CRAN (R 4.1.2)
#>  later                  1.3.0      2021-08-18 [1] CRAN (R 4.1.2)
#>  lattice                0.20-45    2021-09-22 [2] CRAN (R 4.1.2)
#>  lazyeval               0.2.2      2019-03-15 [1] CRAN (R 4.1.2)
#>  lifecycle              1.0.1      2021-09-24 [1] CRAN (R 4.1.2)
#>  listenv                0.8.0      2019-12-05 [1] CRAN (R 4.1.2)
#>  loo                  * 2.4.1      2020-12-09 [1] CRAN (R 4.1.2)
#>  lubridate              1.8.0      2021-10-07 [1] CRAN (R 4.1.2)
#>  magrittr               2.0.1      2020-11-17 [1] CRAN (R 4.1.2)
#>  markdown               1.1        2019-08-07 [1] CRAN (R 4.1.2)
#>  MASS                   7.3-54     2021-05-03 [2] CRAN (R 4.1.2)
#>  Matrix                 1.3-4      2021-06-01 [2] CRAN (R 4.1.2)
#>  matrixStats            0.61.0     2021-09-17 [1] CRAN (R 4.1.2)
#>  memoise                2.0.1      2021-11-26 [1] CRAN (R 4.1.2)
#>  mgcv                   1.8-38     2021-10-06 [2] CRAN (R 4.1.2)
#>  mime                   0.12       2021-09-28 [1] CRAN (R 4.1.2)
#>  miniUI                 0.1.1.1    2018-05-18 [1] CRAN (R 4.1.2)
#>  modelr                 0.1.8      2020-05-19 [1] CRAN (R 4.1.2)
#>  munsell                0.5.0      2018-06-12 [1] CRAN (R 4.1.2)
#>  mvtnorm                1.1-3      2021-10-08 [1] CRAN (R 4.1.2)
#>  nleqslv                3.3.2      2018-05-17 [1] CRAN (R 4.1.2)
#>  nlme                   3.1-153    2021-09-07 [2] CRAN (R 4.1.2)
#>  numDeriv               2016.8-1.1 2019-06-06 [1] CRAN (R 4.1.2)
#>  officedown             0.2.3      2021-11-16 [1] CRAN (R 4.1.2)
#>  officer                0.4.1      2021-11-14 [1] CRAN (R 4.1.2)
#>  openssl                1.4.6      2021-12-19 [1] CRAN (R 4.1.2)
#>  packrat                0.7.0      2021-08-20 [1] CRAN (R 4.1.2)
#>  parallelly             1.30.0     2021-12-17 [1] CRAN (R 4.1.2)
#>  pillar                 1.6.4      2021-10-18 [1] CRAN (R 4.1.2)
#>  pkgbuild               1.3.1      2021-12-20 [1] CRAN (R 4.1.2)
#>  pkgconfig              2.0.3      2019-09-22 [1] CRAN (R 4.1.2)
#>  plyr                   1.8.6      2020-03-03 [1] CRAN (R 4.1.2)
#>  polyclip               1.10-0     2019-03-14 [1] CRAN (R 4.1.2)
#>  posterior              1.2.0      2022-01-05 [1] CRAN (R 4.1.2)
#>  prettyunits            1.1.1      2020-01-24 [1] CRAN (R 4.1.2)
#>  processx               3.5.2      2021-04-30 [1] CRAN (R 4.1.2)
#>  progress               1.2.2      2019-05-16 [1] CRAN (R 4.1.2)
#>  promises               1.2.0.1    2021-02-11 [1] CRAN (R 4.1.2)
#>  ps                     1.6.0      2021-02-28 [1] CRAN (R 4.1.2)
#>  purrr                * 0.3.4      2020-04-17 [1] CRAN (R 4.1.2)
#>  R6                     2.5.1      2021-08-19 [1] CRAN (R 4.1.2)
#>  ragg                   1.2.1      2021-12-06 [1] CRAN (R 4.1.2)
#>  rappdirs               0.3.3      2021-01-31 [1] CRAN (R 4.1.2)
#>  ratlas               * 0.0.0.9000 2022-01-06 [1] Github (atlas-aai/ratlas@267dd5c)
#>  RColorBrewer           1.1-2      2014-12-07 [1] CRAN (R 4.1.2)
#>  Rcpp                 * 1.0.8      2022-01-13 [1] CRAN (R 4.1.2)
#>  RcppArmadillo          0.10.7.5.0 2021-12-17 [1] CRAN (R 4.1.2)
#>  RcppEigen              0.3.3.9.1  2020-12-17 [1] CRAN (R 4.1.2)
#>  RcppParallel           5.1.5      2022-01-05 [1] CRAN (R 4.1.2)
#>  readr                * 2.1.1      2021-11-30 [1] CRAN (R 4.1.2)
#>  readxl                 1.3.1      2019-03-13 [1] CRAN (R 4.1.2)
#>  rematch                1.0.1      2016-04-21 [1] CRAN (R 4.1.2)
#>  rematch2               2.1.2      2020-05-01 [1] CRAN (R 4.1.2)
#>  reprex                 2.0.1      2021-08-05 [1] CRAN (R 4.1.2)
#>  reshape2               1.4.4      2020-04-09 [1] CRAN (R 4.1.2)
#>  rethinking           * 2.21       2022-01-06 [1] Github (rmcelreath/rethinking@783d111)
#>  rlang                  0.4.12     2021-10-18 [1] CRAN (R 4.1.2)
#>  rmarkdown              2.11       2021-09-14 [1] CRAN (R 4.1.2)
#>  rprojroot              2.0.2      2020-11-15 [1] CRAN (R 4.1.2)
#>  rsconnect              0.8.25     2021-11-19 [1] CRAN (R 4.1.2)
#>  rstan                * 2.21.3     2021-12-19 [1] CRAN (R 4.1.2)
#>  rstantools             2.1.1      2020-07-06 [1] CRAN (R 4.1.2)
#>  rstudioapi             0.13       2020-11-12 [1] CRAN (R 4.1.2)
#>  Rttf2pt1               1.3.9      2021-07-22 [1] CRAN (R 4.1.2)
#>  rvest                  1.0.2      2021-10-16 [1] CRAN (R 4.1.2)
#>  rvg                    0.2.5      2020-06-30 [1] CRAN (R 4.1.2)
#>  sass                   0.4.0      2021-05-12 [1] CRAN (R 4.1.2)
#>  scales                 1.1.1      2020-05-11 [1] CRAN (R 4.1.2)
#>  selectr                0.4-2      2019-11-20 [1] CRAN (R 4.1.2)
#>  servr                  0.24       2021-11-16 [1] CRAN (R 4.1.2)
#>  sessioninfo            1.2.2      2021-12-06 [1] any (@1.2.2)
#>  shape                  1.4.6      2021-05-19 [1] CRAN (R 4.1.2)
#>  shiny                  1.7.1      2021-10-02 [1] CRAN (R 4.1.2)
#>  shinyjs                2.1.0      2021-12-23 [1] CRAN (R 4.1.2)
#>  shinystan              2.5.0      2018-05-01 [1] CRAN (R 4.1.2)
#>  shinythemes            1.2.0      2021-01-25 [1] CRAN (R 4.1.2)
#>  sourcetools            0.1.7      2018-04-25 [1] CRAN (R 4.1.2)
#>  StanHeaders          * 2.21.0-7   2020-12-17 [1] CRAN (R 4.1.2)
#>  stringi                1.7.6      2021-11-29 [1] CRAN (R 4.1.2)
#>  stringr              * 1.4.0      2019-02-10 [1] CRAN (R 4.1.2)
#>  svglite                2.0.0      2021-02-20 [1] CRAN (R 4.1.2)
#>  svUnit                 1.0.6      2021-04-19 [1] CRAN (R 4.1.2)
#>  sys                    3.4        2020-07-23 [1] CRAN (R 4.1.2)
#>  systemfonts            1.0.3      2021-10-13 [1] CRAN (R 4.1.2)
#>  tensorA                0.36.2     2020-11-19 [1] CRAN (R 4.1.2)
#>  textshaping            0.3.6      2021-10-13 [1] CRAN (R 4.1.2)
#>  threejs                0.3.3      2020-01-21 [1] CRAN (R 4.1.2)
#>  tibble               * 3.1.6      2021-11-07 [1] CRAN (R 4.1.2)
#>  tidybayes            * 3.0.2      2022-01-05 [1] CRAN (R 4.1.2)
#>  tidybayes.rethinking * 3.0.0      2022-01-06 [1] Github (mjskay/tidybayes.rethinking@7da9946)
#>  tidygraph              1.2.0      2020-05-12 [1] CRAN (R 4.1.2)
#>  tidyr                * 1.1.4      2021-09-27 [1] CRAN (R 4.1.2)
#>  tidyselect             1.1.1      2021-04-30 [1] CRAN (R 4.1.2)
#>  tidyverse            * 1.3.1      2021-04-15 [1] CRAN (R 4.1.2)
#>  tinytex                0.36       2021-12-19 [1] CRAN (R 4.1.2)
#>  tweenr                 1.0.2      2021-03-23 [1] CRAN (R 4.1.2)
#>  tzdb                   0.2.0      2021-10-27 [1] CRAN (R 4.1.2)
#>  utf8                   1.2.2      2021-07-24 [1] CRAN (R 4.1.2)
#>  uuid                   1.0-3      2021-11-01 [1] CRAN (R 4.1.2)
#>  V8                     4.0.0      2021-12-23 [1] CRAN (R 4.1.2)
#>  vctrs                  0.3.8      2021-04-29 [1] CRAN (R 4.1.2)
#>  viridis                0.6.2      2021-10-13 [1] CRAN (R 4.1.2)
#>  viridisLite            0.4.0      2021-04-13 [1] CRAN (R 4.1.2)
#>  vroom                  1.5.7      2021-11-30 [1] CRAN (R 4.1.2)
#>  webshot                0.5.2      2019-11-22 [1] CRAN (R 4.1.2)
#>  withr                  2.4.3      2021-11-30 [1] CRAN (R 4.1.2)
#>  xaringan               0.22       2021-06-23 [1] CRAN (R 4.1.2)
#>  xfun                   0.29       2021-12-14 [1] CRAN (R 4.1.2)
#>  xml2                   1.3.3      2021-11-30 [1] CRAN (R 4.1.2)
#>  xtable                 1.8-4      2019-04-21 [1] CRAN (R 4.1.2)
#>  xts                    0.12.1     2020-09-09 [1] CRAN (R 4.1.2)
#>  yaml                   2.2.1      2020-02-01 [1] CRAN (R 4.1.2)
#>  zip                    2.2.0      2021-05-31 [1] CRAN (R 4.1.2)
#>  zoo                    1.8-9      2021-03-09 [1] CRAN (R 4.1.2)
#> 
#>  [1] /home/runner/work/_temp/Library
#>  [2] /opt/R/4.1.2/lib/R/library
#> 
#> ──────────────────────────────────────────────────────────────────────────────
```
