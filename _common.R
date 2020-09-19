set.seed(2020)
options(digits = 3,
        knitr.kable.NA = "")

suppressPackageStartupMessages(library(rethinking))
suppressPackageStartupMessages(library(brms))
suppressPackageStartupMessages(library(loo))
suppressPackageStartupMessages(library(tidybayes))
suppressPackageStartupMessages(library(tidybayes.rethinking))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(hrbrthemes))
suppressPackageStartupMessages(library(ratlas))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(glue))

ggplot2::theme_set(theme_minimal())

knitr::opts_chunk$set(
  comment = "#>",
  collapse = TRUE,
  cache = FALSE,
  dev = "ragg_png",
  dpi = 96L,
  fig.ext = "png",
  fig.width = 700 / 96,
  fig.retina = 3L,
  fig.asp = 1 / 1.618,
  fig.show = "hold",
  fig.align = "center",
  out.width = "70%"
)

knitr::knit_hooks$set(wrap = function(before, options, envir){
  if (before){
    paste0('<div class="', options$wrap, '">')
  } else {
    paste0('</div>')
  }
})
