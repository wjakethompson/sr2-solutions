set.seed(2020)
options(digits = 3,
        knitr.kable.NA = "",
        brms.backend = "cmdstanr")

suppressPackageStartupMessages(library(rethinking))
suppressPackageStartupMessages(library(brms))
suppressPackageStartupMessages(library(loo))
suppressPackageStartupMessages(library(tidybayes))
suppressPackageStartupMessages(library(ggdist))
suppressPackageStartupMessages(library(tidybayes.rethinking))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(hrbrthemes))
suppressPackageStartupMessages(library(ratlas))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(glue))
suppressPackageStartupMessages(library(geomtextpath))
suppressPackageStartupMessages(library(dagitty))

ggplot2::theme_set(ggplot2::theme_minimal(base_family = "Source Sans Pro"))

ggplot2::update_geom_defaults("text", list(family = "Source Sans Pro"))
ggplot2::update_geom_defaults("label", list(family = "Source Sans Pro"))
ggplot2::update_geom_defaults("textdensity", list(family = "Source Sans Pro"))

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
  out.width = "80%"
)

knitr::knit_hooks$set(wrap = function(before, options, envir){
  if (before){
    paste0('<div class="', options$wrap, '">')
  } else {
    paste0('</div>')
  }
})
