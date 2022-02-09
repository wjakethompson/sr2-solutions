set.seed(2020)
options(digits = 3,
        knitr.kable.NA = "",
        brms.backend = "cmdstanr",
        brms.file_refit = "on_change")

suppressPackageStartupMessages(library(rethinking))
suppressPackageStartupMessages(library(brms))
suppressPackageStartupMessages(library(loo))
suppressPackageStartupMessages(library(ggridges))
suppressPackageStartupMessages(library(tidybayes))
suppressPackageStartupMessages(library(ggdist))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(wjake))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(glue))
suppressPackageStartupMessages(library(geomtextpath))
suppressPackageStartupMessages(library(ggtext))
suppressPackageStartupMessages(library(dagitty))

ggplot2::theme_set(wjake::theme_wjake(base_family = "Source Sans Pro"))

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

### Misc Functions -------------------------------------------------------------
seq_range <- function (x, n) {
  seq(min(x), max(x), length.out = n)
}
