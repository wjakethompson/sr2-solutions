set.seed(2020)
options(digits = 3)

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(rethinking))
suppressPackageStartupMessages(library(brms))
suppressPackageStartupMessages(library(hrbrthemes))
suppressPackageStartupMessages(library(here))

ggplot2::theme_set(theme_minimal())

knitr::opts_chunk$set(
  comment = "#>",
  collapse = TRUE,
  cache = TRUE,
  out.width = "70%",
  fig.align = "center",
  fig.width = 6,
  fig.asp = 0.618,
  fig.show = "hold"
)
