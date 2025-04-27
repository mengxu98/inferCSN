# -*- coding: utf-8 -*-

#' @title inferCSN: inferring cell-type specific gene regulatory network
#'
#' @description
#' An R package for inferring cell-type specific gene regulatory network from single-cell RNA-seq data
#'
#' @author Meng xu (Maintainer), \email{mengxu98@qq.com}
#'
#' @source \url{https://github.com/mengxu98/inferCSN}
#'
#' @md
#' @docType package
#' @name inferCSN-package
"_PACKAGE"

#' @title inferCSN logo
#'
#' @description
#' The inferCSN logo, using ASCII or Unicode characters
#' Use [cli::ansi_strip()] to get rid of the colors.
#' @param unicode Unicode symbols. Default is `TRUE` on UTF-8 platforms.
#'
#' @references
#'  \url{https://github.com/tidyverse/tidyverse/blob/main/R/logo.R}
#'
#' @md
#' @export
#' @examples
#' inferCSN_logo()
inferCSN_logo <- function(
    unicode = cli::is_utf8_output()) {
  logo <- c(
    "           2   1                3
            ____     4    ___________ _   __
    0 ___  / __/__  _____/ ____/ ___// | / /
  / / __ ./ /_/ _ ./ ___/ /    |__ ./  |/ /
 / / / / / __/  __/ /  / /___ ___/ / /|  /
/_/_/ /_/_/  .___/_/   .____//____/_/ |_/
    6      5            7      8     9   "
  )

  hexa <- c("*", ".", "o", "*", ".", "*", ".", "o", ".", "*")
  if (unicode) {
    hexa <- c("*" = "\u2b22", "o" = "\u2b21", "." = ".")[hexa]
  }

  cols <- c(
    "red", "yellow", "green", "magenta", "cyan",
    "yellow", "green", "white", "magenta", "cyan"
  )

  col_hexa <- purrr::map2(
    hexa, cols, ~ cli::make_ansi_style(.y)(.x)
  )

  for (i in 0:9) {
    pat <- paste0("\\b", i, "\\b")
    logo <- sub(pat, col_hexa[[i + 1]], logo)
  }

  structure(cli::col_blue(logo), class = "logo")
}

#' @title print logo
#'
#' @param x Input infromation.
#' @param ... Other parameters.
#'
#' @method print logo
#'
#' @export
print.logo <- function(x, ...) {
  cat(x, ..., sep = "\n")
  invisible(x)
}

.onAttach <- function(libname, pkgname) {
  version <- utils::packageDescription(pkgname, fields = "Version")

  msg <- paste0(
    "---------------------------------------------------
",
    cli::col_blue(pkgname, " version ", version),
    "
This message can be suppressed by:
  suppressPackageStartupMessages(library(inferCSN))
---------------------------------------------------"
  )

  packageStartupMessage(inferCSN_logo())
  packageStartupMessage(msg)
}
