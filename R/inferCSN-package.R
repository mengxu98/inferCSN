# -*- coding: utf-8 -*-

#' @title Inferring Cell-Specific Gene Regulatory Network
#'
#' @useDynLib inferCSN
#'
#' @description
#' An R package for inferring cell-type specific gene regulatory network from single-cell RNA-seq data.
#'
#' @author Meng Xu (Maintainer), \email{mengxu98@qq.com}
#'
#' @source \url{https://mengxu98.github.io/inferCSN/}
#'
#' @md
#' @docType package
#' @name inferCSN-package
"_PACKAGE"

#' @title The logo of inferCSN
#'
#' @description
#' The inferCSN logo, using ASCII or Unicode characters
#' Use [cli::ansi_strip] to get rid of the colors.
#'
#' @md
#' @param unicode Unicode symbols on UTF-8 platforms.
#' Default is [cli::is_utf8_output].
#'
#' @references
#' \url{https://github.com/tidyverse/tidyverse/blob/main/R/logo.R}
#'
#' @export
#' @examples
#' infercsn_logo()
infercsn_logo <- function(unicode = cli::is_utf8_output()) {
  logo <- c(
    "          0          1        2             3     4
           _        ____           ___________ _   __
          (_)____  / __/___  _____/ ____/ ___// | / /
         / // __ ./ /_ / _ ./ ___/ /    .__ ./  |/ /
        / // / / / __//  __/ /  / /___ ___/ / /|  /
       /_//_/ /_/_/   .___/_/   .____//____/_/ |_/
      5               6      7        8          9"
  )

  hexa <- c("*", ".", "o", "*", ".", "o", "*", ".", "o", "*")
  if (unicode) {
    hexa <- c("*" = "\u2b22", "o" = "\u2b21", "." = ".")[hexa]
  }

  cols <- c(
    "red", "yellow", "green", "magenta", "cyan",
    "yellow", "green", "white", "magenta", "cyan"
  )

  col_hexa <- mapply(
    function(x, y) cli::make_ansi_style(y)(x),
    hexa, cols,
    SIMPLIFY = FALSE
  )

  for (i in 0:9) {
    pat <- paste0("\\b", i, "\\b")
    logo <- sub(pat, col_hexa[[i + 1]], logo)
  }

  structure(
    cli::col_blue(logo),
    class = "infercsn_logo"
  )
}

#' @title Print logo
#'
#' @param x Input information.
#' @param ... Other parameters.
#'
#' @return Print the ASCII logo
#'
#' @method print infercsn_logo
#'
#' @export
#'
print.infercsn_logo <- function(x, ...) {
  cat(x, ..., sep = "\n")
  invisible(x)
}

.onAttach <- function(libname, pkgname) {
  verbose <- thisutils::get_verbose()
  if (isTRUE(verbose)) {
    version <- utils::packageDescription(
      pkgname,
      fields = "Version"
    )

    msg <- paste0(
      strrep("-", 60),
      "\n",
      cli::col_blue(pkgname, " version ", version),
      "\n",
      cli::col_grey("This message can be suppressed by:"),
      "\n",
      cli::col_grey("  suppressPackageStartupMessages(library(inferCSN))"),
      "\n",
      cli::col_grey("  or options(log_message.verbose = FALSE)"),
      "\n",
      strrep("-", 60)
    )

    packageStartupMessage(infercsn_logo())
    packageStartupMessage(msg)
  }
}
