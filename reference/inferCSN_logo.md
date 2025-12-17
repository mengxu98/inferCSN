# The logo of inferCSN

The inferCSN logo, using ASCII or Unicode characters Use
[cli::ansi_strip](https://cli.r-lib.org/reference/ansi_strip.html) to
get rid of the colors.

## Usage

``` r
infercsn_logo(unicode = cli::is_utf8_output())
```

## Arguments

- unicode:

  Unicode symbols on UTF-8 platforms. Default is
  [cli::is_utf8_output](https://cli.r-lib.org/reference/is_utf8_output.html).

## References

<https://github.com/tidyverse/tidyverse/blob/main/R/logo.R>

## Examples

``` r
infercsn_logo()
#>           ⬢          .        ⬡             ⬢     .
#>            _        ____           ___________ _   __
#>           (_)____  / __/___  _____/ ____/ ___// | / /
#>          / // __ ./ /_ / _ ./ ___/ /    .__ ./  |/ /
#>         / // / / / __//  __/ /  / /___ ___/ / /|  /
#>        /_//_/ /_/_/   .___/_/   .____//____/_/ |_/
#>       ⬡               ⬢      .        ⬡          ⬢
```
