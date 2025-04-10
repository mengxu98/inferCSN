# Reference: https://github.com/ThinkR-open/prepare-for-cran

# Prepare for CRAN ----

# Run tests and examples
devtools::document()
devtools::check()
devtools::test()
devtools::run_examples()

# Test pkgdown build
pkgdown::build_site()

# Check package as CRAN
rcmdcheck::rcmdcheck(args = c("--no-manual", "--as-cran"))
devtools::check_win_release()
devtools::check_win_devel()

# Update NEWS
# Bump version manually and add list of changes

# Verify you're ready for release, and release
devtools::release()
