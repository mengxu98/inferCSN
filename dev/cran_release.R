# Ref: https://github.com/ThinkR-open/prepare-for-cran

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

# Update NEWS
# Bump version manually and add list of changes

# Add comments for CRAN
usethis::use_cran_comments(open = rlang::is_interactive())

# Verify you're ready for release, and release
devtools::release()
