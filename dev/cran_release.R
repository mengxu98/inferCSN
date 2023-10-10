# Ref: https://github.com/ThinkR-open/prepare-for-cran

# Prepare for CRAN ----

# Run tests and examples
devtools::document()
devtools::test()
devtools::run_examples()

# Check package as CRAN
rcmdcheck::rcmdcheck(args = c("--no-manual", "--as-cran"))

# Check content
if (!require("checkhelper", quietly = TRUE)) {
  install.packages("checkhelper", repos = 'https://thinkr-open.r-universe.dev')
}
checkhelper::find_missing_tags()
# # _Check that you let the house clean after the check, examples and tests
# all_files_remaining <- checkhelper::check_clean_userspace()
# all_files_remaining

# Check spelling
if (!require("spelling", quietly = TRUE)) {
  install.packages("spelling")
}
spelling::spell_check_test(vignettes = TRUE,
                           error = FALSE,
                           skip_on_cran = TRUE)
spelling::spell_check_package()

# check on other distributions
# _win devel
devtools::check_win_devel()
# _rhub
# devtools::check_rhub()
rhub::check_on_windows(check_args = "--force-multiarch")
rhub::check_on_debian(check_args = "--force-multiarch")
# rhub::check_on_solaris()

# Check reverse dependencies
usethis::use_git_ignore("revdep/")
usethis::use_build_ignore("revdep/")

devtools::revdep()

# Update NEWS
# Bump version manually and add list of changes

# Add comments for CRAN
usethis::use_cran_comments(open = rlang::is_interactive())

# Upgrade version number
usethis::use_version(which = c("patch", "minor", "major", "dev")[1])

# Verify you're ready for release, and release
devtools::release()
