# Development dependencies
# Run from the project root after activating renv.
renv::install()
renv::install(c(
  "devtools", "roxygen2", "lintr", "styler", "languageserver", "jsonlite"
))
