# VS Code include paths for the active R environment
# Run from the project root after renv::restore().
local({
  packages <- c("Rcpp", "RcppParallel")
  includes <- vapply(packages, function(package) {
    path <- system.file("include", package = package)
    if (!nzchar(path)) stop("Install ", package, " before setting up VS Code.")
    normalizePath(path, winslash = "/")
  }, character(1))

  compiler <- unname(Sys.which("clang++"))
  if (!nzchar(compiler)) compiler <- unname(Sys.which("g++"))
  if (!nzchar(compiler)) {
    stop("Install a C++ compiler before setting up VS Code.")
  }

  configuration <- list(
    name = "R package",
    includePath = c(
      "${workspaceFolder}/src",
      normalizePath(R.home("include"), winslash = "/"),
      unname(includes)
    ),
    defines = I("NDEBUG"),
    compilerPath = compiler,
    cStandard = "c17",
    cppStandard = "c++20"
  )

  dir.create(".vscode", showWarnings = FALSE)
  jsonlite::write_json(
    list(configurations = list(configuration), version = 4L),
    ".vscode/c_cpp_properties.json", auto_unbox = TRUE, pretty = TRUE
  )
  settings_path <- ".vscode/settings.json"
  settings <- if (file.exists(settings_path)) {
    jsonlite::read_json(settings_path)
  } else {
    list()
  }
  settings[["r.useRenvLibPath"]] <- TRUE
  jsonlite::write_json(
    settings, settings_path, auto_unbox = TRUE, pretty = TRUE
  )
  message("Updated .vscode/c_cpp_properties.json for the active R library.")
})
