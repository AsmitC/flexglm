# flexglm

`flexglm` fits various Bayesian semiparametric generalized linear models in R,
all of which model the baseline distribution nonparametrically, allowing for
more flexible inference. The regression parameters are unchanged.

The package currently supports three models:

- `dirglm()` places a Dirichlet prior on a finite response support.
  See [Alam et al. (2025)](https://doi.org/10.1002/sim.10305)
  for details.
- `dpglm()` places a Dirichlet-process prior on the baseline distribution.
  See [Alam et al. (2025)](https://doi.org/10.48550/arXiv.2502.17827)
  for details.
- `cdpglm()` extends the Dirichlet-process model to clustered or longitudinal
  observations by using a copula. See Yuan et al. (manuscript in preparation)
  for details. Currently, exchangeable and AR(1) correlation structures
  are supported.

All models are fit using Markov chain Monte Carlo. The returned objects
provide posterior draws, acceptance rates, summary methods, and
plotting methods for model parameters.

## Development status

`flexglm` is under active development and has not yet been released.
The statistical methodology and software interface are subject to
change before stable releases are made available.

## Installation

The development version can be installed from GitHub with:

```r
install.packages("remotes")
remotes::install_github("AsmitC/flexglm")
```

Load the package with:

```r
library(flexglm)
```

Installation from GitHub builds the package from source. A working C++
toolchain is therefore required. Some examples include: Rtools on Windows,
Xcode Command Line Tools on macOS, or typical compilers and R development packages on Linux.

## Basic use

Detailed example workflows for fitting each model are contained
in the `vignettes` folder, any of which can be accessed as follows:

```r
vignette("dirglm_example", package = "flexglm")
vignette("dpglm_example", package = "flexglm")
vignette("cdpglm_example", package = "flexglm")
```

## Documentation

Function-level documentation is available through the usual R help system:

```r
help("dirglm", package = "flexglm") # or ?dirglm if the package has been loaded
help("dpglm", package = "flexglm") # or ?dpglm if the package has been loaded
help("cdpglm", package = "flexglm") # or ?cdpglm if the package has been loaded
```

The corresponding control functions are `dirglm.control()`, `dpglm.control()`,
and `cdpglm.control()`. Each function provides extensive customization options
for users, though most users may find the default settings satisfactory.

## Reporting problems

Please report reproducible software defects through the repository's
[issue tracker](https://github.com/AsmitC/flexglm/issues). Include the code and
data needed to reproduce the problem, the complete error or warning message,
and the output of `sessionInfo()`.

## License

`flexglm` is licensed under the [MIT License](LICENSE.md).
