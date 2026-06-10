# SABERPRO — Semi-Analytical Bayesian Estimate Retrieval

A semi-analytically parameterised aquatic radiative transfer model for retrieving posterior distributions of Optically Significant Constituents (OSCs), water depth, and benthic reflectance from remote sensing reflectance (Rrs).

**Creator and lead developer:** Soham Mukherjee  
**Initial packaging:** Raphael Mabit  
**Maintenance and feature integration:** Soham Mukherjee

---

## Mathematical and physics-based formulation

Please refer to the publication:

> *A Semi-Analytical Bayesian Estimate Retrieval (SABER) algorithm for the inversion of Remote Sensing Reflectance in optically deep and shallow waters*  
> <https://doi.org/10.1002/lom3.70004>

---

## Code structure

This package follows the recommendations of <https://r-pkgs.org/> and the [tidyverse style guide](https://style.tidyverse.org/).

The code is written with a functional approach:

- **`fct_*` files** — low-level computations that call compiled C functions in `src/`
- **`utils_*` files** — higher-level, more generic utilities

The central piece is the `objective_factory()` function, which stitches together any forward model with any objective/likelihood function (including user-defined ones).

Forward models, input preparers, and objective functions are stored in registries (see `registry.R`): `.input_preparer_registry`, `.forward_model_registry`, `.objective_function_registry`.

---

## Installation

### System requirements

| Platform | Requirements |
|----------|-------------|
| **R** | ≥ 3.5.0 |
| **Windows** | [Rtools](https://cran.r-project.org/bin/windows/Rtools/) matching your R version |
| **Linux (Ubuntu/Debian)** | `build-essential`, `gfortran`, `cmake`, `libnlopt-dev`, `pkg-config` |
| **macOS** | Xcode Command Line Tools, `pkg-config` (via Homebrew) |

On Ubuntu/Debian, install system dependencies in one step:

```bash
sudo apt install build-essential gfortran cmake libnlopt-dev pkg-config \
                 libssl-dev libcurl4-openssl-dev libxml2-dev
```

### Install from GitHub

```r
# Install devtools if needed
if (!require("devtools")) install.packages("devtools")

# Install SABERPRO
devtools::install_github("homas01123/SABERPRO", dependencies = TRUE)
```

### Verify installation

```r
library(SABERPRO)

# List registered forward models and objective functions
list_forward_models()
list_objective_functions()
```

---

## Running the code

See the [forward_inverse_basics](vignettes/forward_inverse_basics.Rmd) vignette for a worked introduction to forward and inverse modelling with SABERPRO.

---

## Troubleshooting

**Windows:** Ensure Rtools is installed and on the system PATH.

**Linux:** If compilation fails with errors about missing math functions (`cos`, `exp`, etc.), ensure `gfortran` and `libnlopt-dev` are installed (see system requirements above).

**General:** If you see `Error: no benthic classes loaded`, call `select_benthic_classes()` before `make_inversion_params()`.
