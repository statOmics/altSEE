# `iSEEu`

<!-- badges: start -->
<!-- badges: end -->

The `altSEE` package contains material and code that extends the `iSEE` package (https://github.com/iSEE/iSEE).

We welcome contributions from the community, see below for more instructions. 
For example, during the Developer Day at the European Bioconductor 2019 conference (`#EuroBioc2019`, at the UCLouvain, in Brussels, Belgium), we proposed a hackathon-like session, and we focused on the design of "modes", i.e. preconfigured sets of panels and linked content to be used as starting setup when launching `iSEE`.

## Installation

`altSEE` can be easily installed from GitHub, you can use:

```r
BiocManager::install("statomics/altSEE", dependencies = TRUE)
# or alternatively...
remotes::install_github("statomics/altSEE", dependencies = TRUE)
```

Setting `dependencies = TRUE` should ensure that all packages, including the ones in the `Suggests:` field of the `DESCRIPTION`, are installed - this can be essential if you want to reproduce the code in the vignette, for example.


## Expanding the `iSEE` universe with `altSEE`

- install `iSEE` first - the development version is recommended.

```r
BiocManager::install("iSEE", version = "devel")
# or
remotes::install_github("iSEE/iSEE")
```

- fork the `altSEE` repo (https://github.com/statomics/altSEEu) and clone it locally.

```bash
git clone https://github.com/[your_github_username]/altSEE.git
```

- make the desired changes in the files - start from the `R` folder, then document via `roxygen2` - and push to your fork. 

- once your contribution (function, panel, mode) is done, consider adding some information in the package.
  Some examples might be a screenshot of the mode in action (to be placed in the folder `inst/modes_img`), or a well-documented example use case (maybe an entry in the `vignettes` folder). Also add yourself as a contributor (`ctb`) to the DESCRIPTION file.

- make a pull request to the original repo - the GitHub site offers a practical framework to do so, enabling comments, code reviews, and other goodies.

- more on documenting and code guidelines:
  - if possible, please consider adding an example in the dedicated Roxygen preamble to show how to run each function
  - if possible, consider adding one or more unit tests - we use the `testthat` framework
