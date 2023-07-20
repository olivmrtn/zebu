# zebu: Local Association Measures

## Description

The `zebu` R package implements different tools related to local association measures in three main functions. It replaces the previous web application written using `shiny`.

- `lassie` estimates local (and global) association measures: Ducher's Z, pointwise mutual information, normalized pointwise mutual information and chi-squared residuals.

- `permtest` accesses the significance of local (and global) association values using p-values estimated by permutations.

- `chisqtest` accesses the significance for two dimensional chi-squared analysis.

## Installation

Install [R](https://www.r-project.org/) and [RStudio](https://posit.co/download/rstudio-desktop/).

Get the released version from [CRAN](https://CRAN.R-project.org/package=zebu).

```R
install.packages("zebu")
```

## Usage and Theory

Please read the [vignette](https://CRAN.R-project.org/package=zebu/vignettes/zebu.html).

```R
vignette("zebu")
```

## Contact

* Olivier Martin
[oliviermfmartin@tutanota.com](mailto:oliviermfmartin@tutanota.com)

* Michel Ducher
[michel.ducher@chu-lyon.fr](mailto:michel.ducher@chu-lyon.fr)

## License

zebu 0.2.1
Copyright (C) 2023

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

The GNU General Public License is available at https://gnu.org/licenses/

The source code can be found at https://github.com/oliviermfmartin/zebu
