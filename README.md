# zebu

## Description

The `zebu` R package implements different tools related to local association measures in three main functions. It replaces the previous web application written using `shiny`.

- `lassie` estimates local (and global) association measures: Ducher's Z, pointwise mutual information, normalized pointwise mutual information and chi-squared residuals.

- `permtest` accesses the significance of local (and global) association values using p-values estimated by permutations.

- `subgroups` identifies if the association between variables is dependent on the value of another variable.

## Installation

Install [R](https://r-project.org/) and [RStudio](https://rstudio.com/).

Get the released version from CRAN:
```R
install.packages("zebu")
```

or the development version from Github:

```R
# install.packages("devtools")
devtools::install_github("oliviermfmartin/zebu")
```

## Usage and Theory

Please read the vignette.

```R
vignette("zebu")
```

## Contact

* Olivier Martin
[oliviermfmartin@tutanota.com](mailto:oliviermfmartin@tutanota.com)

* Michel Ducher
[michel.ducher@chu-lyon.fr](mailto:michel.ducher@chu-lyon.fr)

## License

zebu 0.1
Copyright (C) 2017

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

The GNU General Public License is available at https://gnu.org/licenses/

The source code can be found at https://github.com/oliviermfmartin/zebu
