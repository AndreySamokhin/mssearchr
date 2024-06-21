# mssearchr

The primary goal of the `mssearchr` package is to enhance the capabilities of R
users for conducting library searches against electron ionization mass spectral
databases.

The `mssearchr` package offers the following tools:
- custom implementation of the Identity (EI Normal) algorithm;
- custom implementation of the Similarity (EI Simple) algorithm;
- library search via NIST API (i.e., by calling the *nistms$.exe* file);
- reading/writing *msp* files.



## Installation

``` r
# Install 'mssearchr' from CRAN:
install.packages("mssearchr")

# Install 'mssearchr' from GitHub:
library(devtools)
install_github("https://github.com/AndreySamokhin/mssearchr")
```

