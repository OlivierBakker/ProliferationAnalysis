# ProliferationAnalysis

# Installation

Depends on
- devtools
- minpack.lm https://cran.r-project.org/web/packages/minpack.lm/index.html

```
devtools::install_github("https://github.com/OlivierBakker/ProliferationAnalysis/tree/main")
```

# Example

The following simulates some proliferation data and fits peaks to it. 

Or see the vingette for more detaills: <a href="https://html-preview.github.io/?url=https://github.com/OlivierBakker/ProliferationAnalysis/blob/main/vignettes/fitting_proliferation_model.html" >Fitting a proliferation model with complex mixture</a>
Note: download the vingette html in the vingettes folder for it to render properly

```
library(ProliferationAnalysis)

# Simulate proliferation data
y <- 10 ^ rnorm(1000, mean=10, sd=0.05)
y <- c(y/4, y/2, y)

# Fit peaks using least squares
peaks <- fit.peaks(y, peak.0.lower.bound=9.8)
get.prolif.stats(peaks)


# Simulate proliferation data #2
y <- 10 ^ rnorm(1000, mean=10, sd=0.05)
y <- c((y/8)[1:500], y/4, y/2, y[1:500])

# Fit peaks using MLE
peaks <- fit.peaks(y, peak.0.lower.bound=9.8, mode="MLE")

get.prolif.stats(peaks)
```
