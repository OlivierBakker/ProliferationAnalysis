# ProliferationAnalysis

# Installation
To install, simply run
```
devtools::install_github("https://github.com/OlivierBakker/ProliferationAnalysis/tree/main")
```
Code is all in base R except for the Levanberg-Marquad algorithm which uses
the implementation in minpack.lm

# Example
The following simulates some proliferation data and fits peaks to it. 

Or see the vingette for more detaills: <a href="https://html-preview.github.io/?url=https://github.com/OlivierBakker/ProliferationAnalysis/blob/main/vignettes/fitting_proliferation_model.html" >Fitting a proliferation model with complex mixture</a>

```
library(ProliferationAnalysis)

# Simulate proliferation data
y <- 10 ^ rnorm(1000, mean=10, sd=0.05)
y <- c(y/4, y/2, y)

# Fit peaks using least squares
peaks <- fit.peaks(y, peak.0.lower.bound=9.8)
get.prolif.stats(peaks)

# Fit peaks using MLE
peaks <- fit.peaks(y, peak.0.lower.bound=9.8, mode="MLE")

get.prolif.stats(peaks)
```

