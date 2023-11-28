# ProliferationAnalysis

# Installation

Depends on
- devtools
- minpack.lm https://cran.r-project.org/web/packages/minpack.lm/index.html

```
install.packages(minpack.lm)
devtools::install_github("https://github.com/OlivierBakker/ProliferationAnalysis/tree/main")
```

# Example

The following simulates some proliferation data and fits peaks to it
```
# Simulate proliferation data
y <- 10 ^ rnorm(1000, mean=10, sd=0.05)
y <- c(y/4, y/2, y)

# Fit peaks
peaks <- fit.peaks(y, 0, 0, peak.0.lower.bound=9.8)
get.prolif.stats(peaks)
```
