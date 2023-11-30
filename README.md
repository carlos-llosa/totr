# totr

Tensor-on-Tensor Regression with Kronecker Separable Covariance 

This is an R package to do linear regression with tensor responses and tensor covariates [1]. This implementation allows modeling the voxel-to-voxel correlation in the response tensor, and allows for four different types of low-rank regression coefficients.

The `totr` package can be directly installed from GitHub, using the devtools library:

```
install.packages("devtools")
library(devtools)
install_github("carlos-llosa/totr")
library(totr)
``` 

[1] Llosa-Vite, C., & Maitra, R. (2022). Reduced-Rank Tensor-on-Tensor Regression and Tensor-variate Analysis of Variance. IEEE Transactions on Pattern Analysis and Machine Intelligence, 45(2), 2282 - 2296.[(link to article)](https://doi.org/10.1109/TPAMI.2022.3164836)
