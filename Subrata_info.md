
Sun Feb  5 09:26:21 PM CST 2023
===============================

tensorA package is used for some purpose in the backend. 
Somehow this package was creating problem while installing:
```{bash}
* installing to library ‘/home/subrata/R/x86_64-pc-linux-gnu-library/4.2’
ERROR: dependency ‘tensorA’ is not available for package ‘totr’
* removing ‘/home/subrata/R/x86_64-pc-linux-gnu-library/4.2/totr’
```
However, it can be installed using devtools for now.


## About tensorA package:
The functions which would be possibly needed for tensor:
```{R}
to.tensor
level.tensor
ftable

%e%
$r%
mul.tensor
einstein.tensor (?)
diagmul.tensor

diag.tensor
delta.tensor
bind.tensor
slice.tensor
tripledelta.tensor
undrop.tensor
to.matrix.tensor
to.tensor
toPos.tensor
trace.tensor

mean.tensor
var.tensor
chol.tensor
svd.tensor
inv.tensor
solve.tensor
```


```{R}
A <- to.tensor( 1:20, c(a=2,b=2,c=5) )
A
ftable(A)
B <- to.tensor( c(0,1,1,0) , c(a=2,"a'"=2))
A %e% B
drag.tensor( A , B, c("a","b"))
A %e% one.tensor(c(c=5))/5 # a mean of matrices
reorder.tensor(A,c("c","b","a"))
A - reorder.tensor(A,c("c","b","a")) # =0 since sequence is irrelevant
inv.tensor(A,"a",by="c")
```