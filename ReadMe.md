### EqSolvR
Package for solving chemical equilibria for a given set of reactants and products
#### Installation

Install R by following the instructions on [https://cran.r-project.org](https://cran.r-project.org/) 

Open R and install the devtools package

```r
install.packages("devtools")
```

And finally install EqSolvR

```r
devtools::install_github("shearwavesplitter/EqSolvR")
```

To use EqSolvR

```r
library('EqSolvR')
```

#### Example
```r
library('EqSolvR')
## Run the solver with the default species of Na, Cl, K, SO4, Ca, and Mg 
## Add H2SO4 as an additional product
d <- chemsolve(Tc=300,Nat=0.2,Kt=0.2,Clt=0.4,SO4t=0.2,Cat=0.1,Mgt=0.1,exprod="H2SO4",exconstit=c("H","H","SO4"),exnumz=3,excharges=0,exa=0,exK=-6)
```