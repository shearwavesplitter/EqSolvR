### EqSolvR
 Package for solving chemical equilibria for a given set of reactants and products. The motivation for this program was to enable the calculation, between 300°C and 400°C at 0.5 kb, of pH and speciation given a simple mix of salts. This package has been writting in such a manner that an advanced user can easily set their own reactants, products and temperature in the generic version (chemsolve_generic). chemsolve_generic is the workhorse for chemsolve. Equations are solved numerically using the multiroot function from rootSolve package.
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
## Defaults to a charge balance against Cl (bal="Cl")
## Add H2SO4 and MgCl2 as an additional products (the K values, exK, here are just examples)
d <- chemsolve(exprod=c("H2SO4","MgCl2"),exconstit=c("H","H","SO4","Mg","Cl","Cl"),exnumz=c(3,3),excharges=c(0,0),exa=c(0,0),exK=c(-6,-3))
```