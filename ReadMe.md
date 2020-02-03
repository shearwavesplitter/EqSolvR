### EqSolvR
 Package for solving chemical equilibria for a given set of reactants and products. The motivation for this program was to enable the calculation, between 300°C and 400°C at 0.5 kb, of pH and speciation given a simple mix of salts. This package has been written in such a manner that an advanced user can easily set their own reactants, products and log K at any temperature in the generic version (chemsolve_generic). chemsolve_generic is the workhorse for chemsolve. Equations are solved numerically using the multiroot function from rootSolve package.
 
The **[eq_creator.xlsx](https://github.com/shearwavesplitter/EqSolvR/blob/master/eq_creator.xlsx)** excel spreadsheet, created by Lucjan Sajkowski, can be used easily format the R inputs.

#### [Manual](https://github.com/shearwavesplitter/EqSolvR/blob/master/EqSolvR.pdf)

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
d <- chemsolve(Tc = 400, Nat = 0.4, Kt = 0.2, Clt = 0.6, SO4t = 0.2,Cat = 0.1, Mgt = 0.1, start = c(1e-06, 1e-05, 0.3, 0.1, 0.3,0.01, 0.001, 0.02), maxitr = 100, exprod = NULL, exconstit = NULL,exnumz = NULL, excharges = NULL, exa = NULL, exK = NULL, bal = "Cl")

```
#### To cite this package
Stefan Mroczek and Ed Mroczek (2017). EqSolvR: Chemical Equilibrium Solver. R package version 1.2.5.

```latex
  @Manual{,
    title = {EqSolvR: Chemical Equilibrium Solver},
    author = {Stefan Mroczek and Ed Mroczek},
    year = {2017},
    note = {R package version 1.2.5},
  }
```
