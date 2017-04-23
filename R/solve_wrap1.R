#' @title Mass balance and charge solver
#' @description Mass balance and charge balance solver for chemical equilibria. 
#' @param Tc Temperature (degrees C)
#' @param Nat Sodium concentration (mol/kg)
#' @param Kt Potassium (mol/kg)
#' @param Clt Chlorine (mol/kg)
#' @param SO4t Sulphate (mol/kg)
#' @param Cat Calcium (mol/kg)
#' @param Mgt Magnesium (mol/kg)
#' @param start Initial guess for the calculated equalibrium concentration of the basis species
#' @param maxitr Maximum number of iteration
#' @param exprod A vector of the names of the additional product species which react to form the basis species
#' @param exconstit  A vector of the chemical symbols of the additional product species in terms of the basis species
#' @param exnumz A vector of the number of basis consituents for each of the additional product species given by the equilibrium equation
#' @param excharges A vector of the total charage of the additional product species
#' @param exK A vector of the log K values of the additional product species
#' @param exa A vector of the ion size paramters for the additional species
#' @details A wrapper for the chemsolve_generic function that allow easy addition of product species. If you want to add additional reactant species (i.e. basis species) then the chemsolve_generic function must be used
#' @return A list containing the concentrations, gamma values, and pH at equilibrium
#' @export
#' @examples
#' ## Add H2SO4 as an extra product species
#'
#' chemsolve(exprod="H2SO4",exconstit="H","H","SO4",exnumz=3,excharges=0,exa=0,exa=-6)
#' 
chemsolve <- function(Tc=300,Nat=0.2,Kt=0.2,Clt=0.4,SO4t=0.2,Cat=0.1,Mgt=0.1,start=c(0.00001,0.00001,0.15,0.15,0.15,0.104756881,0.05,0.05),maxitr=100,exprod=NULL,exconstit=NULL,exnumz=NULL,excharges=NULL,exa=NULL,exK=NULL) {

spec <- c("Na","K","Cl","SO4","Ca","Mg")
concz <- c(Nat,Kt,Clt,SO4t,Cat,Mgt)
speccharges <- c(1,1,-1,-2,2,2)
speca <- c(4,3,3.5,4,6,8)


products <- c("NaSO4","HSO4","KSO4","NaCl","KCl","HCl","KOH","NaOH","CaSO4","MgSO4","MgCl","CaCl","CaCl2","MgOH","CaOH")
constit <- c("Na","SO4","H","SO4","K","SO4","Na","Cl","K","Cl","H","Cl","K","OH","Na","OH","Ca","SO4","Mg","SO4","Mg","Cl","Ca","Cl","Ca","Cl","Cl","Mg","OH","Ca","OH")
numz <- c(rep(2,12),3,2,2)
charges <- c(-1,-1,-1,0,0,0,0,0,0,0,1,1,0,1,1)
as <- c(4,4,4,0,0,0,0,0,0,0,8,6,0,8,6)
kt <- subset(ktable, T== Tc)
ks <- sapply(1:length(products),function(k) kt[[which(colnames(kt) == products[k])]])

if(!is.null(exprod)){
	products <- c(products,exprod)
	constit <- c(constit,exconstit)
	numz <- c(numz,exnumz)
	charges <- c(charges,excharges)
	as <- c(as,exa)
	ks <- c(ks,exK)
}
prodz <- prods(names=products,number=numz,species=constit,K=ks,a=as)

res <- chemsolve_generic(species=spec,conc=concz,a=speca,charges=speccharges,prod=prodz,Tc=Tc,start=start,maxitr=maxitr)

return(res)
 

}