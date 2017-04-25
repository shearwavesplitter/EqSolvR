#' @title Mass balance and charge solver
#' @description Mass balance and charge balance solver for chemical equilibria. 
#' @param Tc Temperature (degrees C - between 300 and 400)
#' @param Nat Sodium concentration (mol/kg); total 
#' @param Kt Potassium (mol/kg); total
#' @param Clt Chloride (mol/kg); total
#' @param SO4t Sulphate (mol/kg); total
#' @param Cat Calcium (mol/kg); total
#' @param Mgt Magnesium (mol/kg); total
#' @param start Initial guess for the calculated equalibrium concentration of the basis species
#' @param maxitr Maximum number of iteration
#' @param exprod A vector of the names of the complexes which dissociate to form the basis species
#' @param exconstit  A vector of the chemical symbol names of the the basis species that are the products of the dissociation equilibium for each complex
#' @param exnumz A vector of the stiochiometery given by the equilibrium reaction for each of the complexes
#' @param excharges A vector of the charge of the complex species
#' @param exK A vector of the log K of the dissociation constants
#' @param exa A vector of the ion size paramters for the complexes
#' @details A wrapper for the chemsolve_generic function that allow easy addition of product species. If you want to add additional reactant species (i.e. basis species) then the chemsolve_generic function must be used.
#' The basis species are:  Na+, K+, Mg2+, Ca2+, Cl-, SO42-. The default complexes are:  NaCl°, KCl°, HCl°, KOH°, NaOH°, KSO4-, NaSO4-,HSO4-,CaSO4°,MgSO4°, MgCl+,CaCl+,CaCl2°,MgOH+,CaOH+.Additional complexes based on the existing basis species are easily added. \cr
#' Use the generic function (chemsolve_generic) if new basis species need to be added or if the log K/temperature range is extended (up or down). \cr
#' Charge balance is fixed on H+. \cr
#' Normally total initial moles anions = total moles cations  but excess anions will be balanced by more H+ and vice versa. It is important to choose good initial starting values; for H+, OH- and equilibrium concentrations of the basis species. \cr
#' Complex dissociation constants (Log K)  are from SupCrt 92  slop98.dat \url{http://geopig.asu.edu/?q=tools} \cr
#' The Debye_Hückel parameters (A, B & Bdot) equations are polynomial fits to data from tables in Helgeson (1969)  Helgeson & Kirkham (1974) by Nellie Olsen (Note Bdot not used at temperatures  greater than 300°C). \cr
#' Helgeson H. C. (1969) Thermodynamics of hydrothermal systems at elevated temperatures and pressures. American Journal of Science 267, 729-804. \cr
#' Helgeson H. C. and Kirkham D. H. (1974) Theoretical prediction of the thermodynamic behavior of aqueous electrolytes at high pressures and temperatures: II. Debye-Hückel parameters for activity coefficients and relative partial molar properties American Journal of Science 274, 1199-1261.
#' @return A list containing the concentrations, activity coefficients, and pH at equilibrium
#' @importFrom rootSolve multiroot
#' @export
#' @examples
#' ## Add H2SO4 as an additional complex given the existing list of basis species
#'
#' chemsolve(exprod="H2SO4",exconstit="H","H","SO4",exnumz=3,excharges=0,exa=0,exK=-6)
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


## Determine K
tab <- ktable
temps <- tab[[1]]
if(Tc %in% tab[[1]]){
}else{
	print("Interpolating constants")
	if (Tc < min(tab[[1]])){print("Caution: Temperature out of bounds for interpolation")}
	if (Tc > max(tab[[1]])){print("Caution: Temperature out of bounds for interpolation")}
	
		Kws <- lm(tab[[2]] ~ poly(temps,3))
		K_2s <- lm(tab[[3]] ~ poly(temps,3))
		K_3s <-  lm(tab[[4]] ~ poly(temps,3))
		K_4s <-  lm(tab[[5]] ~ poly(temps,3))
		K_5s <-  lm(tab[[6]] ~ poly(temps,3))
		K_6s <-  lm(tab[[7]] ~ poly(temps,3))
		K_7s <-  lm(tab[[8]] ~ poly(temps,3))
		K_8s <-  lm(tab[[9]] ~ poly(temps,3))
		K_9s <-  lm(tab[[10]] ~ poly(temps,3))
		K_10s <- lm(tab[[11]] ~ poly(temps,3))
		K_11s <- lm(tab[[12]] ~ poly(temps,3))
		K_12s <- lm(tab[[13]] ~ poly(temps,3))
		K_13s <- lm(tab[[14]] ~ poly(temps,3))
		K_14s <- lm(tab[[15]] ~ poly(temps,3))
		K_15s <- lm(tab[[16]] ~ poly(temps,3))
		K_16s <- lm(tab[[17]] ~ poly(temps,3))

		Kw <- as.numeric(predict(Kws,newdata=data.frame(temps=Tc)))
		K_2 <- as.numeric(predict(K_2s,newdata=data.frame(temps=Tc)))
		K_3 <- as.numeric(predict(K_3s,newdata=data.frame(temps=Tc)))
		K_4 <- as.numeric(predict(K_4s,newdata=data.frame(temps=Tc)))
		K_5 <- as.numeric(predict(K_5s,newdata=data.frame(temps=Tc)))
		K_6 <- as.numeric(predict(K_6s,newdata=data.frame(temps=Tc)))
		K_7 <- as.numeric(predict(K_7s,newdata=data.frame(temps=Tc)))
		K_8 <- as.numeric(predict(K_8s,newdata=data.frame(temps=Tc)))
		K_9 <- as.numeric(predict(K_9s,newdata=data.frame(temps=Tc)))
		K_10 <- as.numeric(predict(K_10s,newdata=data.frame(temps=Tc)))
		K_11 <- as.numeric(predict(K_11s,newdata=data.frame(temps=Tc)))
		K_12 <- as.numeric(predict(K_12s,newdata=data.frame(temps=Tc)))
		K_13 <- as.numeric(predict(K_13s,newdata=data.frame(temps=Tc)))
		K_14 <- as.numeric(predict(K_14s,newdata=data.frame(temps=Tc)))
		K_15 <- as.numeric(predict(K_15s,newdata=data.frame(temps=Tc)))
		K_16 <- as.numeric(predict(K_16s,newdata=data.frame(temps=Tc)))
	line <- c(Tc,Kw,K_2,K_3,K_4,K_5,K_6,K_7,K_8,K_9,K_10,K_11,K_12,K_13,K_14,K_15,K_16)
	tab <- rbind(tab,line)
}


	kt <- subset(tab, T== Tc)
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