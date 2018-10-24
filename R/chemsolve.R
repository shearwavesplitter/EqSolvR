#' @title Mass balance and charge solver
#' @description Mass balance and charge balance solver for chemical equilibria at 0.5 kb. This is a wrapper function for chemsolve_generic and prods.
#' @param Tc Temperature (degrees C - between 300 and 400)
#' @param Nat Sodium concentration (mol/kg); total 
#' @param Kt Potassium (mol/kg); total
#' @param Clt Chloride (mol/kg); total
#' @param SO4t Sulphate (mol/kg); total
#' @param Cat Calcium (mol/kg); total
#' @param Mgt Magnesium (mol/kg); total
#' @param start Initial guess for the calculated equilbrium concentrations of the basis species (in the order of H+, OH-, Na+, K+, Cl-, SO42+, Ca2+, Mg2+).
#' @param maxitr Maximum number of iterations
#' @param exprod A vector of the names of the additional complex(es)
#' @param exconstit  A vector of the chemical symbol names of the the basis species that are constitute each of the the additional complexes
#' @param exnumz A vector of the stiochiometry given by the equilibrium reaction for each of the additional complexes
#' @param excharges A vector of the charge of the additional complex species
#' @param exK A vector of the log K of the dissociation constants of the additional complex species 
#' @param exa A vector of the ion size parameters for the complexes of the additional complex species 
#' @param bal Species to balance against (Defaults to Cl, can also be set to NULL for none)
#' @details A wrapper for the chemsolve_generic function that allow easy addition of product species.
#' Use the generic function (chemsolve_generic) if new basis species need to be added or if the log K/temperature range is extended (up or down). \cr
#' Normally total moles anions = total moles cations. The charge balance (without any speciation) is adjusted to zero 
#' by balancing against Cl- (default) as otherwise the calculation can be too sensitive to H+. The balancing species can be easily changed to any of the other basis species (or none). \cr
#' Choose reasonable starting values; for H+, OH- and equilibrium concentrations of the basis species. \cr
#' If negative concentrations are calculated,  choose better initial starting values. \cr
#' If bal = NULL and a large excess of choride (or sulphate) is present this will simimulate the addition of acid HCl. In this case use a a lower starting value of H+.\cr
#' To exclude a  basis value set the basis concentration to zero and the concentrations of this and the derived species 
#' will be vanishing small and can be ignored. In the generic version of the program the basis species are simply left out. \cr
#' The basis species are:  Na+, K+, Mg2+, Ca2+, Cl-, SO42-. The default complexes are:  NaCl°, KCl°, HCl°, KOH°, NaOH°, KSO4-, NaSO4-,HSO4-,CaSO4°,MgSO4°, MgCl+,CaCl+,CaCl2°,MgOH+,CaOH+. 
#' Additional complexes based on the existing basis species can easily be added (see example below). \cr
#' To set up a problem in chemsolve or chemsolve_generic it is relatively trivial to set up and manipulate a sequence of commands in Excel and cut and past into R.
#' Complex dissociation constants (log K)  are from SupCrt 92  slop98.dat \url{http://geopig.asu.edu/?q=tools} \cr
#' The Debye_Hückel parameters (A, B) equations are polynomial fits to data at 0.5 kb from tables in Helgeson & Kirkham (1974).\cr
#' Note Bdot is not used.\cr
#' Helgeson H. C. and Kirkham D. H. (1974) Theoretical prediction of the thermodynamic behavior of aqueous electrolytes at high pressures and temperatures: II. Debye-Hückel parameters for activity coefficients and relative partial molar properties American Journal of Science 274, 1199-1261. \cr
#' @return A list containing the concentrations, activity coefficients, and pH at equilibrium
#' @importFrom rootSolve multiroot
#' @export
#' @examples
#' ## Add KHSO4° as an additional complex given the existing list of basis species and 
#' ## calculate the equilibrium concentrations and pH at 400°C.
#' chemsolve(exprod =  c("KHSO4"), exconstit = c("K", "H", "SO4"),
#' exnumz =  c(3), excharges =  c(0), exa =  c(0),  exK =  c(-8.701), bal = "Cl")
#' ## Determine the equilibria at a range of temperatures. 
#' ## Additional complexes can be added per the previous example.
#' temps <- seq(300,400,10) #A vector of temperatures repeating every 10 degrees from 300 to 400
#' ## Na concentration is changed slightly from the default. Further parameters can be added.
#' r <- lapply(temps,chemsolve,Nat=0.45) #Creates a list of the results
#' r[[1]] #Display results from first temperature
#' r[[10]] #Display the results of the 10th temperature
chemsolve <- function(Tc=400,Nat=0.4,Kt=0.2,Clt=0.6,SO4t=0.2,Cat=0.1,Mgt=0.1,start=c(1e-6,1e-5,0.3,0.1,0.3,0.01,0.001,0.02),maxitr=100,exprod=NULL,exconstit=NULL,exnumz=NULL,excharges=NULL,exa=NULL,exK=NULL,bal="Cl") {
spec <- c("Na","K","Cl","SO4","Ca","Mg")
concz <- c(Nat,Kt,Clt,SO4t,Cat,Mgt)
speccharges <- c(1,1,-1,-2,2,2)
speca <- c(4,3,3.5,4,6,8)

ABt <- seq(250,450,25)
As <- c(0.8822,0.9595,1.0529,1.1705,1.3267,1.5464,1.8789,2.4301,3.3553)
Bs <- c(0.3729,0.3787,0.3850,0.3921,0.4004,0.4104,0.4230,0.4386,0.4548)

ABtab <- as.data.frame(cbind(ABt,As,Bs),stringsAsFactors=FALSE)

products <- c("NaSO4","HSO4","KSO4","NaCl","KCl","HCl","KOH","NaOH","CaSO4","MgSO4","MgCl","CaCl","CaCl2","MgOH","CaOH")
constit <- c("Na","SO4","H","SO4","K","SO4","Na","Cl","K","Cl","H","Cl","K","OH","Na","OH","Ca","SO4","Mg","SO4","Mg","Cl","Ca","Cl","Ca","Cl","Cl","Mg","OH","Ca","OH")
numz <- c(rep(2,12),3,2,2)
charges <- c(-1,-1,-1,0,0,0,0,0,0,0,1,1,0,1,1)
as <- c(4,4,4,0,0,0,0,0,0,0,8,6,0,8,6)


## Determine A and B
temps <- ABtab$ABt
if(Tc %in% temps){
Aspec <- ABtab$As[which(temps == Tc)]
Bspec <- ABtab$Bs[which(temps == Tc)]

}else{
	print("Interpolating A and B")
	if (Tc < min(temps)){warning("Caution: Temperature out of bounds for interpolation of A and B")}
	if (Tc > max(temps)){warning("Caution: Temperature out of bounds for interpolation of A and B")}
	Az <- lm(ABtab$As ~ poly(temps,3))
	Bz <- lm(ABtab$Bs ~ poly(temps,3))
	Aspec <- as.numeric(predict(Az,newdata=data.frame(temps=Tc)))
	Bspec <- as.numeric(predict(Bz,newdata=data.frame(temps=Tc)))
}

## Determine K
tab <- ktable
temps <- tab[[1]]
if(Tc %in% tab[[1]]){
Kw <- tab[which(tab[[1]] == Tc),which(colnames(tab) == "H2O")]
}else{
	print("Interpolating constants")
	if (Tc < min(tab[[1]])){warning("Caution: Temperature out of bounds for interpolation for K")}
	if (Tc > max(tab[[1]])){warning("Caution: Temperature out of bounds for interpolation for K")}
	
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

res <- chemsolve_generic(species=spec,conc=concz,a=speca,charges=speccharges,prod=prodz,start=start,maxitr=maxitr,bal=bal,A=Aspec,B=Bspec,Ksoln=Kw)

res$Temperature <- Tc
return(res)
 

}
