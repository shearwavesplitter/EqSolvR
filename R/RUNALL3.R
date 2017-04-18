#' export
htes3 <- function(Tc=300,Nat=0.2,Kt=0.2,Clt=0.4,SO4t=0.2,Cat=0.1,Mgt=0.1,start=c(0.00001,0.00001,0.15,0.15,0.15,0.104756881,0.05,0.05),maxitr=100,exprod=NULL,exconstit=NULL,exnumz=NULL,excharges=NULL,exa=NULL,exK=NULL) {

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

res <- htes2(species=spec,conc=concz,a=speca,charges=speccharges,prod=prodz,Tc=Tc,start=start,maxitr=maxitr)

return(res)
 

}