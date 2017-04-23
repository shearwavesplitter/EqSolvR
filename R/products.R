#' @title Create prod dataframe
#' @description Creates the dataframe of the derived species for use in chemsolve
#' @param names A vector of names of the species which react to form the basis species
#' @param number A vector of the number of basis consituents for each of the product species given by the equilibrium equation
#' @param species A vector of the chemical symbols of the product species in terms of the basis species
#' @param K A vector of log K values for the product species
#' @param a A vector of ion size parameters for the product species
#' @export
prods <- function(names=c("NaCl","KCl"),number=c(2,2),species=c("Na","Cl","K","Cl"),K=c(-6.68,0.001),a=c(0,0)){
exfull <- list()
	for (i in 1:length(names)){
		ex <- list()
		ex$name <- names[i]
		cnt1 <- sum(number[1:i])-(number[i]-1)
		cnt2 <- sum(number[1:i])
		ex$species <- species[cnt1:cnt2]
		ex$K <- K[i]
		ex$a <- a[i]

		exfull[[length(exfull)+1]] <- ex
	}
return(exfull)
}