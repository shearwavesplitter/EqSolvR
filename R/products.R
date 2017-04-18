#' export
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