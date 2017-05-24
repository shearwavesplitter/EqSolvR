#' @title Mass balance and charge solver for general cases
#' @description Mass balance and charge balance solver for chemical equilibria
#' @param species Chemical symbols of the basis species
#' @param conc Total concentrations of the basis species (mol/kg)
#' @param a Ion size parameters for the basis species
#' @param prod Dataframe detailing the derived species
#' @param Tc Temperature (degrees centigrade)
#' @param start Initial guess for the calculated equalibrium concentration of the basis species
#' @param maxitr Maximum number of iterations
#' @param solvent Symbols for solvent species (should not be changed)
#' @param solvcharge Charges for solvent species (should not be changed)
#' @param solva Ion size parameters (should not be changed)
#' @param Ksoln log K of the solvent (should not be changed)
#' @return A list containing the concentrations, gamma values, and pH at equilibrium
#' @details A generic function to add any basis species, product species or if the log K/temperature range need to be extended. Requires all parameters (e.g. log K at the given temperature). 
#' @importFrom rootSolve multiroot
#' @export
#' @examples
#' ## Add H2SO4 as an additional complex given the existing list of basis species
#'
#' ## Define the product species NaCl and KCl
#' products <- prods(names=c("NaCl","KCl"),number=c(2,2),
#' + species=c("Na","Cl","K","Cl"),K=c(-6.68,0.001),a=c(0,0))
#' Run chemsolve with Na, K, and Cl basis species at 300 degrees
#' Chemsolve_generic(species=c("Na","K","Cl"), conc=c(0.2,0.2,0.4),
#' + a=c(4,3,3.5),charges=c(1,1,-1),prod,Tc=300,start=c(0.00001,0.00001,0.15),prod=products)
chemsolve_generic <- function(solvent=c("H","OH"),solvcharge=c("1","-1"),solva=c("9","4"),Ksoln=-10.908,species=c("Na","K","Cl","SO4","Ca","Mg"), conc=c(0.2,0.2,0.4,0.2,0.1,0.1),a=c(4,3,3.5,4,6,8),charges=c(1,1,-1,-2,2,2),prod,Tc=300,start=c(0.00001,0.00001,0.15,0.15,0.15,0.104756881,0.05,0.05),maxitr=100,bal=NULL){
	Tk <- Tc+273.15

	if (Tc < 300){
		Bdot <- -0.00000030611*Tk*Tk+0.0001027*Tk+0.038646;-0.0000028612*Tk*Tk+0.00094299*Tk-0.026621
	} else {
		Bdot <- 0
	}


	## Adjust one thing so that initial concentrations charge balance to zero
	if(!is.null(bal)){
		whb <- which(species == bal)
		cnc2 <- conc[-whb]
		chrg2 <- charges[-whb]
		ccprod <- cnc2*chrg2
		sumprod <- -sum(ccprod)
		newconc <- sumprod/charges[whb]
		if(newconc < 0){print("Negative concentration is not possible: Continuing with original concentration");print("Consider changing balance species");warning("Initial concentrations not charge balanced")}else{
		print(paste0("Charge balance on ",species[whb]))
		print(paste0("New ",species[whb]," initial concentration of ",newconc))
		conc[whb] <- newconc
		}
	}
	##

	pl <- length(prod)	
	rl <- length(species)
	tl <- pl+rl+2


	k <- 1:2
	vecsolv <- sapply(k,function (i) c(solvent[i],solvcharge[i],solva[i],paste0('x[',i,']'),paste0('g[',i,']')))
	vecsolv <- t(vecsolv)
	vec <- as.data.frame(vecsolv)


	k <- 1:rl
	vec <- sapply(k,function (i) c(species[i],charges[i],a[i],paste0('x[',i+2,']'),paste0('g[',i+2,']')))
	vec <- t(vec)
	vec <- as.data.frame(vec)

	base <- rbind(vecsolv,vec)
	colnames(base) <- c("species","charge","a","eqn","gamma")

	for (i in 1:pl){
		pr <- prod[[i]]
		spec <- pr$species
		unspec <- unique(spec)
		chrg <- 0
		for (j in 1:length(unspec)){
			fr <- unspec[j]
			num <- length(which(fr == spec))
			wh <- which(fr == base$species)
			basex <- base$eqn[wh]
			baseg <- base$gamma[wh]
			basech <- as.numeric(as.character(base$charge[wh]))*as.numeric(as.character(num))
			eqnp1 <- paste0(baseg,"^",num,"*",basex,"^",num)
			chrg <- chrg + basech
			if(j ==1){fulleqn <- eqnp1}else{fulleqn <- paste0(fulleqn,"*",eqnp1)}

			}
			
			fulleqn <- paste0(fulleqn,"/(g[",i+rl+2,"]*",10^pr$K,")")
			
			vecprod <- cbind(pr$name,chrg,pr$a,fulleqn,paste0("g[",i+rl+2,"]"))
			if(i == 1){prodlist <- vecprod}else{prodlist <- rbind(prodlist,vecprod)}
			
			
	}
	prodlist <- as.data.frame(prodlist)
	colnames(prodlist) <- c("species","charge","a","eqn","gamma")


	#initial g
	g <- rep(1,tl)
	#intial guess
	x <- start 
	for (i in 1:pl){
		eq <- parse(text=as.character(prodlist$eqn[i]))
		guess <- eval(eq)
		x <- c(x,guess)
	}

	full <- rbind(base,prodlist)
		

	### Form equations

	#Charge balance
	for (i in 1:length(full$species)){
		crg <- full$charge[i]
		eqn <- full$eqn[i]
		baleqn <- paste0(crg,"*",eqn)
		if (i ==1){fullbaleqn <- baleqn}else{fullbaleqn <- paste0(fullbaleqn,"+",baleqn) }

	}

	eqbal <- parse(text=paste0("F1=",fullbaleqn))


	#The water one
	soleqn <- parse(text=paste0("F2=",full$gamma[1],"*",full$eqn[1],"*",full$gamma[2],"*",full$eqn[2],"/",10^Ksoln,"-1"))

	#Balances
	for (i in 1:rl){	
		spec <- species[i]
		cnc <- conc[i]
		j <- i+2
		wh <- which(full$species == spec)
		first <- full$eqn[wh]
		for(l in 1:length(prod)){
			pr <- prod[[l]]
			tf <- spec %in% pr$species
			tfn <- length(which(pr$species == spec))
			nam <- pr$name

			if(tf){
				wh2 <- which(full$species == nam)
				eqcon <- paste0(tfn,"*",full$eqn[wh2])
				first <- paste0(first,"+",eqcon)
			}
			
		}
		last <- paste0("F",j,"=",first,"-",cnc)
		last <- as.data.frame(last)
		if(i ==1){coneqls <- last}else{coneqls <- rbind(coneqls,last)}
	}
	coneqls$last <- as.character(coneqls$last)
	for(i in 1:length(coneqls$last)){
		assign(paste0("f",i+2),parse(text=coneqls$last[i]))

	}
	
	f1 <- eqbal
	f2 <- soleqn

	
	vctr <- 'c('
	for(i in 1:(rl+2)){
		if (i == 1){vctr <- paste0(vctr,"eval(f",i,")")}else{vctr <- paste0(vctr,",eval(f",i,")")}
		#print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
		#print(get(paste0("f",i)))
	}

	vctr <- paste0(vctr,")")
	vctrp <- parse(text=vctr)
	
	full$charge <- as.numeric(as.character(full$charge))
	full$a <- as.numeric(as.character(full$a))

	#variables for loop
	cr <- 1
	it <- 0
	while (cr != tl){
		it <- it +1 #iterations 
		if(it > maxitr){stop("Max iterations reached")}
		if(exists('ss')){
			
			x[1:(rl+2)] <- ss$root

			cs <- x+1
			
			cr <- cs/(xs+1)
			cr <- sum(cr)
		}

		xs <- x
		I <- 0.5*sum(xs*full$charge^2)
		a <- full$a
		B <- 0.00025587*Tk+0.2487
		A <- 0.00000000041775*((Tk)^4)-0.00000069009*((Tk)^3)+0.00042737*((Tk)^2)-0.11558*(Tk)+11.975
		logg <- ((A*full$charge^2*sqrt(I))/(1+a*B*sqrt(I)))+Bdot*I
		g <- 10^(-logg)
		#Recalculate  X_9-23

		for (i in 1:pl){
			eq <- parse(text=as.character(prodlist$eqn[i]))
			guess <- eval(eq)
			x[i+rl+2] <- guess
		}

		xs <- x
		I <- 0.5*sum(xs*full$charge^2)
		B <- 0.00025587*Tk+0.2487
		A <- 0.00000000041775*((Tk)^4)-0.00000069009*((Tk)^3)+0.00042737*((Tk)^2)-0.11558*(Tk)+11.975
		logg <- ((A*full$charge^2*sqrt(I))/(1+a*B*sqrt(I)))+Bdot*I	
		g <- 10^(-logg)


	model <- function(x) eval(vctrp)
		ss <- multiroot(f=model,start=x[1:(2+rl)],positive=TRUE,maxiter=maxitr)

	}

	f <- c(ss$root,x[(rl+3):tl])
	names <- full$species
	f <- t(f)
	f <- as.data.frame(f)
	colnames(f) <- names
	g <- t(g)
	g <- as.data.frame(g)
	colnames(g) <- names
	res <- list()
	res$conc <- f
	res$gamma <- g
	res$pH <- -log10(as.numeric(f[1]*as.numeric(g[1])))
	#res$logK <- k
	res$estim.precis <- ss$estim.precis
	res$f.root <- ss$f.root
	print(paste0("Iterations: ",it))

	#nat <- f[5]+f[9]+f[12]+f[16]
	#kt <- f[4]+f[11]+f[13]+f[15]
	#clt <- f[3]+f[12]+f[13]+f[14]+f[19]+f[20]+f[21]*2

	#so4t <- f[6]+f[9]+f[10]+f[11]+f[17]+f[18]
	#cat <- f[7]+f[17]+f[20]+f[21]+f[23]
	#mgt <- f[8]+f[18]+f[19]+f[22]

	decimalplaces <- function(x) {
 	   if ((x %% 1) != 0) {
   	     nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed=TRUE)[[1]][[2]])
   	 } else {
   	     return(0)
   	 }
	}
	
 	#nat <- as.numeric(round(nat,decimalplaces(Nat)+1))
 	#kt <- as.numeric(round(kt,decimalplaces(Kt)+1))
 	#clt <- as.numeric(round(clt,decimalplaces(Clt)+1))
 	#so4t <- as.numeric(round(so4t,decimalplaces(SO4t)+1))
 	#cat <- as.numeric(round(cat,decimalplaces(Cat)+1))
 	#mgt <- as.numeric(round(mgt,decimalplaces(Mgt)+1))


	#if(nat != Nat | kt != Kt | clt != Clt | so4t != SO4t | cat != Cat | mgt != Mgt){print("WARNING: Inconsistent concentrations")}



return(res)
}
