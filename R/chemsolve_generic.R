#' @title Mass balance and charge solver for general cases
#' @description Mass balance and charge balance solver for chemical equilibria. This is the workhorse for chemsolve and requires the prods function
#' @param species Chemical symbols of the basis species
#' @param conc Total concentrations of the basis species (mol/kg)
#' @param a Ion size parameters for the basis species
#' @param prod Dataframe detailing the derived species (Output from prods function)
#' @param A A value 
#' @param B B value 
#' @param Bdot Bdot value is zero by default
#' @param start Initial guess for the calculated equilibrium concentration of the basis species (in the same order as the solvent and then the species vectors)
#' @param maxitr Maximum number of iterations
#' @param solvent Symbols for solvent species (should not be changed for water)
#' @param solvcharge Charges for solvent species (should not be changed for water)
#' @param solva Ion size parameters (should not be changed for water)
#' @param Ksoln log K of the solvent
#' @param bal The species to charge balance against (e.g. the default of "Cl") or NULL for none
#' @return A list containing the concentrations, gamma values, and pH at equilibrium
#' @details A generic function to add any basis species, product species or if the log K/temperature range need to be extended. Requires all parameters (e.g. log K at the given temperature). The temperature is indirectly set through the log K, A, B & Bdot values. These parameters need to be reinitialised each time, together with reactants and products, if a calculation across a range of temperatures is required.  A useful upgrade would be to carry over defaults and use a lookup table and interpolation to initialise the parameters across a range of temperatures. This is similar to the wrapper function except that there the defaults are built in and cannot be changed by the casual user. 
#' For more details see chemsolve documentation.
#' @importFrom rootSolve multiroot
#' @export
#' @examples
#' ## Define all the product species including KHSO4° FeSO4° and FeCl+ 
#' ## which are not included as default species in chemsolve 
#' ## LogK are at 400°C and 0.5kb. Firsly parameters for the products are defined
#' prd=c("NaSO4","HSO4","KSO4","NaCl","KCl","HCl","KOH","NaOH","CaSO4",
#' "MgSO4","MgCl","CaCl","CaCl2","MgOH","CaOH","FeCl","FeSO4","KHSO4")
#' prdnms=c(2,2,2,2,2,2,2,2,2,2,2,2,3,2,2,2,2,3)
#' prdconstit=c("Na","SO4","H","SO4","K",
#' "SO4","Na","Cl","K","Cl","H","Cl","K","OH","Na","OH","Ca","SO4","Mg","SO4","Mg","Cl",
#' "Ca","Cl","Ca","Cl","Cl","Mg","OH","Ca","OH","Fe","Cl","Fe","SO4","K","H","SO4")
#' prdK=c(-3.549,-7.444,-3.899,-1.737,-1.236,-2.689,-1.446,-1.164,-6.173,-6.014,-3.18,
#' -3.692,-4.783,-6.149,-5.635,-5.745,-3.814, -8.701)
#' prda=c(4,4,4,0,0,0,0,0,0,0,8,6,0,8,6,6,0,0)
#' ## product input is created with the prods function
#' products <- prods(names=prd,number=prdnms,species=prdconstit,K=prdK,a=prda)	
#' 														
#' ## Starting species are defined for chemsolve
#' chspec=c("Na", "K", "Cl", "SO4", "Ca", "Mg", "Fe")
#' chconc=c(0.4, 0.2, 0.8, 0.2, 0.1, 0.1,0.1)
#' cha=c(4, 3,3.5, 4, 6, 8,6)
#' chc=c(1, 1, -1, -2, 2, 2,2)
#' 
#' ## chemsolve is now run with the previously defined products
#' ## Defaults include water as the solvent and a charge balance against Cl
#' ## A & B are at 400°C and 0.5kb
#' chemsolve_generic(species = chspec, conc = chconc, a = cha, charges = chc, A = 1.8789,B = 0.423, Bdot = 0, 
#' start = c(1e-06, 1e-05, 0.3, 0.1, 0.3,0.01, 0.001, 0.02,1e-8),  prod = products)

chemsolve_generic <- function(solvent=c("H","OH"),solvcharge=c("1","-1"),solva=c("9","4"),Ksoln=-11.356,species, conc,a,charges,prod,A,B,Bdot=0,start,maxitr=100,bal="Cl"){

    lsolv <- length(solvent)
    lstart <- length(start)
    lspecies <- length(species)
    if(lstart != (lsolv+lspecies)){stop('start vector must be equal to length of solvent plus species')}

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


			x[1:length(ss$root)] <- ss$root

			cs <- x+1

			cr <- cs[1:length(xs)]/(xs+1)
			cr <- sum(cr)
		}

		xs <- x
		I <- 0.5*sum(xs*full$charge^2)
		a <- full$a
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
		logg <- ((A*full$charge^2*sqrt(I))/(1+a*B*sqrt(I)))+Bdot*I	
		g <- 10^(-logg)


	model <- function(x) eval(vctrp)
		ss <- multiroot(f=model,start=x[1:(2+rl)],positive=TRUE,maxiter=maxitr)

	}

	f <- c(ss$root,x[(rl+3):tl])
	names <- paste0(full$species,"(",full$charge,")")
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
	res$A <- A
	res$B <- B
	res$Bdot <- Bdot
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
