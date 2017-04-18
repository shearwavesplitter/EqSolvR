#' @title EqSolv
#' @description Find the equilibrium concentrations
#' @param Tc Temperature (degrees C)
#' @param Nat Sodium
#' @param Kt Potassium
#' @param Clt Chlorine
#' @param SO4t Sulphate
#' @param start Initial guess
#' @param maxitr Maximum number of iterations
#' @return A list containing the concentrations, gamma values, and pH at equilibrium
#' @export
htes <- function(Tc=300,Nat=0.2,Kt=0.2,Clt=0.4,SO4t=0.2,Cat=0.1,Mgt=0.1,start=c(0.00001,0.00001,0.15,0.15,0.15,0.104756881,0.05,0.05),maxitr=100){


	Tk <- Tc+273.15

	if (Tc < 300){
		Bdot <- -0.00000030611*Tk*Tk+0.0001027*Tk+0.038646;-0.0000028612*Tk*Tk+0.00094299*Tk-0.026621
	} else {
		Bdot <- 0
	}

	z1 <- 1; z2 <- -1; z3 <- -1; z4 <- 1; z5 <- 1; z6 <- -2;z7 <- -1; z8 <- 2; z9 <- -1; z10 <- -1; z11 <- -1; z12 <-0; z13 <-0; z14 <-0; z15 <-0; z16 <- 0; z17 <-0; z18 <-0; z19 <-1; z20 <-1; z21 <-0; z22 <-1; z23 <-1;
	z <- c(z1,z2,z3,z4,z5,z6,z7,z8,z9,z10,z11,z12,z13,z14,z15,z16,z17,z18,z19,z20,z21,z22,z23)

	#Constants
	tab <- ktable
	temps <- tab[[1]]

	if(Tc %in% tab[[1]]){
		wh <- which(temps == Tc)

		Kw <- 10^tab[[2]][wh]
		K_2 <- 10^tab[[3]][wh]
		K_3 <- 10^tab[[4]][wh]
		K_4 <- 10^tab[[5]][wh]
		K_5 <- 10^tab[[6]][wh]
		K_6 <- 10^tab[[7]][wh]
		K_7 <- 10^tab[[8]][wh]
		K_8 <- 10^tab[[9]][wh]
		K_9 <- 10^tab[[10]][wh]
		K_10 <- 10^tab[[11]][wh]
		K_11 <- 10^tab[[12]][wh] 
		K_12 <- 10^tab[[13]][wh]
		K_13 <- 10^tab[[14]][wh]
		K_14 <- 10^tab[[15]][wh]
		K_15 <- 10^tab[[16]][wh]
		K_16 <- 10^tab[[17]][wh]
	} else {
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

		Kw <- 10^as.numeric(predict(Kws,newdata=data.frame(temps=Tc)))
		K_2 <- 10^as.numeric(predict(K_2s,newdata=data.frame(temps=Tc)))
		K_3 <- 10^as.numeric(predict(K_3s,newdata=data.frame(temps=Tc)))
		K_4 <- 10^as.numeric(predict(K_4s,newdata=data.frame(temps=Tc)))
		K_5 <- 10^as.numeric(predict(K_5s,newdata=data.frame(temps=Tc)))
		K_6 <- 10^as.numeric(predict(K_6s,newdata=data.frame(temps=Tc)))
		K_7 <- 10^as.numeric(predict(K_7s,newdata=data.frame(temps=Tc)))
		K_8 <- 10^as.numeric(predict(K_8s,newdata=data.frame(temps=Tc)))
		K_9 <- 10^as.numeric(predict(K_9s,newdata=data.frame(temps=Tc)))
		K_10 <- 10^as.numeric(predict(K_10s,newdata=data.frame(temps=Tc)))
		K_11 <- 10^as.numeric(predict(K_11s,newdata=data.frame(temps=Tc)))
		K_12 <- 10^as.numeric(predict(K_12s,newdata=data.frame(temps=Tc)))
		K_13 <- 10^as.numeric(predict(K_13s,newdata=data.frame(temps=Tc)))
		K_14 <- 10^as.numeric(predict(K_14s,newdata=data.frame(temps=Tc)))
		K_15 <- 10^as.numeric(predict(K_15s,newdata=data.frame(temps=Tc)))
		K_16 <- 10^as.numeric(predict(K_16s,newdata=data.frame(temps=Tc)))
#plot fit
#v1 <- as.numeric(predict(Kws,newdata=data.frame(temps=seq(300,400,1))))
#v2 <- seq(300,400,1)
#plot(v2,v1,type="l")
#points(tab[[1]],tab[[2]])

	}

	k <- c(Kw,K_2,K_3,K_4,K_5,K_6,K_7,K_8,K_9,K_10,K_11,K_12,K_13,K_14,K_15,K_16)
	k <- log10(k)
	k <- t(k)
 	k <- as.data.frame(k)
	colnames(k) <- names(ktable[2:17])

	#initial g
	g <- rep(1,23)
	#intial guess
	X_1 <- start[1]
	X_2 <- start[2]
	X_3 <- start[3]
	X_4 <- start[4]
	X_5 <- start[5]
	X_6 <- start[6]
	X_7 <- start[7]
	X_8 <- start[8]
	X_9 <- g[5]*X_5*g[6]*X_6/(g[9]*K_7)
	X_10 <- g[1]*X_1*g[6]*X_6/(g[10]*K_8)
	X_11 <- g[4]*X_4*g[6]*X_6/(g[11]*K_9)
	X_12 <- g[5]*X_5*g[3]*X_3/K_2
	X_13 <- g[4]*X_4*g[3]*X_3/K_3
	X_14 <- g[1]*X_1*g[3]*X_3/K_4
	X_15 <- g[4]*X_4*g[2]*X_2/K_5
	X_16 <- g[5]*X_5*g[2]*X_2/K_6
	X_17 <- X_7*g[7]*X_6*g[6]/K_15
	X_18 <- X_8*g[8]*X_6*g[6]/K_16
	X_19 <- X_8*g[8]*X_3*g[3]/(g[19]*K_10)
	X_20 <- X_7*g[7]*X_3*g[3]/(g[20]*K_11)
	X_21 <- X_7*g[7]*X_3^2*g[3]^2/K_12
	X_22 <- X_8*g[8]*X_2*g[2]/(g[22]*K_13)
	X_23 <- X_7*g[7]*X_2*g[2]/(g[23]*K_14)

	#variables for loop
	cr <- 1
	it <- 0
	while (cr != 23){
		it <- it +1 #iterations 
		if(it > maxitr){stop("Max iterations reached")}
		if(exists('ss')){
			X_1 <- ss$root[1]
			X_2 <- ss$root[2]
			X_3 <- ss$root[3]
			X_4 <- ss$root[4]
			X_5 <- ss$root[5]
			X_6 <- ss$root[6]
			X_7 <- ss$root[7]
			X_8 <- ss$root[8]
			cs <- c(X_1,X_2,X_3,X_4,X_5,X_6,X_7,X_8,X_9,X_10,X_11,X_12,X_13,X_14,X_15,X_16,X_17,X_18,X_19,X_20,X_21,X_22,X_23)+1
			
			cr <- cs/(xs+1)
			#cr[is.nan(cr)] <- 1
			cr <- sum(cr)
		}

		xs <- c(X_1,X_2,X_3,X_4,X_5,X_6,X_7,X_8,X_9,X_10,X_11,X_12,X_13,X_14,X_15,X_16,X_17,X_18,X_19,X_20,X_21,X_22,X_23)
		I <- 0.5*sum(xs*z^2)
		a <- c(9,4,3.5,3,4,4,6,8,4,4,4,0,0,0,0,0,0,0,8,6,0,8,6)
		B <- 0.00025587*Tk+0.2487
		A <- 0.00000000041775*((Tk)^4)-0.00000069009*((Tk)^3)+0.00042737*((Tk)^2)-0.11558*(Tk)+11.975
		logg <- ((A*z^2*sqrt(I))/(1+a*B*sqrt(I)))+Bdot*I
		g <- 10^(-logg)
		#Recalculate  X_9-23
	X_9 <- g[5]*X_5*g[6]*X_6/(g[9]*K_7)
	X_10 <- g[1]*X_1*g[6]*X_6/(g[10]*K_8)
	X_11 <- g[4]*X_4*g[6]*X_6/(g[11]*K_9)
	X_12 <- g[5]*X_5*g[3]*X_3/K_2
	X_13 <- g[4]*X_4*g[3]*X_3/K_3
	X_14 <- g[1]*X_1*g[3]*X_3/K_4
	X_15 <- g[4]*X_4*g[2]*X_2/K_5
	X_16 <- g[5]*X_5*g[2]*X_2/K_6
	X_17 <- X_7*g[7]*X_6*g[6]/K_15
	X_18 <- X_8*g[8]*X_6*g[6]/K_16
	X_19 <- X_8*g[8]*X_3*g[3]/(g[19]*K_10)
	X_20 <- X_7*g[7]*X_3*g[3]/(g[20]*K_11)
	X_21 <- 2*X_7*g[7]*X_3^2*g[3]^2/K_12
	X_22 <- X_8*g[8]*X_2*g[2]/(g[22]*K_13)
	X_23 <- X_7*g[7]*X_2*g[2]/(g[23]*K_14)

		xs <- c(X_1,X_2,X_3,X_4,X_5,X_6,X_7,X_8,X_9,X_10,X_11,X_12,X_13,X_14,X_15,X_16,X_17,X_18,X_19,X_20,X_21,X_22,X_23)
		I <- 0.5*sum(xs*z^2)
		B <- 0.00025587*Tk+0.2487
		A <- 0.00000000041775*((Tk)^4)-0.00000069009*((Tk)^3)+0.00042737*((Tk)^2)-0.11558*(Tk)+11.975
		logg <- ((A*z^2*sqrt(I))/(1+a*B*sqrt(I)))+Bdot*I	
		g <- 10^(-logg)

#Super function 1
sfn1 <- function(eqn,charge,g,name){
	ex <- list()
	ex$eqn <- eqn
	ex$charge <- charge
	ex$g <- g
	ex$name <- name
return(ex)
}



	H <- sfn1('x[1]',1,g[1],'H') #H+
	OH <- sfn1('x[2]',-1,g[2],'OH') #OH-
	Cl <- sfn1('x[3]',-1,g[3],'Cl') #Cl-
	K <- sfn1('x[4]',1,g[4],'K') #K+
	Na <- sfn1('x[5]',1,g[5],'Na') #NA+
	SO4 <- sfn1('x[6]',-2,g[6],'SO4') #SO4=
	Ca <- sfn1('x[7]',2,g[7],'Ca') #Ca2+
	Mg <- sfn1('x[8]',2,g[8],'Mg') #Mg2+

sfn2 <- function(species1,species2,number1,number2,gamma,k,name){
	eqn <- paste0(species1$g,"*",species1$eqn,"*",species2$g,"*",species2$eqn,"/(",gamma,"*",k,")")
	charge <- species1$charge*number1+species2$charge*number2
	ex2 <- list()
	ex2$eqn <- eqn
	ex2$charge <- charge
	ex2$name <- name
return(ex2)
}


	NaSO4 <- sfn2(Na,SO4,1,1,g[9],K_7,"NaSO4") #NaSO4-
	HSO4 <- sfn2(H,SO4,1,1,g[10],K_8,"HSO4") #HSO4-
	KSO4 <- sfn2(K,SO4,1,1,g[11],K_9,"KSO4") #KSO4-
	NaCl <- sfn2(Na,Cl,1,1,1,K_2,"NaCl")
	KCl <- sfn2(K,Cl,1,1,1,K_3,"KCl")
	HCl <- sfn2(H,Cl,1,1,1,K_4,"HCl")
	KOH <- sfn2(K,OH,1,1,1,K_5,"KOH")
	NaOH <- sfn2(Na,OH,1,1,1,K_6,"NaOH")
	CaSO4F <- sfn2(Ca,SO4,1,1,g[6],K_15,"CaSO4F")#CaSO4Ḟ
	MgSO4F <- sfn2(Mg,SO4,1,1,g[6],K_16,"MgSO4F") #MgSO4Ḟ
	MgCl <- sfn2(Mg,Cl,1,1,g[19],K_10,"MgCl")  #MgCl+
	CaCl <- sfn2(Ca,Cl,1,1,g[20],K_11,"CaCl") #CaCl+
	CaCl2 <-  sfn2(Ca,Cl,1,2,g[3],K_12,"CaCl")#CaCl2Ḟ
	MgOH <- sfn2(Mg,OH,1,1,g[22],K_13,"MgOH") #MgOH+
	CaOH <- sfn2(Ca,OH,1,1,g[23],K_14,"CaOH") #MgOH+
	
	#f1 <- paste0(NaSO4$charge*NaSO4$eqn,"+",


	f1 <- parse(text='F1=x[1]-x[2]-x[3]+x[4]+x[5]-2*x[6]+2*x[7]+2*x[8]+-g[4]*x[4]*g[6]*x[6]/(g[11]*K_9)-g[5]*x[5]*g[6]*x[6]/(g[9]*K_7)-g[1]*x[1]*g[6]*x[6]/(g[10]*K_8)+x[8]*g[8]*x[3]*g[3]/(g[19]*K_10)+x[7]*g[7]*x[3]*g[3]/(g[20]*K_11)+x[8]*g[8]*x[2]*g[2]/(g[22]*K_13)+x[7]*g[7]*x[2]*g[2]/(g[23]*K_14)')
	f2 <- parse(text='F2=x[5]+g[5]*x[5]*g[3]*x[3]/(K_2)+g[5]*x[5]*g[2]*x[2]/(K_6)+g[5]*x[5]*g[6]*x[6]/(g[9]*K_7) -Nat')	
	f4 <- parse(text='F4=x[3]+ g[5]*x[5]*g[3]*x[3]/(K_2) + g[4]*x[4]*g[3]*x[3]/(K_3)+ g[1]*x[1]*g[3]*x[3]/(K_4)+x[8]*g[8]*x[3]*g[3]/(g[19]*K_10)+x[7]*g[7]*x[3]*g[3]/(g[20]*K_11)+2*x[7]*g[7]*x[3]^2*g[3]^2/K_12-Clt')
	f3 <- parse(text='F3=x[4]+ g[4]*x[4]*g[3]*x[3]/(K_3)+ g[4]*x[4]*g[2]*x[2]/(K_5) + g[4]*x[4]*g[6]*x[6]/(g[11]*K_9) -Kt')
	f5 <- parse(text='F5=g[1]*g[2]*x[1]*x[2]/Kw-1')
	f6 <- parse(text='F6=x[6] + g[4]*x[4]*g[6]*x[6]/(g[11]*K_9) + g[5]*x[5]*g[6]*x[6]/(g[9]*K_7)+g[1]*x[1]*g[6]*x[6]/(g[10]*K_8)+x[7]*g[7]*x[6]*g[6]/K_15+x[8]*g[8]*x[6]*g[6]/K_16 - SO4t')
	f7 <- parse(text='F7=x[7]+x[7]*g[7]*x[6]*g[6]/K_15+x[7]*g[7]*x[3]*g[3]/(g[20]*K_11)+x[7]*g[7]*x[3]^2*g[3]^2/K_12+x[7]*g[7]*x[2]*g[2]/(g[23]*K_14) -Cat')
	f8 <- parse(text='F8=x[8]+x[8]*g[8]*x[6]*g[6]/K_16+x[8]*g[8]*x[3]*g[3]/(g[19]*K_10)+x[8]*g[8]*x[2]*g[2]/(g[22]*K_13) - Mgt')	
	
		vctr <- parse(text='c(eval(f1),eval(f2),eval(f3),eval(f4),eval(f5),eval(f6),eval(f7),eval(f8))')

				#model <- function(x) c(eval(f1),
					#eval(f2),
					#eval(f3),
					#eval(f4),
					#eval(f5),
					#eval(f6),
					#eval(f7),
					#eval(f8))

		model <- function(x) eval(vctr)


	#model <- function(x) c(F1=x[1]-x[2]-x[3]+x[4]+x[5]-2*x[6]+2*x[7]+2*x[8]+-g[4]*x[4]*g[6]*x[6]/(g[11]*K_9)-g[5]*x[5]*g[6]*x[6]/(g[9]*K_7)-g[1]*x[1]*g[6]*x[6]/(g[10]*K_8)+x[8]*g[8]*x[3]*g[3]/(g[19]*K_10)+x[7]*g[7]*x[3]*g[3]/(g[20]*K_11)+x[8]*g[8]*x[2]*g[2]/(g[22]*K_13)+x[7]*g[7]*x[2]*g[2]/(g[23]*K_14),
					#F2=x[5]+g[5]*x[5]*g[3]*x[3]/(K_2)+g[5]*x[5]*g[2]*x[2]/(K_6)+g[5]*x[5]*g[6]*x[6]/(g[9]*K_7) -Nat,
					#F3=x[4]+ g[4]*x[4]*g[3]*x[3]/(K_3)+ g[4]*x[4]*g[2]*x[2]/(K_5) + g[4]*x[4]*g[6]*x[6]/(g[11]*K_9) -Kt,
					#F4=x[3]+ g[5]*x[5]*g[3]*x[3]/(K_2) + g[4]*x[4]*g[3]*x[3]/(K_3)+ g[1]*x[1]*g[3]*x[3]/(K_4)+x[8]*g[8]*x[3]*g[3]/(g[19]*K_10)+x[7]*g[7]*x[3]*g[3]/(g[20]*K_11)+2*x[7]*g[7]*x[3]^2*g[3]^2/K_12-Clt,
					#f,
					#F6=x[6] + g[4]*x[4]*g[6]*x[6]/(g[11]*K_9) + g[5]*x[5]*g[6]*x[6]/(g[9]*K_7)+g[1]*x[1]*g[6]*x[6]/(g[10]*K_8)+x[7]*g[7]*x[6]*g[6]/K_15+x[8]*g[8]*x[6]*g[6]/K_16 - SO4t,
					#F7=x[7]+x[7]*g[7]*x[6]*g[6]/K_15+x[7]*g[7]*x[3]*g[3]/(g[20]*K_11)+x[7]*g[7]*x[3]^2*g[3]^2/K_12+x[7]*g[7]*x[2]*g[2]/(g[23]*K_14) -Cat,
					#F8=x[8]+x[8]*g[8]*x[6]*g[6]/K_16+x[8]*g[8]*x[3]*g[3]/(g[19]*K_10)+x[8]*g[8]*x[2]*g[2]/(g[22]*K_13) - Mgt)
		ss <- multiroot(f=model,start=c(X_1,X_2,X_3,X_4,X_5,X_6,X_7,X_8),positive=TRUE,maxiter=maxitr)

	}

	f <- c(ss$root,X_9,X_10,X_11,X_12,X_13,X_14,X_15,X_16,X_17,X_18,X_19,X_20,X_21,X_22,X_23)
	names <- c("H+","OH-","Cl-","K+","Na+","SO4=","Ca2+","Mg2+","NaSO4-","HSO4-","KSO4-","NaCl","KCl","HCl","KOH","NaOH","CaSO4F","MgSO4F","MgCl+","CaCl+","CaCl2F","MgOH","CaOH+")
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
	res$logK <- k
	res$estim.precis <- ss$estim.precis
	res$f.root <- ss$f.root
	print(paste0("Iterations: ",it))

	nat <- f[5]+f[9]+f[12]+f[16]
	kt <- f[4]+f[11]+f[13]+f[15]
	clt <- f[3]+f[12]+f[13]+f[14]+f[19]+f[20]+f[21]*2

	so4t <- f[6]+f[9]+f[10]+f[11]+f[17]+f[18]
	cat <- f[7]+f[17]+f[20]+f[21]+f[23]
	mgt <- f[8]+f[18]+f[19]+f[22]

	decimalplaces <- function(x) {
 	   if ((x %% 1) != 0) {
   	     nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed=TRUE)[[1]][[2]])
   	 } else {
   	     return(0)
   	 }
	}
	
 	nat <- as.numeric(round(nat,decimalplaces(Nat)+1))
 	kt <- as.numeric(round(kt,decimalplaces(Kt)+1))
 	clt <- as.numeric(round(clt,decimalplaces(Clt)+1))
 	so4t <- as.numeric(round(so4t,decimalplaces(SO4t)+1))
 	cat <- as.numeric(round(cat,decimalplaces(Cat)+1))
 	mgt <- as.numeric(round(mgt,decimalplaces(Mgt)+1))


	if(nat != Nat | kt != Kt | clt != Clt | so4t != SO4t | cat != Cat | mgt != Mgt){print("WARNING: Inconsistent concentrations")}



return(res)
}



