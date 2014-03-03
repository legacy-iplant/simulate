setwd("/users/dustin/desktop")
load("Data.Sim.50.Rdata")
require(mvtBinaryEP)

data1 <- Data.Sim[[1]][[1]]
data2 <- Data.Sim[[1]][[2]]
data3 <- Data.Sim[[1]][[3]]

generate <- function(minor.allele.freq=0.5,cor=0.5,
		sim=10,y.data=data1,marker.data=data3,
		beta.pm=c(47.71,8.96),
		sigma.sq.y.pm=22.25,
		sigma.sq.gamma.pm=1.44,
		sigma.sq.alpha.pm=11.45) {
	
	y <- y.data$length
	x.y <- model.matrix(~treat,data=y.data)
	n <- length(y)
	p <- ncol(t(marker.data))

	tab.id <- table(y.data$id)
	n.id <- length(tab.id)
	alpha.id <- rep(1:n.id,tab.id)
	tab.line <- table(y.data$line)
	n.line <- length(tab.line)
	gamma.id <- rep(1:n.line,tab.line)
	o.y <- is.na(y)
	delta <- c(10:1,rep(0,p-10))

	data.sim <- vector(mode="list",length=(sim+1))
	data.sim[[1]] <- vector(mode="list",length=3)
	data.sim[[1]][[1]] <- y.data
	data.sim[[1]][[2]] <- x.y
	data.sim[[1]][[3]] <- marker.data

	for (i in 2:(sim+1)) {

		x.l <- ep(mu=minor.allele.freq,rho=cor,n=p,nRep=n.line)[[1]]
		#colnames(x.l) <- colnames(t(data3))
		#rownames(x.l) <- rownames(t(data3))
	
		gamma <- rnorm(n.line,x.l%*%delta,sqrt(sigma.sq.gamma.pm))
		alpha <- rnorm(n.id,0,sqrt(sigma.sq.alpha.pm))
		eta <- x.y%*%beta.pm + alpha[alpha.id] + gamma[gamma.id]
		y.sim <- rnorm(n,eta,sqrt(sigma.sq.y.pm))
		y.sim[o.y==TRUE] <- NA

		pheno.sim <- y.data
		pheno.sim$length <- y.sim
		names(pheno.sim)[1] <- "y"
		data.sim[[i]][[1]] <- pheno.sim
		data.sim[[i]][[2]] <- t(x.l)

	}

	return(data.sim)

}