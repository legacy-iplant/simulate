require(mvtBinaryEP)
require(getopt)

seperator <- ","
ydata <- "data1.csv"
markerdata <- "data2.csv"
maf <- 0.5
dcor <- 0.5
dsim <- 3
yname <- "length"
dmu <- 11
dbeta <- 8
dsigmay <- 5
dsigmagamma <- 4
dsigmaalpha <- 4

ydata <- read.table(file=ydata,sep=seperator,header=TRUE)
markerdata <- read.table(file=markerdata,sep=seperator,header=TRUE)

#main <- function(
		minor.allele.freq=maf
		cor=dcor
		sim=dsim
		y.data=ydata
		marker.data=markerdata
		beta.pm=c(dmu,dbeta)
		sigma.sq.y.pm=dsigmay
		sigma.sq.gamma.pm=dsigmagamma
		sigma.sq.alpha.pm=dsigmaalpha
		y.name=yname
	
	y.data <- as.data.frame(y.data)
	y <- eval(parse(text=paste("y.data",y.name,sep="$")))
	x.y <- model.matrix(~treat,data=y.data)
	n <- length(y)
	p <- ncol(t(marker.data))
	marker.names <- markerdata[,1]
	line.names <- colnames(markerdata)
	my.cor <- cor(t(markerdata[,2:ncol(markerdata)]))
	colnames(my.cor) <- marker.names
	rownames(my.cor) <- marker.names

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
	data.sim[[1]][[4]] <- sim

	for (i in 2:(sim+1)) {

		x.l <- ep(mu=minor.allele.freq,rho=cor,n=p,nRep=n.line)[[1]]
		colnames(x.l) <- colnames(t(markerdata))
		rownames(x.l) <- rownames(t(markerdata))
	
		gamma <- rnorm(n.line,x.l%*%delta,sqrt(sigma.sq.gamma.pm))
		alpha <- rnorm(n.id,0,sqrt(sigma.sq.alpha.pm))
		eta <- x.y%*%beta.pm + alpha[alpha.id] + gamma[gamma.id]
		y.sim <- rnorm(n,eta,sqrt(sigma.sq.y.pm))
		y.sim[o.y==TRUE] <- NA

		pheno.sim <- y.data
		pheno.sim$length <- y.sim
		#names(pheno.sim)[1] <- "y"
		data.sim[[i]][[1]] <- pheno.sim
		data.sim[[i]][[2]] <- t(x.l)

	}

	#return(data.sim)

#}