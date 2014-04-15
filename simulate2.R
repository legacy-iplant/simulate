#setwd("/users/dustin/documents/simulate")
#options(error = quote({dump.frames(to.file=TRUE); q()}))
require(mvtBinaryEP)
require(getopt)
args <- commandArgs(TRUE)

options <- matrix(c("ydata","a",1,"character",
			"markerdata","b",1,"character",
			"maf","c",1,"double",
			"cor","d",1,"double",
			"sim","e",1,"integer",
			"yname","f",1,"character",
			"mu","g",1,"double",
			"beta","h",1,"double",
			"sigmay","i",1,"double",
			"sigmagamma","j",1,"double",
			"sigmaalpha","k",1,"double",
			"sep","l",1,"character"),
			ncol=4,byrow=TRUE)

ret.opts <- getopt(options,args)
seperator <- ret.opts$sep
ydata <- ret.opts$ydata
markerdata <- ret.opts$markerdata
maf <- ret.opts$maf
dcor <- ret.opts$cor
dsim <- ret.opts$sim
yname <- ret.opts$yname
dmu <- ret.opts$mu
dbeta <- ret.opts$beta
dsigmay <- ret.opts$sigmay
dsigmagamma <- ret.opts$sigmagamma
dsigmaalpha <- ret.opts$sigmaalpha

ydata <- read.table(file=ydata,sep=seperator,header=TRUE)
markerdata <- read.table(file=markerdata,sep=seperator,header=TRUE)

main <- function(minor.allele.freq=maf,cor=dcor,
		sim=dsim,y.data=ydata,marker.data=markerdata,
		beta.pm=c(dmu,dbeta),
		sigma.sq.y.pm=dsigmay,
		sigma.sq.gamma.pm=dsigmagamma,
		sigma.sq.alpha.pm=dsigmaalpha,
		y.name=yname) {
	
	y.data <- as.data.frame(y.data)
	y <- eval(parse(text=paste("y.data",y.name,sep="$")))
	x.y <- model.matrix(~treat,data=y.data)
	n <- length(y)
	p <- ncol(t(marker.data))
	print(head(marker.data[,1:6]))

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
		names(pheno.sim)[1] <- "y"
		data.sim[[i]][[1]] <- pheno.sim
		data.sim[[i]][[2]] <- t(x.l)

	}

	return(data.sim)

}

thedata <- main()

for (i in 2:(thedata[[1]][[4]]+1)) {
	write.csv(x=thedata[[i]][[1]],row.names=FALSE,
		file=paste("sim_",i,"_ydata.csv",sep=""))
	write.csv(x=thedata[[i]][[2]],row.names=FALSE,
		file=paste("sim_",i,"_markerdata.csv",sep=""))
}