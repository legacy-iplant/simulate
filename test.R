data(iris)
X <- iris[,1:4]
fit <- lm(Sepal.Length~.,data=X)
print(summary(fit))
MSE <- anova(fit)[[3]][4]
X <- cbind(1,X[,2:4])
X <- as.matrix(X)
print(sqrt(diag(MSE * solve(t(X)%*%X))))

X <- cbind(iris[,1:4])
N = nrow(X)
B = 1000
num.coef <- 4

## Bootstrapping
#bootstrap <- function(model,X,B,num.coef) {
	N <- length(X)
	sim = list()
	for (r in 1:B) {
		newX <- list()
		newX <- sample(1:N,N,replace=TRUE)
		newX <- as.data.frame(X[newX,]) 
		fit <- eval(parse(text=paste(model)))
		sim[[r]] <- fit$coefficients
	}

	se <- list()
	test <- list()
	for (j in 1:4) {
		for (i in 1:B) {
			test[[i]] <- sim[[i]][j]
		}
		test <- unlist(test)
		se <- append(se,sd(test))
	}
	se <- unlist(se)
	return(se)
#}

standard.errors <- bootstrap("lm(Sepal.Width~.,data=newX)",X,1000,4)
print(standard.errors)