from sklearn import datasets
from sklearn import linear_model
from sklearn import metrics
from numpy import *

iris = datasets.load_iris()
X, Y = iris.data[:,1:4], iris.data[:,0]

clf = linear_model.LinearRegression()
clf.fit(X,Y)

MSE = mean((Y - clf.predict(X))**2)
var_est = MSE * diag(linalg.pinv(dot(X.T,X)))
SE_est = sqrt(var_est)
print(SE_est)

## Bootstrapping
from numpy import *
from random import *
N = len(X)
B = 1000
sim = list()
num_coef = 3
	
for r in range(B):
	newX = list()
	for i in range(N):
		newX.append(randint(0,N-1))
	newX = X[newX,]
	clf = linear_model.LinearRegression()
	clf.fit(newX,Y)
	sim.append(clf.coef_)

se = list()
for j in range(num_coef):
	test = list()
	for i in range(len(sim)):
		test.append(sim[i][j])
	se.append(sqrt(var(test)))