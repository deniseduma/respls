#!/usr/bin/Rscript

library(elasticnet)

source("enet_func.R")

data(diabetes)
attach(diabetes)
##fit the lasso model (treated as a special case of the elastic net)
print("dim(x)")
print(dim(x))
print("length(y)")
print(length(y))

obj1 <- enet2(x,y,lambda=1)
plot(obj1)
