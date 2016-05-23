library("glmnet")
source("admm_lasso.r")

m <- 1500
n <- 5000
p <- 100

rho <- 1
alpha <- 1.8

x0 <- Matrix(0, n , 1)
x0[runif(p, 1, n)] <- rnorm(p) 
A <- Matrix(rnorm(m*n), nrow=m, ncol=n)
A <- A %*% Diagonal(n, 1/sqrt(colSums(A^2)))
b <- A %*% x0 + Matrix(rnorm(m, 0, sqrt(0.001)))

lam_max <- max(abs(t(A)%*%b))
lam <- 0.01*lam_max

res <- factor(A, rho)
L <- res$L 
U <- res$U

start <- proc.time()
init_x <- Matrix(0, n , 1)
res <- admm_lasso(A,b,init_x,L,U,lam,rho,alpha,F)
elapsed <- proc.time() - start
print(paste("Elapsed time ", elapsed[3], sep=""))

cat("\n")
print(paste("norm(x0)=", sqrt(sum(x0^2)), ", norm(x)=", sqrt(sum(res$x^2)), ", length(x)=",length(res$x), sep=""))

print("glmnet solution")
r <- glmnet(A, b, family="gaussian", nlambda=1, lambda=lam, standardize=FALSE, intercept=FALSE)
print(paste("norm(beta)=", sqrt(sum(r$beta^2)),", length(r$beta)=",length(r$beta),sep=""))


