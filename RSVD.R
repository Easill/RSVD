library(rsvd)
?rsvd

#Simulate a general matrix with 1000 rows and 1000 columns
vy= rnorm(1000*1000,0,1)
y= matrix(vy,1000,1000,byrow=TRUE)

#Compute the randSVD for the first hundred components of the matrix y and measure the time
start.time <- Sys.time()
prova1= rsvd(y,k=1000)
Sys.time()- start.time

#Compare with a classical SVD
start.time <- Sys.time()
prova2= svd(y)
Sys.time()- start.time

P1 <- prova1$u %*% diag(prova1$d) %*% t(prova1$v)

print(100 * norm( y - P1, 'F') / norm( y,'F')) # percentage error

P2 <- prova2$u %*% diag(prova2$d) %*% t(prova2$v)

print(100 * norm( y - P2, 'F') / norm( y,'F')) # percentage error
