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


# Importation d'une image

library(raster)
library(rgdal)
photo <- as.matrix(raster("https://husson.github.io/img/Lena.png"))
dim(photo)
n <- nrow(photo)
p <- ncol(photo)

#SVD sur l'image

photo.svd <- svd(photo)
# UDV
u <- photo.svd$u
d <- photo.svd$d
v <- photo.svd$v

start.time <- Sys.time()
par(mar=c(0,0,2,0),mfrow=c(2,3),xaxt="n",yaxt="n")
image(photo, col = grey(seq(0, 1, length = 256)),asp=1)
title(paste("Original, 100%"),cex.main=.9)
r=1 ; image((u[,1]*d[1])%*%t(v[,1]), col = grey(seq(0, 1, length = 256)),asp=1)
title(paste("r=",r, ", ",round(((n+1+p)*r)/(n*p)*100,1),"%"),cex.main=0.9)
for (r in c(10, 20, 50, 100)){
  image(u[,1:r]%*%diag(d[1:r])%*%t(v[,1:r]), col = grey(seq(0, 1, length = 256)),asp=1)
  title(paste("r=",r, ", ",round(((n+1+p)*r)/(n*p)*100,1),"%"),cex.main=0.9)
}
Sys.time()- start.time


#RSVD sur l'image

k = 512
photo.rsvd <- rsvd(photo, k=k)

Phto <- photo.rsvd$u %*% diag(photo.rsvd$d) %*% t(photo.rsvd$v)

print(100 * norm( photo - Phto, 'F') / norm( photo,'F')) # percentage error

start.time <- Sys.time()
par(mar=c(0,0,2,0),mfrow=c(2,3),xaxt="n",yaxt="n")
image(photo, col = grey(seq(0, 1, length = 256)),asp=1)
title(paste("Original, 100%"),cex.main=.9)
r=1 ; image((photo.rsvd$u[,1]*photo.rsvd$d[1])%*%t(photo.rsvd$v[,1]), col = grey(seq(0, 1, length = 256)),asp=1)
title(paste("r=",r, ", ",round(((n+1+p)*r)/(n*p)*100,1),"%"),cex.main=0.9)
for (r in c(10, 20, 50, 100)){
  image(photo.rsvd$u[,1:r]%*%diag(photo.rsvd$d[1:r])%*%t(photo.rsvd$v[,1:r]), col = grey(seq(0, 1, length = 256)),asp=1)
  title(paste("r=",r, ", ",round(((n+1+p)*r)/(n*p)*100,1),"%"),cex.main=0.9)
}
Sys.time()- start.time
