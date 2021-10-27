library(rsvd)

?rsvd #On regarde tout d'abord ce que fait la fonction rsvd

#####################################
#####   Comparaison "à froid"   #####
#####################################

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

#Séparation de l'onglet "Plot" pour obtenir 6 graphiques et les comparer
par(mar=c(0,0,2,0),mfrow=c(2,3),xaxt="n",yaxt="n")

#Image originale en noir et blanc
image(photo, col = grey(seq(0, 1, length = 256)),asp=1)
title(paste("Original, 100%"),cex.main=.9)

#Image conservant uniquement le rang r égal 1
r=1 ; image((u[,1]*d[1])%*%t(v[,1]), col = grey(seq(0, 1, length = 256)),asp=1)
title(paste("r=",r, ", ",round(((n+1+p)*r)/(n*p)*100,1),"%"),cex.main=0.9)

### Pour r égal à 10 puis 20 puis 50 puis 100 ###
### Représentation graphique ###

for (r in c(10, 20, 50, 100)){
  image(u[,1:r]%*%diag(d[1:r])%*%t(v[,1:r]), col = grey(seq(0, 1, length = 256)),asp=1)
  title(paste("r=",r),cex.main=0.9)
}

# r = 50
# image(u[,1:r]%*%diag(d[1:r])%*%t(v[,1:r]), col = grey(seq(0, 1, length = 256)),asp=1)
# title(paste("r=",r),cex.main=0.9)



#RSVD sur l'image

k = 50
photo.rsvd <- rsvd(photo, k=k)

Phto <- photo.rsvd$u %*% diag(photo.rsvd$d) %*% t(photo.rsvd$v)

print(100 * norm( photo - Phto, 'F') / norm( photo,'F')) # percentage error

#Séparation de l'onglet "Plot" pour obtenir 6 graphiques et les comparer
par(mar=c(0,0,2,0),mfrow=c(2,3),xaxt="n",yaxt="n") 

#Image originale en noir et blanc
image(photo, col = grey(seq(0, 1, length = 256)),asp=1)
title(paste("Original, 100%"),cex.main=.9)

#Image conservant uniquement le rang k égal 1
k=1 ; image((photo.rsvd$u[,1]*photo.rsvd$d[1])%*%t(photo.rsvd$v[,1]), col = grey(seq(0, 1, length = 256)),asp=1)
title(paste("k=",k, ", ",round(((n+1+p)*k)/(n*p)*100,1),"%"),cex.main=0.9)

### Pour k égal à 10 puis 20 puis 50 puis 100 ###
### Représentation graphique ###

for (k in c(10, 20, 50, 100)){
  image(photo.rsvd$u[,1:k]%*%diag(photo.rsvd$d[1:k])%*%t(photo.rsvd$v[,1:k]), col = grey(seq(0, 1, length = 256)),asp=1)
  title(paste("k=",k),cex.main=0.9)
}

# par(xaxt="n",yaxt="n")
# k = 50
# image(photo.rsvd$u[,1:k]%*%diag(photo.rsvd$d[1:k])%*%t(photo.rsvd$v[,1:k]), col = grey(seq(0, 1, length = 256)),asp=1)
# title(paste("k=",k),cex.main=0.9)


###################################
##### Test du computing time  #####
###################################

# Importation d'une image

photo <- as.matrix(raster("https://husson.github.io/img/Lena.png"))
dim(photo)
n <- nrow(photo)
p <- ncol(photo)

#SVD sur l'image

start.time <- Sys.time()
photo.svd <- svd(photo)
compsvd <- Sys.time() - start.time
compsvd

#RSVD sur l'image
k = 50

start.time <- Sys.time()
photo.rsvd <- rsvd(photo, k=k)
comprsvd <- Sys.time() - start.time
comprsvd

test1 <- rsvd(photo)
test <- rsvd(photo, p=0, q=1)
# computetime <- rep(0,ncol(photo))

# for (i in 1:ncol(photo)){
#   start.time <- Sys.time()
#   photo.rsvd <- rsvd(photo, k=k)
#   comprsvd <- Sys.time() - start.time
#   computetime[i] <- comprsvd
# }

# write.table(computetime,"computetime.txt")
computetime <- read.table("computetime.txt", skip = 1, col.names = c("rank","temps"))
length(computetime)

head(computetime)
tail(computetime)

#Peu importe le rang indiqué pour effectuer la rsvd, le computing time reste trois fois inférieur à celui de la SVD


#############################################
##### oversampling et power iterations  #####
#############################################

#### oversampling ####

par(mar=c(0,0,2,0),mfrow=c(2,3),xaxt="n",yaxt="n")

k = 10
photo.rsvd <- rsvd(photo, k=k)
image(photo.rsvd$u[,1:k]%*%diag(photo.rsvd$d[1:k])%*%t(photo.rsvd$v[,1:k]), col = grey(seq(0, 1, length = 256)),asp=1)

Photo.recons.rsvdP <- list()

for (p in c(10,20,50,100)){
  photo.rsvdP <- rsvd(photo, k=k, p = p)
  #image(photo.rsvdP$u%*%diag(photo.rsvdP$d)%*%t(photo.rsvdP$v), col = grey(seq(0, 1, length = 256)),asp=1)
  if (p == 10){
    Photo.recons.rsvdP[[1]] <- photo.rsvdP$u %*% diag(photo.rsvdP$d) %*% t(photo.rsvdP$v)
  }
  else if (p == 20){
    Photo.recons.rsvdP[[2]] <- photo.rsvdP$u %*% diag(photo.rsvdP$d) %*% t(photo.rsvdP$v)
  }
  else if (p == 50){
    Photo.recons.rsvdP[[3]] <- photo.rsvdP$u %*% diag(photo.rsvdP$d) %*% t(photo.rsvdP$v)
  }
  else {
    Photo.recons.rsvdP[[4]] <- photo.rsvdP$u %*% diag(photo.rsvdP$d) %*% t(photo.rsvdP$v)
  }
}

Photo.rsvdP <- list()
for (p in c(10,20,50,100)){
  if (p == 10){
    Photo.rsvdP[[1]] <- rsvd(photo, k=k, p = p)
  }
  else if (p == 20){
    Photo.rsvdP[[2]] <- rsvd(photo, k=k, p = p)
  }
  else if (p == 50){
    Photo.rsvdP[[3]] <- rsvd(photo, k=k, p = p)
  }
  else {
    Photo.rsvdP[[4]] <- rsvd(photo, k=k, p = p)
  }
}

table(Photo.rsvdP[[2]]$d==Photo.rsvdP[[1]]$d) #preuve d'une différence entre les matrices


table(Photo.recons.rsvdP[[2]]==Photo.recons.rsvdP[[1]]) #preuve d'une différence entre les matrices

Photo.recons.rsvd <- photo.rsvd$u %*% diag(photo.rsvd$d) %*% t(photo.rsvd$v)

#Les deux matrices sont bel et bien différentes, 
#même si aucune différence n'est visible sur l'image

for (i in 1:4){
  print(table(Photo.recons.rsvd==Photo.recons.rsvdP[[i]]))
}


#### power iterations ####

#La même chose avec le power iteration, on évite la partie graphique pour maximiser la rapidité

Photo.recons.rsvdQ <- list()

photo.rsvdQ5 <- rsvd(photo, k=10, q = 5)
photo.rsvdQ0 <- rsvd(photo, k=10, q = 1)

res5 <- photo.rsvdQ5$d/max(photo.rsvdQ5$d)
res0 <- photo.rsvdQ0$d/max(photo.rsvdQ0$d)

plot(res5, type = "l", col = "red",
     xlab = "rang", ylab = expression(sigma), main = "Valeurs singulières en fonction du rang")
lines(res0, type = "l", col = "blue")
legend(8, 0.4, c("q = 5", "q = 1"), fill = c("red", "blue"))


# 
# plot(res_10_5, type = "l", col = "red")
# lines(res_10_0, type = "l", col = "blue")


for (q in c(3,4,5,10)){
  photo.rsvdQ <- rsvd(photo, k=k, q = q)
  if (q == 3){
    Photo.recons.rsvdQ[[1]] <- photo.rsvdQ$u %*% diag(photo.rsvdQ$d) %*% t(photo.rsvdQ$v)
  }
  else if (q == 4){
    Photo.recons.rsvdQ[[2]] <- photo.rsvdQ$u %*% diag(photo.rsvdQ$d) %*% t(photo.rsvdQ$v)
  }
  else if (q == 5){
    Photo.recons.rsvdQ[[3]] <- photo.rsvdQ$u %*% diag(photo.rsvdQ$d) %*% t(photo.rsvdQ$v)
  }
  else {
    Photo.recons.rsvdQ[[4]] <- photo.rsvdQ$u %*% diag(photo.rsvdQ$d) %*% t(photo.rsvdQ$v)
  }
}

print(table(Photo.recons.rsvdQ[[2]]==Photo.recons.rsvdQ[[1]])) #preuve d'une différence entre les matrices

Photo.recons.rsvd <- photo.rsvd$u %*% diag(photo.rsvd$d) %*% t(photo.rsvd$v)

for (i in 1:4){
  print(table(Photo.recons.rsvd==Photo.recons.rsvdQ[[i]]))
}

### Test des valeurs d'erreurs ### (éloignement de la vérité)

#Oversampling

for (i in 1:4){
  print(100 * norm( photo - Photo.recons.rsvdP[[i]], 'F') / norm( photo,'F')) # percentage error
}

#Power iteration

for (i in 1:4){
  print(100 * norm( photo - Photo.recons.rsvdQ[[i]], 'F') / norm( photo,'F')) # percentage error
}

#Base RSVD

print(100 * norm( photo - Photo.recons.rsvd, 'F') / norm( photo,'F')) # percentage error

#L'augmentation du score dépend des colonnes aléatoires prises avec l'oversampling et à la valeur attribuée lors de la power iteration


#Base SVD

Photo.recons.svd <- u[,1:k]%*%diag(d[1:k])%*%t(v[,1:k])

print(100 * norm( photo - Photo.recons.svd, 'F') / norm( photo,'F')) # percentage error

table(Photo.recons.svd==Photo.recons.rsvd)

#Le pourcentage d'erreur est similaire en svd et en rsvd pour des mêmes paramètres



###################################
##### Test avec un jeu généré #####
###################################

#Simulate a general matrix with 1000 rows and 1000 columns
vy= rnorm(1000*1000,0,1)
y= matrix(vy,1000,1000,byrow=TRUE)

#Compute the randSVD for the first hundred components of the matrix y and measure the time
start.time <- Sys.time()
prova1= rsvd(y,k=100)
Sys.time()- start.time

start.time <- Sys.time()
prova3= rsvd(y,k=100,q=5,p=20)
Sys.time()- start.time

#Compare with a classical SVD
start.time <- Sys.time()
prova2= svd(y)
Sys.time()- start.time

P1 <- prova1$u %*% diag(prova1$d) %*% t(prova1$v)

print(100 * norm( y - P1, 'F') / norm( y,'F')) # percentage error

P3 <- prova3$u %*% diag(prova3$d) %*% t(prova3$v)

print(100 * norm( y - P3, 'F') / norm( y,'F')) # percentage error

P2 <- prova2$u[,1:k] %*% diag(prova2$d[1:k]) %*% t(prova2$v[,1:k])

print(100 * norm( y - P2, 'F') / norm( y,'F')) # percentage error


rsvdP <- list()
for (p in c(10,20,50,100)){
  if (p == 10){
    rsvdP[[1]] <- rsvd(y, k=k, p = p)
  }
  else if (p == 20){
    rsvdP[[2]] <- rsvd(y, k=k, p = p)
  }
  else if (p == 50){
    rsvdP[[3]] <- rsvd(y, k=k, p = p)
  }
  else {
    rsvdP[[4]] <- rsvd(y, k=k, p = p)
  }
}
for (i in 1:4){
  cat(c("5 premières valeurs singulières pour rsvd",i, ": ",round(head(rsvdP[[i]]$d, 5),2),"\n"))
}

table(rsvdP[[2]]$d==rsvdP[[1]]$d) #preuve d'une différence entre les matrices



rsvdQ5 <- rsvd(y, k=50, q = 5)
rsvdQ0 <- rsvd(y, k=50, q = 1)

res5 <- rsvdQ5$d/max(rsvdQ5$d)
res0 <- rsvdQ0$d/max(rsvdQ0$d)

plot(res5, type = "l", col = "red", ylim = c(0, 1),
     xlab = "rang", ylab = expression(sigma), main = "Valeurs singulières en fonction du rang")
lines(res0, type = "l", col = "blue")
legend(40, 0.4, c("q = 5", "q = 1"), fill = c("red", "blue"))



#################
### Curiosité ###
#################

photo <- as.matrix(raster("https://husson.github.io/img/Lena.png"))

#Et si on prend les r dernières valeurs ?

#SVD

for (r in c(10, 20, 50, 100)){
  image(u[,(512-r):512]%*%diag(d[(512-r):512])%*%t(v[,(512-r):512]), col = grey(seq(0, 1, length = 256)),asp=1)
  title(paste("r=",r, ", ",round(((n+1+p)*r)/(n*p)*100,1),"%"),cex.main=0.9)
}


#RSVD

par(new=TRUE)
for (k in c(10, 20, 50, 100)){
  image(photo.rsvd$u[,(512-k):512]%*%diag(photo.rsvd$d[(512-k):512])%*%t(photo.rsvd$v[,(512-k):512]), col = grey(seq(0, 1, length = 256)),asp=1)
  title(paste("k=",k, ", ",round(((n+1+p)*k)/(n*p)*100,1),"%"),cex.main=0.9)
}


par(new=TRUE)

photo.rsvd <- rsvd(photo, k = 200)
dern <- 20

image(photo.rsvd$u[,(k-dern):k]%*%diag(photo.rsvd$d[(k-dern):k])%*%t(photo.rsvd$v[,(k-dern):k]), col = grey(seq(0, 1, length = 256)),asp=1)

photo.rsvdP <- rsvd(photo, k = 200, p = 30)
image(photo.rsvdP$u[,(k-dern):k]%*%diag(photo.rsvdP$d[(k-dern):k])%*%t(photo.rsvdP$v[,(k-dern):k]), col = grey(seq(0, 1, length = 256)),asp=1)
