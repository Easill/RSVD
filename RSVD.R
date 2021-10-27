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
  title(paste("r=",r, ", ",round(((n+1+p)*r)/(n*p)*100,1),"%"),cex.main=0.9)
}


#RSVD sur l'image

k = 512
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
  title(paste("k=",k, ", ",round(((n+1+p)*k)/(n*p)*100,1),"%"),cex.main=0.9)
}

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

table(Photo.recons.rsvdP[[2]]==Photo.recons.rsvdP[[1]]) #preuve d'une différence entre les matrices

Photo.recons.rsvd <- photo.rsvd$u %*% diag(photo.rsvd$d) %*% t(photo.rsvd$v)

for (i in 1:4){
  print(table(Photo.recons.rsvd==Photo.recons.rsvdP[[i]]))
}
#Les deux matrices sont bel et bien différentes, 
#même si aucune différence n'est visible sur l'image

#### power iterations ####

#La même chose avec le power iteration, on évite la partie graphique pour maximiser la rapidité

Photo.recons.rsvdQ <- list()

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

k = 100

start.time <- Sys.time()
photo.rsvd <- rsvd(photo, k=k)
comprsvd <- Sys.time() - start.time
comprsvd


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
prova3= rsvd(y,k=100,q=5,p=50)
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
