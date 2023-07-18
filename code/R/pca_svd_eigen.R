princomp(USArrests, cor = TRUE)$loadings

prcomp(USArrests, scale = TRUE)
prcomp(~ Murder + Assault + Rape, data = USArrests, scale = TRUE)
plot(prcomp(USArrests))
summary(prcomp(USArrests, scale = TRUE))
biplot(prcomp(USArrests, scale = TRUE))

svd(cor(USArrests))$u
eigen(cor(USArrests))$vectors
