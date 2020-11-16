library(kinship2)
library(MCMCglmm)
library(pedigree)
library(sqldf)
library(coda)

data("minnbreast")

subset <- sqldf("SELECT * FROM minnbreast WHERE bcpc = '1'")
fam <- sqldf("SELECT distinct famid FROM subset GROUP BY famid")
minnbreast <- sqldf("SELECT a.* FROM minnbreast as a INNER JOIN fam as b WHERE a.famid=b.famid")

# Ordering the pedigree data

Ped <- minnbreast[,c(1,3,4)]
(ord <- orderPed(Ped))
pedigree <- Ped[nrow(minnbreast):1,]
(ord <- orderPed(Ped))
pedigree <- Ped[order(ord),]
pedigree <- as.data.frame(pedigree)


# Renaming variables

colnames(minnbreast)[which(names(minnbreast) == "id")] <- "animal"
colnames(pedigree)[which(names(pedigree) == "id")] <- "ID"
colnames(pedigree)[which(names(pedigree) == "fatherid")] <- "FATHER"
colnames(pedigree)[which(names(pedigree) == "motherid")] <- "MOTHER"

pedigree$FATHER[pedigree$FATHER == 0] <- NA
pedigree$MOTHER[pedigree$MOTHER == 0] <- NA


# Bivariate model preparation:

# split male and female


minnbreast$prostate <- minnbreast$sex == "M" & minnbreast$cancer == 1
minnbreast$breast <- minnbreast$sex == "F" & minnbreast$cancer == 1

minnbreast$prostate[minnbreast$prostate == TRUE] <- 1
minnbreast$prostate[minnbreast$prostate == 0 & minnbreast$sex=="F"] <- NA

minnbreast$breast[minnbreast$breast == TRUE] <- 1
minnbreast$breast[minnbreast$breast == 0 & minnbreast$sex=="M"] <- NA

prior.bivariate <- list(G = list(G1 = list(V = diag(2), n = 0.02)),R = list(V = diag(2), n = 0.02))

model.bivariate <-
  MCMCglmm(
    cbind(prostate,breast) ~ trait -1,
    random = ~ us(trait):animal,
    rcov = ~us(trait):units,
    family = rep("gaussian",2),
    data = minnbreast,
    pedigree = pedigree,
    prior = prior.bivariate,
    nitt= 600000,
    thin = 1000
  )


summary(model.bivariate)
plot(model.bivariate$Sol)
plot(model.bivariate$VCV)

herit1<-model.bivariate$VCV[,'traitprostate:traitprostate.animal']/
  (model.bivariate$VCV[,'traitprostate:traitprostate.animal']+model.bivariate$VCV[,'traitprostate:traitprostate.units'])
herit2<-model.bivariate$VCV[,'traitbreast:traitbreast.animal']/
  (model.bivariate$VCV[,'traitbreast:traitbreast.animal']+model.bivariate$VCV[,'traitbreast:traitbreast.units'])
mean(herit1)
mean(herit2)

autocorr(model.bivariate$Sol)
# Correlation for the variance:
autocorr(model.bivariate$VCV)

# Liability Threshold Model
h2_o_prostate <- mean(herit1)
K_prostate <- mean(model.bivariate$Sol[,1])
z_prostate <- dnorm(qnorm(K_prostate))
h2_l_prostate <- h2_o_prostate * K_prostate * (1-K_prostate)/(z_prostate^2)
h2_l_prostate

h2_o_breast <- mean(herit2)
K_breast <- mean(model.bivariate$Sol[,2])
z_breast <- dnorm(qnorm(K_breast))
h2_l_breast <- h2_o_breast * K_breast * (1-K_breast)/(z_breast^2)
h2_l_breast
