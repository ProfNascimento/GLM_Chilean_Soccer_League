rm(list = ls())
require(rstan)
require(matrixStats)
set.seed(1) #set seed
options(mc.cores=3)

data = read.csv("Ligachie2020.csv",
                col.names = c("Home","score1","score2", "Away"), 
                stringsAsFactors = FALSE, sep = ";",header = FALSE)
head(data)

ng = nrow(data) # Number of games
nt = length(unique(data$Home)) # Number of teams

#Convert team names for each match into numbers
teams = unique(data$Away)
ht = unlist(sapply(1:ng, function(g) which(teams == data$Home[g])))
at = unlist(sapply(1:ng, function(g) which(teams == data$Away[g])))

np=5 #Save the last 5 games to predict
ngob = ng-np #number of games to fit

# Data Wrangling 
my_data = list(
  nt = nt,
  ng = ngob ,
  ht = ht[1:ngob],
  at = at[1:ngob],
  s1 = data$score1[1:ngob],
  s2 = data$score2[1:ngob],
  np = np,
  htnew = ht[(ngob+1):ng],
  atnew = at[(ngob+1):ng]
)

#################################################
##      MODEL 1 - Non-Hierarchical model       ##
## Traditional method (Non-informative Priors) ##
#################################################
nhpoolfit = stan(file = "~/STAN/non_hier_model.stan", 
                 data = my_data, iter = 2000, chains = 4, seed = 123456)
print(nhpoolfit)

#Plot the predicted scores of the last 5 matches
nhpoolparams = extract(nhpoolfit)
pred_scores = c(colMedians(nhpoolparams$s1new),colMedians(nhpoolparams$s2new))
pred_errors = c(sapply(1:np, function(x) sd(nhpoolparams$s1new[,x])),
                sapply(1:np, function(x) sd(nhpoolparams$s1new[,x])))
true_scores = c(data$score1[(ngob+1):ng],data$score2[(ngob+1):ng] )
plot(true_scores , pred_scores , xlim=c(0,4), ylim=c(0,2), pch=20,
     ylab="predicted scores", xlab="true scores")
abline(a=0, b=1, lty="dashed")
arrows(true_scores , pred_scores+pred_errors , true_scores ,
       pred_scores -pred_errors , length = 0.05 , angle = 90, code = 3)

#Teams " attack and defense abilities are all close to 0.
attack = colMedians(nhpoolparams$att)
defense = colMedians(nhpoolparams$def)

plot(attack ,defense ,xlim=c( -1.7 ,1),ylim=c( -1,1.5))
abline(h=0)
abline(v=0)
text(attack ,defense , labels=teams , cex=0.7 , pos=4)

##################################
## MODEL 2 - Hierarchical model ##
##################################
hfit = stan(file = "~/STAN/hier_model.stan", 
            data = my_data, iter = 2000, chains = 4, seed = 123456)

print(hfit)
pairs(hfit , pars=c("mu_att", "tau_att", "mu_def", "tau_def"))

hparams = extract(hfit)
pred_scoresH = c(colMedians(hparams$s1new),colMedians(hparams$s2new))
pred_errorsH = c(sapply(1:np, function(x) sd(hparams$s1new[,x])),
                sapply(1:np, function(x) sd(hparams$s1new[,x])))
true_scores = c(data$score1[(ngob+1):ng],data$score2[(ngob+1):ng] )
plot(true_scores , pred_scoresH , xlim=c(0,4), ylim=c(0,2), pch=20,
     ylab="predicted scores", xlab="true scores")
abline(a=0, b=1, lty="dashed")
arrows(true_scores , pred_scores+pred_errors , true_scores ,
       pred_scores -pred_errors , length = 0.05 , angle = 90, code = 3, rgb(0 ,0 ,0 ,0.3))

attack = colMedians(hparams$att)
attacksd = sapply(1:nt, function(x) sd(hparams$att[,x]))
defense = colMedians(hparams$def)
defensesd = sapply(1:nt, function(x) sd(hparams$def[,x]))

plot(attack ,defense , xlim=c(-0.1,0.7), ylim=c( -0.01 ,0.005) , pch=20)
arrows(attack -attacksd , defense , attack+attacksd , defense , code=3,
       angle = 90, length = 0.04 , col=rgb(0 ,0 ,0 ,0.2))
arrows(attack , defense -defensesd , attack , defense+defensesd , code=3,
       angle = 90, length = 0.04 ,col=rgb(0 ,0 ,0 ,0.2))
abline(h=0)
abline(v=0)
text(attack ,defense , labels=teams , cex=0.7 , adj=c( -0.05 , -0.8) )
points(x=attack[2],y=defense[2],col="red",pch=16)
points(x=attack[4],y=defense[4],col="red",pch=16)
points(x=attack[5],y=defense[5],col="red",pch=16)
points(x=attack[6],y=defense[6],col="red",pch=16)
points(x=attack[9],y=defense[9],col="red",pch=16)
points(x=attack[10],y=defense[10],col="red",pch=16)
points(x=attack[11],y=defense[11],col="red",pch=16)
points(x=attack[12],y=defense[12],col="red",pch=16)
points(x=attack[13],y=defense[13],col="red",pch=16)
points(x=attack[15],y=defense[15],col="red",pch=16)
points(x=attack[17],y=defense[17],col="red",pch=16)
points(x=attack[18],y=defense[18],col="red",pch=16)

# Checking the Computational Inferences (CHAINS)
require(shinystan)
summary(hfit)
teams
launch_shinystan(hfit)
a<-summary(hfit)
tabla<-data.frame(a)
write.csv2(tabla, "tabla3.csv")

######################
## MODEL COMPARISON ##
######################
library("loo")

log_lik <- extract_log_lik(nhpoolfit, merge_chains = FALSE)
r_eff <- relative_eff(exp(log_lik))
loo <- loo(log_lik, r_eff = r_eff, cores = 2)
print(loo)

log_lik_h <- extract_log_lik(hfit, merge_chains = FALSE)
r_eff_h <- relative_eff(exp(log_lik_h))
loo_h <- loo(log_lik_h, r_eff = r_eff_h, cores = 2)
print(loo_h)

# LOO Model Comparison
loo_compare(x = list(loo, loo_h))