.libPaths("~/rlibs_R461")
if(!require(HeckmanEM))   install.packages("HeckmanEM");   library(HeckmanEM)
if(!require(mvtnorm))   install.packages("mvtnorm");   library(mvtnorm)
if(!require(loo))   install.packages("loo"); library("loo")
if(!require(rstan))   install.packages("rstan"); library(rstan, quietly = T)
if(!require(shinystan))   install.packages("shinystan"); library(shinystan)

source("UtilitariesHeckmanNew.R")

#data1 <- read.csv("statadata2.csv")

#x<- cbind(1, data1$age, data1$female, data1$educ, data1$blhisp,  data1$totchr, data1$ins)
#w<- cbind(x, data1$income)
#data1$lambexp[is.na(data1$lambexp)]<-0
#cc<- data1$dambexp                                                                                                            
#y<- data1$lnambx
#n<-length(cc)


library(AER)
data("PSID1976")

cc         <- ifelse(PSID1976$participation=="yes", 1, 0) # cc == 0  missing data
wage      <- PSID1976$wage
hwage      <- PSID1976$hwage
youngkids  <- PSID1976$youngkids
tax        <- PSID1976$tax
feducation <- PSID1976$feducation
education  <- PSID1976$education
city       <- PSID1976$city

# Transformado as variáveis "factors" em quantitativas
city     <- ifelse(city=="yes", 1, 0)     # yes == 1 and no == 0

x1= hwage
x2= youngkids
x3= tax
x4= feducation
x5= education
x6= city

#library(HeckmanEM)

x <- cbind(1,x5,x6)
w <- cbind(1,x1,x2,x3,x4,x5,x6)
y <- ifelse(wage>0, log(wage), 0)


n<-length(cc)


data = list(N = n, N_y = sum(cc==1), p = ncol(x), q = ncol(w), X = x[cc > 0, ], Z = w, D = cc, y = y[cc > 0])

fit.n_stan <- stan(
  file = "HeckmanNormal.stan",
  data = data,
  thin = 5,
  chains = 1,
  iter = 10000,
  warmup = 1000
)



fit.sn_stan <- stan(
  file = "HeckmanSkewNormallast.stan",
  data = data,
  thin = 5,
  chains = 1,
  iter = 10000,
  warmup = 1000
)

# --------------------------------------------------
# Results
# --------------------------------------------------

fit.sn_stan

print(
  fit.sn_stan,
  pars = c(
    "beta",
    "gamma",
    "sigma2",
    "rho",
    "lambda"
  )
)
