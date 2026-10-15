SLsn<- readRDS("C:\\Users\\heeju\\OneDrive\\바탕 화면\\RCH\\SkewHeckman\\snsn.rds")
SLn<- readRDS("C:\\Users\\heeju\\OneDrive\\바탕 화면\\RCH\\SkewHeckman\\nnnn.rds")
SLsn[2]

library(rstan)
library(bayesplot)
library(ggplot2)
fit <- as.array(SLsn[[1]])

mcmc_trace(fit, pars=c("beta[1]","beta[2]","beta[3]","gamma[1]","gamma[2]","gamma[3]","gamma[4]","gamma[5]"
                       ,"gamma[6]","gamma[7]","rho","sigma","lambda" ))
           
##

# fit: centered 모델로 돌린 stanfit 객체
pars <- c(paste0("beta[", 1:3, "]"),
          paste0("gamma[", 1:7, "]"),
          "rho", "sigma", "alpha")

# R-hat을 패널 제목에 붙이기
rh  <- summary(SLsn[[1]], pars = pars)$summary[, "Rhat"]
lab <- setNames(sprintf("%s  (R-hat = %.3f)", pars, rh), pars)

# 체인을 합쳐서 한 색으로 (그림처럼)
draws <- as.matrix(SLsn[[1]], pars = pars)

dimnames(draws)$parameters[
  dimnames(draws)$parameters == "lambda"
] <- "alpha"


color_scheme_set("red")
p <- mcmc_trace(draws, pars = pars,
                facet_args = list(ncol = 3, labeller = as_labeller(lab))) +
  xlab("iteration") +
  theme_classic(base_family = "serif") +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(size = 8))

ggsave("app1_convergence.pdf", p, width = 8, height = 11)
ggsave("app1_convergence.eps", p, width = 8, height = 11)

SLsn[[1]]
