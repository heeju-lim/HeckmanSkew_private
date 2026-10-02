

load("C:\\Users\\heeju\\OneDrive\\바탕 화면\\RCH\\Stan\\application\\fit.n_stan.RData")

fit.n_stan[[2]]

fit.sn_stan <- readRDS("C:\\Users\\heeju\\OneDrive\\바탕 화면\\Bayesian Heckman codes\\app2_fit_sn_stan.rds")
print(
  fit.sn_stan,
  pars = c(
    "beta",
    "gamma",
    "sigma2",
    "rho",
    "lambda",
    "log_lik_total"
  )
)
str(fit.n_stan)

library(rstan)
pdf("fit.n_stan2.pdf")
plot(fit.n_stan[[1]],pars=c("beta[2]","beta[3]","beta[4]","beta[5]","beta[6]","beta[7]", "gamma[2]", "gamma[3]","gamma[4]","gamma[5]","gamma[6]","gamma[7]","gamma[8]"))+ggtitle("SLn")+ theme(plot.title = element_text(hjust = 0.5))
dev.off()


library(rstan)
pdf("fit.n_stan2.pdf")
plot(fit.sn_stan,pars=c("beta[2]","beta[3]","beta[4]","beta[5]","beta[6]","beta[7]", "gamma[2]", "gamma[3]","gamma[4]","gamma[5]","gamma[6]","gamma[7]","gamma[8]"))+ggtitle("SLn")+ theme(plot.title = element_text(hjust = 0.5))
dev.off()



library(ggplot2)

pars <- c(paste0("beta[", 2:7, "]"), paste0("gamma[", 2:8, "]"))

get_ci <- function(fit, model) {
  draws <- as.matrix(fit, pars = pars)
  data.frame(
    param = pars,
    model = model,
    med   = apply(draws, 2, median),
    l95   = apply(draws, 2, quantile, 0.025),
    u95   = apply(draws, 2, quantile, 0.975),
    l80   = apply(draws, 2, quantile, 0.10),
    u80   = apply(draws, 2, quantile, 0.90)
  )
}

df <- rbind(get_ci(fit.n_stan[[1]], "SLn"),
            get_ci(fit.sn_stan,      "SLsn"))

df$param <- factor(df$param, levels = rev(pars))          # 위에서부터 beta[2] ...
df$model <- factor(df$model, levels = c("SLsn", "SLn"))   # 각 파라미터에서 SLn이 위쪽

pd <- position_dodge(width = 0.6)

ggplot(df, aes(x = param, y = med, colour = model)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_linerange(aes(ymin = l95, ymax = u95), position = pd, linewidth = 0.5) +
  geom_linerange(aes(ymin = l80, ymax = u80), position = pd, linewidth = 1.6) +
  geom_point(position = pd, size = 2.2) +
  coord_flip() +
  scale_colour_manual(values = c(SLn = "#1f4e79", SLsn = "#b2182b"),
                      guide = guide_legend(reverse = TRUE)) +
  labs(x = NULL, y = "Posterior estimate", colour = NULL) +
  theme_bw() +
  theme(legend.position = "top",
        axis.text.y = element_text(face = "bold")) +
  scale_x_discrete(labels = c(
    "beta[2]" = "age (O)", "beta[3]" = "female (O)", "beta[4]" = "educ (O)",
    "beta[5]" = "blhisp (O)", "beta[6]" = "totchr (O)", "beta[7]" = "ins (O)",
    "gamma[2]" = "age (S)", "gamma[3]" = "female (S)", "gamma[4]" = "educ (S)",
    "gamma[5]" = "blhisp (S)", "gamma[6]" = "totchr (S)", "gamma[7]" = "ins (S)",
    "gamma[8]" = "income (S)")) 



library(ggplot2)

pars <- c(paste0("beta[", 2:7, "]"), paste0("gamma[", 2:8, "]"))
vars_out <- c("age", "female", "educ", "blhisp", "totchr", "ins")
vars_sel <- c("age", "female", "educ", "blhisp", "totchr", "ins", "income")

get_ci <- function(fit, model) {
  draws <- as.matrix(fit, pars = pars)
  data.frame(
    param = pars,
    model = model,
    med = apply(draws, 2, median),
    l95 = apply(draws, 2, quantile, 0.025),
    u95 = apply(draws, 2, quantile, 0.975),
    l80 = apply(draws, 2, quantile, 0.10),
    u80 = apply(draws, 2, quantile, 0.90)
  )
}

df <- rbind(get_ci(fit.n_stan[[1]], "SLn"),
            get_ci(fit.sn_stan,      "SLsn"))

# 모형 구분 + 변수 이름 + 세로 위치
df$eq  <- factor(ifelse(grepl("beta", df$param), "Outcome", "Selection"),
                 levels = c("Outcome", "Selection"))
key <- data.frame(
  param = pars,
  var   = c(vars_out, vars_sel),
  pos   = c(14:9, 7:1)          # 위에서부터 순서대로
)
df <- merge(df, key, by = "param")
df$y <- df$pos + ifelse(df$model == "SLn", 0.18, -0.18)   # SLn 위, SLsn 아래

ggplot(df, aes(y = y, colour = model)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_segment(aes(x = l95, xend = u95, yend = y), linewidth = 0.5) +
  geom_segment(aes(x = l80, xend = u80, yend = y), linewidth = 1.6) +
  geom_point(aes(x = med), size = 2.2) +
  scale_y_continuous(breaks = key$pos, labels = key$var,
                     expand = expansion(add = 0.6)) +
  scale_colour_manual(values = c(SLn = "#1f4e79", SLsn = "#b2182b")) +
  facet_grid(eq ~ ., scales = "free_y", space = "free_y", switch = "y") +
  labs(x = "Posterior estimate", y = NULL, colour = NULL) +
  theme_bw() +
  theme(legend.position  = "top",
        strip.placement  = "outside",
        strip.background = element_rect(fill = "grey90"),
        strip.text.y.left = element_text(angle = 90, face = "bold"),
        axis.text.y      = element_text(face = "bold"),
        panel.grid.minor = element_blank()) +
  theme_bw(base_size = 15) +
  theme(legend.position   = "top",
        legend.text       = element_text(size = 15),
        strip.placement   = "outside",
        strip.background  = element_rect(fill = "grey90"),
        strip.text.y.left = element_text(angle = 90, face = "bold", size = 15),
        axis.text.y       = element_text(face = "bold", size = 14, colour = "black"),
        axis.text.x       = element_text(size = 13, colour = "black"),
        axis.title.x      = element_text(size = 15),
        panel.grid.minor  = element_blank())

ggsave("app2_compare.pdf", width = 7, height = 6.5)
