## ------------------------------------------------------------------
## (a) Posterior estimates SLn vs SLsn   (b) Q-Q plot of quantile residuals
## Both panels share colours, shapes, legend and font sizes.
## ------------------------------------------------------------------
library(rstan)
library(ggplot2)
library(patchwork)
library(pbivnorm)

load("C:\\Users\\heeju\\OneDrive\\바탕 화면\\RCH\\Stan\\application\\fit.n_stan2.RData")
fit.sn_stan <- readRDS("C:\\Users\\heeju\\OneDrive\\바탕 화면\\Bayesian Heckman codes\\app1_fit_sn_stan.rds")

fit_n  <- SLn[[1]]
fit_sn <- SLsn[[1]]

## ---- common style (SLn: blue circle, SLsn: red triangle) ----
cols   <- c(SLn = "#1f4e79", SLsn = "#b2182b")
shapes <- c(SLn = 16, SLsn = 17)

base_theme <- theme_bw(base_size = 15) +
  theme(legend.position  = "top",
        legend.text      = element_text(size = 15),
        axis.text.y      = element_text(face = "bold", size = 14, colour = "black"),
        axis.text.x      = element_text(size = 13, colour = "black"),
        axis.title.x     = element_text(size = 15),
        axis.title.y     = element_text(size = 15),
        panel.grid.minor = element_blank())

## =================================================================
## (a) Posterior estimates
## =================================================================
pars     <- c(paste0("beta[", 2:3, "]"), paste0("gamma[", c(2:3, 5:7), "]"))
vars_out <- c("educ", "city")                                  # beta[2], beta[3]
vars_sel <- c("hwage", "youngkids", "feduc", "educ", "city")   # gamma[2,3,5,6,7]

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

df <- rbind(get_ci(fit_n, "SLn"), get_ci(fit_sn, "SLsn"))
df$model <- factor(df$model, levels = c("SLn", "SLsn"))
df$eq <- factor(ifelse(grepl("beta", df$param), "Outcome", "Selection"),
                levels = c("Outcome", "Selection"))

key <- data.frame(param = pars,
                  var   = c(vars_out, vars_sel),
                  pos   = c(8:7, 5:1))
df <- merge(df, key, by = "param")
df$y <- df$pos + ifelse(df$model == "SLn", 0.18, -0.18)   # SLn above, SLsn below

p_est <- ggplot(df, aes(y = y, colour = model, shape = model)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_segment(aes(x = l95, xend = u95, yend = y), linewidth = 0.5) +
  geom_segment(aes(x = l80, xend = u80, yend = y), linewidth = 1.6) +
  geom_point(aes(x = med), size = 3) +
  scale_y_continuous(breaks = key$pos, labels = key$var,
                     expand = expansion(add = 0.6)) +
  scale_colour_manual(values = cols) +
  scale_shape_manual(values = shapes) +
  facet_grid(eq ~ ., scales = "free", space = "free_y", switch = "y") +
  labs(x = "Posterior estimate", y = NULL, colour = NULL, shape = NULL) +
  base_theme +
  theme(strip.placement   = "outside",
        strip.background  = element_rect(fill = "grey90"),
        strip.text.y.left = element_text(angle = 90, face = "bold", size = 15))

## =================================================================
## (b) Quantile residuals  r_i = Phi^{-1}( F(y_i | D_i = 1) )
## =================================================================
## data: replace with the objects you passed to Stan if different
library(sampleSelection)
data(Mroz87)
d   <- Mroz87
sel <- d$lfp == 1
Zs  <- cbind(1, d$huswage, d$kids5, d$mtr, d$fatheduc, d$educ, d$city)[sel, ]
X   <- cbind(1, d$educ, d$city)[sel, ]
y   <- log(d$wage[sel])

dens_sel <- function(res, zg, sigma, lambda, rho) {
  omega <- sqrt(sigma^2 + lambda^2)
  lt    <- -lambda * rho / omega
  am    <- sigma * rho / omega^2
  ak    <- lambda / (sigma * omega)
  O11   <- 1 - rho^2 + lt^2
  ra    <- -lt / sqrt(O11)
  a <- pmin(pmax((zg + am * res) / sqrt(O11), -37), 37)
  b <- pmin(pmax(ak * res, -37), 37)
  2 * dnorm(res, 0, omega) * pbivnorm(a, b, -ra)
}

cdf_sel <- function(r, zg, sigma, lambda, rho) {
  pz <- pnorm(zg)
  if (r <= 0) {
    num <- integrate(dens_sel, -Inf, r, zg = zg, sigma = sigma,
                     lambda = lambda, rho = rho, rel.tol = 1e-8)$value
  } else {
    num <- pz - integrate(dens_sel, r, Inf, zg = zg, sigma = sigma,
                          lambda = lambda, rho = rho, rel.tol = 1e-8)$value
  }
  min(max(num / pz, 1e-12), 1 - 1e-12)
}

qres <- function(fit, skew, ndraw = 200) {
  pd  <- rstan::extract(fit)
  S   <- length(pd$rho)
  idx <- round(seq(1, S, length.out = min(ndraw, S)))
  U   <- matrix(NA_real_, length(y), length(idx))
  for (k in seq_along(idx)) {
    s      <- idx[k]
    sigma  <- sqrt(pd$sigma2[s])
    lambda <- if (skew) pd$lambda[s] else 0
    res    <- drop(y - X %*% pd$beta[s, ])
    zg     <- drop(Zs %*% pd$gamma[s, ])
    U[, k] <- mapply(cdf_sel, res, zg,
                     MoreArgs = list(sigma = sigma, lambda = lambda, rho = pd$rho[s]))
  }
  qnorm(rowMeans(U))
}

# skip if r_n / r_sn already exist in the session
if (!exists("r_n"))  r_n  <- qres(fit_n,  skew = FALSE)
if (!exists("r_sn")) r_sn <- qres(fit_sn, skew = TRUE)

rd <- data.frame(r     = c(r_n, r_sn),
                 model = factor(rep(c("SLn", "SLsn"), each = length(y)),
                                levels = c("SLn", "SLsn")))

p_qq <- ggplot(rd, aes(sample = r, colour = model, shape = model)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "grey40") +
  stat_qq(size = 2, alpha = 0.7) +
  scale_colour_manual(values = cols) +
  scale_shape_manual(values = shapes) +
  labs(x = "Theoretical quantiles", y = "Quantile residuals",
       colour = NULL, shape = NULL) +
  base_theme

## =================================================================
## combine: one shared legend on top, equal-width (near-square) panels
## =================================================================
p_all <- (p_est | p_qq) +
  plot_layout(guides = "collect", widths = c(1, 1)) &
  theme(legend.position = "top")

p_all
ggsave("app1_compare_residual.pdf", p_all, width = 12, height = 6.5)
