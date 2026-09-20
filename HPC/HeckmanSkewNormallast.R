functions {
  // 10-point Gauss–Legendre nodes/weights on [-1, 1]
  vector gl10_x() {
    vector[10] x;
    x[1]  = -0.9739065285171717;
    x[2]  = -0.8650633666889845;
    x[3]  = -0.6794095682990244;
    x[4]  = -0.4333953941292472;
    x[5]  = -0.1488743389816312;
    x[6]  =  0.1488743389816312;
    x[7]  =  0.4333953941292472;
    x[8]  =  0.6794095682990244;
    x[9]  =  0.8650633666889845;
    x[10] =  0.9739065285171717;
    return x;
  }

  vector gl10_w() {
    vector[10] w;
    w[1]  = 0.0666713443086881;
    w[2]  = 0.1494513491505806;
    w[3]  = 0.2190863625159820;
    w[4]  = 0.2692667193099963;
    w[5]  = 0.2955242247147529;
    w[6]  = 0.2955242247147529;
    w[7]  = 0.2692667193099963;
    w[8]  = 0.2190863625159820;
    w[9]  = 0.1494513491505806;
    w[10] = 0.0666713443086881;
    return w;
  }

  // Approximation to standard bivariate normal CDF:
  // Phi2(a,b;rho) = P(Z1<=a, Z2<=b), where corr(Z1,Z2)=rho
  // Uses Gauss–Legendre quadrature (no integrate_1d and no MVN lcdf/cdf needed).
  real my_bvnorm_std(real a, real b, real rho) {
    real rho_c = fmin(0.999999, fmax(-0.999999, rho));

    // Handle negative correlation via identity:
    // Phi2(a,b;rho) = Phi(a) - Phi2(a,-b;-rho), for rho<0
    if (rho_c < 0)
      return Phi(a) - my_bvnorm_std(a, -b, -rho_c);

    // Independence
    if (abs(rho_c) < 1e-12)
      return Phi(a) * Phi(b);

    real theta = asin(rho_c);  // in (0, pi/2)
    real base  = Phi(a) * Phi(b);

    vector[10] x = gl10_x();
    vector[10] w = gl10_w();

    // Gauss–Legendre on [0, theta]
    real half = 0.5 * theta;
    real sum_int = 0.0;

    for (i in 1:10) {
      real phi  = half * (x[i] + 1);  // map [-1,1] -> [0,theta]
      real s    = sin(phi);
      real c    = cos(phi);
      real expo = - (a*a + b*b - 2*a*b*s) / (2 * c*c);
      sum_int  += w[i] * exp(expo);
    }

    return base + (half * sum_int) / (2 * pi());
  }

  // General bivariate normal CDF with means, sds, and correlation
  real my_bvnorm(real y1, real y2,
                  real mu1, real mu2,
                  real s1,  real s2,
                  real rho) {
    real a = (y1 - mu1) / s1;
    real b = (y2 - mu2) / s2;
    return my_bvnorm_std(a, b, rho);
  }
}


data {
  // dimensions
  int<lower=1> N;
  int<lower=1, upper=N> N_y;
  int<lower=1> p;
  int<lower=1> q;
  // covariates
  matrix[N_y, p] X;
  matrix[N, q] Z;
  // responses
  array[N] int<lower=0, upper=1> D; //int<lower=0, upper=1> D[N];
  vector[N_y] y;
}
parameters {
  vector[p] beta;
  vector[q] gamma;
  real rho_raw;
  real<lower=0> sigma;
  real lambda;
}

transformed parameters {
  real<lower=0> sigma2 = sigma^2;
  real rho = tanh(rho_raw);
  real lambdat = -lambda * rho / sqrt(sigma2 + lambda^2);
  real auxmt = sigma * rho / (sigma2 + lambda^2);
  real auxkt = lambda / (sigma * sqrt(sigma2 + lambda^2));
  real<lower=0> Omega11 =
      (1 - rho^2) + lambdat^2;
  real rhoaaux =
      -lambdat / sqrt(Omega11);
}

model {
  // naive (truncated) priors
  beta ~ multi_normal(rep_vector(0, p), diag_matrix(rep_vector(100, p)));
  gamma ~ multi_normal(rep_vector(0, q), diag_matrix(rep_vector(100, q)));
  rho_raw ~ normal(0, 1);
  sigma ~ normal(0, 2);
  lambda ~ normal(0, 5);
  //rho ~ normal(murh,sigmarh);
  //sigma ~ cauchy(0, 4);
  //lambda ~ cauchy(0, 4);
  // murh ~ normal(0, 1);
  //sigmarh ~ cauchy(0, 4);
  {
    // log-likelihood
    vector[N_y] Xb = X * beta;
    vector[N] Zg = Z * gamma;
    int ny = 1;
    for(n in 1:N) {
      if(D[n] > 0) {
        real mut1 =  Zg[n]+auxmt*(y[ny]- Xb[ny]);
        real mut2 = auxkt*(y[ny]- Xb[ny]);
        real pp = my_bvnorm(-mut1, mut2, 0, 0, sqrt(Omega11), 1, rhoaaux);
        real aa = Phi(mut2);
        target += normal_lpdf(y[ny] | Xb[ny], sqrt(sigma2+lambda^2))+log(aa-pp)+log(2);
        ny += 1;
      }
      else {
        target += log(Phi(-Zg[n])); //log(Phi(-Zg[n]));
      }
    }
  }
}
generated quantities {

  vector[N] log_lik;
  vector[N_y] Xb = X * beta;
  vector[N] Zg = Z * gamma;

  real log_lik_total;
  real AIC;
  real BIC;

  int ny = 1;

  for (n in 1:N) {

    if (D[n] > 0) {

      real mut1 =
        Zg[n] + auxmt * (y[ny] - Xb[ny]);

      real mut2 =
        auxkt * (y[ny] - Xb[ny]);

      real pp =
        my_bvnorm(
          -mut1,
          mut2,
          0,
          0,
          sqrt(Omega11),
          1,
          rhoaaux
        );

      real aa = Phi(mut2);

      log_lik[n] =
        normal_lpdf(
          y[ny] |
          Xb[ny],
          sqrt(sigma2 + lambda^2)
        )
        + log(aa - pp)
        + log(2);

      ny += 1;

    } else {

      log_lik[n] =
        log(Phi(-Zg[n]));
    }
  }

  log_lik_total = sum(log_lik);

  AIC =
    -2 * log_lik_total
    + 2 * (p + q + 3);

  BIC =
    -2 * log_lik_total
    + log(N) * (p + q + 3);
}
