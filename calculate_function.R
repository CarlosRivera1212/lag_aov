# --------------------------- #
# Function to calculate various statistics based on input data and matrices
f_calc <- function(X, Z, W, y, trt, blq, ntot) {
  # Calculate the reparametrized design matrix X0
  X0 = X %*% Z
  
  # Calculate the perpendicular projection matrix M0
  M0 = X0 %*% solve(t(X0) %*% X0) %*% t(X0)
  
  # Perform ANOVA
  mod_aov = aov(y ~ trt + blq)
  smod_aov = summary(mod_aov)[[1]]
  F_aov = smod_aov[1, 4]
  pv_aov = smod_aov[1, 5]
  
  # Calculate statistics related to ALFA
  H = W %*% M0 %*% y
  normH = norm(H)
  MH = H %*% solve(t(H) %*% H) %*% t(H)
  
  I = diag(1, ntot)
  f1v = MH %*% (I - M0) %*% y   # as a vector
  f1e = norm(f1v)               # as a scalar
  
  I_M0y = (I - M0) %*% y
  
  fee = norm(I_M0y)
  
  vrz = (1 / ntot) * t(y) %*% I_M0y   # variance
  
  yM0W = t(y) %*% M0 %*% W
  
  # Estimator of overlap
  mu = y
  
  alfa = t(y) %*% M0 %*% t(W) %*% (I - M0) %*% (y) / norm(W %*% M0 %*% y)
  
  U = M0 %*% W %*% (I - M0)
  U_ = (t(U) + U) / 2
  Us = forceSymmetric(U, uplo = 'U')
  Ul = forceSymmetric(U, uplo = 'L')
  
  V = M0 %*% W %*% W %*% M0
  V_ = (t(V) + V) / 2
  Vs = forceSymmetric(V, uplo = 'U')
  Vl = forceSymmetric(V, uplo = 'L')
  
  media_As = t(mu) %*% Us %*% mu
  media_Bs = vrz * matrix.trace(W %*% M0 %*% W) + t(mu) %*% Vs %*% mu
  
  var_B_m = M0 %*% W %*% W %*% M0 %*% W %*% W %*% M0
  var_B = 2 * vrz ** 2 * matrix.trace(var_B_m) + 4 * vrz * t(mu) %*% var_B_m %*% mu
  
  media_alfas = media_As * (media_Bs ** 2 + var_B) / media_Bs ** 3
  var_alfas = (media_As ** 2 / media_Bs ** 4) * var_B
  
  media_Al = t(mu) %*% Ul %*% mu
  media_Bl = vrz * matrix.trace(W %*% M0 %*% W) + t(mu) %*% Vl %*% mu
  
  media_alfal = media_Al * (media_Bl ** 2 + var_B) / media_Bl ** 3
  var_alfal = (media_Al ** 2 / media_Bl ** 4) * var_B
  
  media_A_ = t(mu) %*% U_ %*% mu
  media_B_ = vrz * matrix.trace(W %*% M0 %*% W) + t(mu) %*% V_ %*% mu
  
  media_alfa_ = media_A_ * (media_B_ ** 2 + var_B) / media_B_ ** 3
  var_alfa_ = (media_A_ ** 2 / media_B_ ** 4) * var_B
  
  # Calculate Moran's I statistics
  IM = Moran.I(y, W)
  pmoran = IM$p.value
  IMaov = Moran.I(mod_aov$residuals, W)
  pmoran_aov = IMaov$p.value
  
  # Calculate the OT statistic
  OT = f_calc_OT(X0, y, W)
  
  # Calculate the chi-squared p-value for the OT statistic
  prchi = pchisq(q = OT, df = 1, lower.tail = FALSE)
  
  # Create a list of computed statistics
  sal = c(
    'OT' = OT,
    'alfa' = alfa,
    'media_As' = media_As[1,1],
    'media_Bs' = media_Bs[1,1],
    'media_alfas' = media_alfas[1,1],
    'var_alfas' = var_alfas[1,1],
    'media_Al' = media_Al[1,1],
    'media_Bl' = media_Bl[1,1],
    'media_alfal' = media_alfal[1,1],
    'var_alfal' = var_alfal[1,1],
    'media_A_' = media_A_[1,1],
    'media_B_' = media_B_[1,1],
    'media_alfa_' = media_alfa_[1,1],
    'var_alfa_' = var_alfa_[1,1],
    'F_aov' = F_aov,
    'pv_aov' = pv_aov,
    'pmoran' = pmoran,
    'pmoran_aov' = pmoran_aov,
    'prchi' = prchi
  )
  
  return(sal)
}

# --------------------------- #
# Function to calculate the OT statistic
f_calc_OT <- function(X0, y, W) {
  # Create an identity matrix
  I = diag(1, length(y))
  
  # Calculate the perpendicular projection matrix M0
  M0 = X0 %*% solve(t(X0) %*% X0) %*% t(X0)
  
  # Calculate I_M0y
  I_M0y = (I - M0) %*% y
  
  # Calculate vrz (variance)
  vrz = (1 / length(y)) * t(y) %*% I_M0y
  
  # Calculate H
  H = W %*% M0 %*% y
  
  # Recalculate M0
  M0 = X0 %*% solve(t(X0) %*% X0) %*% t(X0)
  
  # Calculate MH
  MH = H %*% solve(t(H) %*% H) %*% t(H)
  
  # Calculate the correlation matrix COR
  COR = I + X0 %*% solve(t(X0) %*% (I - MH) %*% X0) %*% t(X0)
  
  # Calculate P1
  P1 = (1 / vrz) %*% t(y) %*% (I - M0) %*% MH
  
  # Calculate P2
  P2 = MH %*% (I - M0) %*% y
  
  # Calculate OT statistic
  OT = P1 %*% COR %*% P2
  
  return(OT[1, 1])
}
