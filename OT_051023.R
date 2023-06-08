library(ggplot2)
library(dplyr)
library(ape)
library(matrixcalc)
library(Matrix)
library(gstat)
library(parallel)

setwd('~/CARLOS/UNAL/D/articulos/rezago_espacial')

f_data <- function(ntrt, nblq, nrep, nsim) {
  ntot = ntrt * nblq * nrep
  
  dataxy = expand.grid(x = seq(ntrt * nblq), y = seq(nrep))
  dataxy = sample_frac(dataxy)
  
  data0 = data.frame(
    trt = gl(ntrt, ntot / ntrt, ntot, paste0('T', 1:ntrt)),
    blq = gl(nblq, nrep, ntot, paste0('B', 1:nblq)),
    rep = seq(nrep)
  )
  
  g.dummy = gstat(formula=z~1, locations=~x+y, dummy=T, beta=2.0,
                  model=vgm(psill=0.2, range=5, model='Gau'), nmax=20)
  # rto_1 = predict(g.dummy, newdata=dataxy[,c('x','y')], nsim=nsim)
  rto = predict(g.dummy, newdata=dataxy, nsim=nsim)
  # rto[,-c(1,2)] = 100*rto[, -c(1,2)]
  
  data_exp = data0 %>% 
    sample_frac() %>% 
    bind_cols(rto) %>% 
    arrange(trt, blq, rep)
  
  return(data_exp)
}

f_data_ord <- function(data, nr){
  data0 = select(data, -c(x,y))
  dataxy = select(data, c(x,y))
  
  row_sel = sample(1:nrow(data), nr)
  row_ord = sample(row_sel)
  dataxy[row_sel, ] = dataxy[row_ord,]
  
  data_s = bind_cols(dataxy, data0)
  return(data_s)
}

f_W <- function(data_xy) {
  papa.dists = as.matrix(dist(data_xy))
  # selección de vecinos más cercanos
  # papa.dists[papa.dists > sqrt(2)] = 0
  # papa.dists[papa.dists > 10*sqrt(2)] = 0
  # inverso de las distancias
  papa.dists.inv = 1 / papa.dists
  # anulando algunas distancias
  diag(papa.dists.inv) = 0
  # anulando infinitos
  papa.dists.inv[is.infinite(papa.dists.inv)] = 0
  # matriz de distancias inversas
  W0 = papa.dists.inv
  # estandarizando la matriz para generar pesos
  #sumaf=apply(W ,1,sum); W=W/sumaf
  W = W0 / sum(W0)
  # matriz de pesos
  
  # mapeando la matriz de pesos
  # image(W, asp = 1)
  
  return(W)
}

f_X <- function(ntrt, nblq, nrep){
  m_zero = matrix(0, nrep*ntrt*nblq, 1+ntrt+nblq)
  m_zero[,1] = 1
  
  nrep_trt = nrep*nblq
  
  # tratamientos
  for(i in seq(ntrt)){
    pos_trt = 1:nrep_trt+((i-1)*nrep_trt)
    m_zero[pos_trt, i+1] = 1
  }
  
  # bloques
  nrep_blq = nrep
  for(i in seq(nblq)){
    pos_blq = 1:nrep_blq+(i-1)*nrep_blq+rep(0:(ntrt-1), each=nrep_blq)*nrep_trt
    m_zero[pos_blq, i+ntrt+1] = 1
  }
  
  return(m_zero)
}

f_Z <- function(ntrt, nblq, rnk){
  f1 = matrix(1, 1, rnk)
  
  val = seq(-ntrt*nblq, ntrt*nblq)
  f_trt = matrix(sample(val, rnk*(ntrt-1), F), ntrt-1, byrow = T)
  f_trt = rbind(f_trt, -colSums(f_trt))
  
  f_blq = matrix(sample(val, rnk*(nblq-1), F), nblq-1, byrow = T)
  f_blq = rbind(f_blq, -colSums(f_blq))
  
  mz = rbind(f1, f_trt, f_blq)
  
  if (rankMatrix(mz)[1] != rnk) {
    mz = f_Z(ntrt, nblq, rnk)
  }
  else {
    return(mz)
  }
}
f_calc_OT <- function(X0, y, W){
  I = diag(1, length(y))
  
  M0 = X0 %*% solve(t(X0) %*% X0) %*% t(X0)
  I_M0y = (I - M0) %*% y
  vrz = (1 / length(y)) * t(y) %*% I_M0y
  
  H = W %*% M0 %*% y
  M0 = X0 %*% solve(t(X0) %*% X0) %*% t(X0)
  MH = H %*% solve(t(H) %*% H) %*% t(H)
  
  COR = I + X0 %*% solve(t(X0) %*% (I - MH) %*% X0) %*% t(X0)
  P1 = (1 / vrz) %*% t(y) %*% (I - M0) %*% MH
  P2 = MH %*% (I - M0) %*% y
  OT = P1 %*% COR %*% P2
  return(OT[1,1])
}
f_calc <- function(X, Z, W, y, trt, blq, ntot){
  # la matriz de diseño reparametrizada
  X0 = X %*% Z
  # la matriz de proyección perpendicular de X0
  M0 = X0 %*% solve(t(X0) %*% X0) %*% t(X0)
  
  # ANOVA
  mod_aov = aov(y ~ trt + blq)
  smod_aov = summary(mod_aov)[[1]]
  F_aov = smod_aov[1, 4]
  pv_aov = smod_aov[1, 5]
  
  # ALFA
  # una matriz intermedia H
  H = W %*% M0 %*% y
  # la norma de H
  normH = norm(H)
  # la matriz de proyección perpendicular de H
  MH = H %*% solve(t(H) %*% H) %*% t(H)
  
  # calculos intermedios
  I = diag(1, ntot)
  f1v = MH %*% (I - M0) %*% y   #como vector
  f1e = norm(f1v)               #como escalar
  
  I_M0y = (I - M0) %*% y
  
  fee = norm(I_M0y)
  
  vrz = (1 / ntot) * t(y) %*% I_M0y   #varianza
  
  yM0W = t(y) %*% M0 %*% W
  
  # estimador de solapamiento
  # mu = rep(mean(y), length(y))
  mu = y
  
  alfa = t(y) %*% M0 %*% t(W) %*% (I - M0) %*% (y) / norm(W %*% M0 %*% y)
  # media_A = vrz * matrix.trace(U) + t(y) %*% U %*% y
  # media_B = vrz * matrix.trace(V) + t(y) %*% V %*% y
  
  U = M0 %*% W %*% (I - M0)
  U_ = (t(U) + U) / 2
  Us = forceSymmetric(U, uplo = 'U')
  Ul = forceSymmetric(U, uplo = 'L')
  
  V = M0 %*% W %*% W %*% M0
  V_ = (t(V) + V) / 2
  Vs = forceSymmetric(V, uplo = 'U')
  Vl = forceSymmetric(V, uplo = 'L')
  
  media_As = t(mu) %*% Us %*% mu
  # media_Bs = var(y) * matrix.trace(W %*% M0 %*% W) + t(mu) %*% V %*% mu
  media_Bs = vrz * matrix.trace(W %*% M0 %*% W) + t(mu) %*% Vs %*% mu
  
  var_B_m = M0 %*% W %*% W %*% M0 %*% W %*% W %*% M0
  # var_B = 2 * var(y) ** 2 * matrix.trace(var_B_m) + 4 * var(y) * t(mu) %*% var_B_m %*% mu
  var_B = 2 * vrz ** 2 * matrix.trace(var_B_m) + 4 * vrz * t(mu) %*% var_B_m %*% mu
  
  media_alfas = media_As * (media_Bs ** 2 + var_B) / media_Bs ** 3
  var_alfas = (media_As**2 / media_Bs**4) * var_B
  
  
  media_Al = t(mu) %*% Ul %*% mu
  media_Bl = vrz * matrix.trace(W %*% M0 %*% W) + t(mu) %*% Vl %*% mu

  media_alfal = media_Al * (media_Bl ** 2 + var_B) / media_Bl ** 3
  var_alfal = (media_Al**2 / media_Bl**4) * var_B
  
  
  media_A_ = t(mu) %*% U_ %*% mu
  media_B_ = vrz * matrix.trace(W %*% M0 %*% W) + t(mu) %*% V_ %*% mu
  
  media_alfa_ = media_A_ * (media_B_ ** 2 + var_B) / media_B_ ** 3
  var_alfa_ = (media_A_**2 / media_B_**4) * var_B
  
  
  # MORAN
  IM = Moran.I(y, W)
  pmoran = IM$p.value
  IMaov = Moran.I(mod_aov$residuals, W)
  pmoran_aov = IMaov$p.value

  ########
  # H = W %*% M0 %*% y
  # P1 = t(y) %*% (I - M0) %*% H
  # la matriz G
  # G = solve(t(X0) %*% (I - MH) %*% X0)
  # MH = H %*% solve(t(H) %*% H) %*% t(H)
  
  # el estadístico de solapamiento
  # OT = (1 / (vrz * normH)) * P1**2 * (1 + normH**(-2) * t(H) %*% X0 %*% G %*% t(X0) %*% H)
  OT = f_calc_OT(X0, y, W)
  # H'*(I+X0*G*X0')*H
  
  # la distribución del estadístcio
  prchi = pchisq(q = OT, df = 1, lower.tail = F)
  
  sal = c(
    'OT' = OT,
    'alfa' = alfa,
    # 'var_B' = var_B,
    
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

ntrt = 3
nblq = 2
nrep = 10
ntot = ntrt*nblq*nrep
nsim = 1000

set.seed(123)
data = f_data(ntrt, nblq, nrep, nsim)
plt1 = ggplot(data)+
  aes(x,y, fill=sim1,
      label = paste(trt, blq, rep, sep = '.'))+
  geom_tile()+
  geom_text()+
  scale_fill_gradientn(colours = terrain.colors(10, rev=T))
# 
# ggsave('plt1.tif', plt1, device = 'tiff',
#        scale = 1, width = 16, height = 9, units = 'cm')

num_ord = seq(2, nrow(data), 2)
data_ord = lapply(num_ord, function(nr){
  f_data_ord(data, nr)
})
length(data_ord)
dim(data_ord[[1]])

data_ord_tot = c(list(data), data_ord)
class(data_ord_tot)
length(data_ord_tot)

# la matriz de diseño
X = f_X(ntrt, nblq, nrep)
rnk = rankMatrix(X)[1]
# la matriz de reparametrización
# Z=matrix(c(1,0,-1,0,0,-1,1,-2,0,0,-2,1,0,1,1,1,0,1,0,-1,0,-1,0,1),ncol=4,byrow=T)
Z = f_Z(ntrt, nblq, rnk)

# # # # # # # # # # # # # # # # # # # # # # # # # # # 
# data_i_1 = data
data_i_1 = data_ord_tot[[1]]
tibble(data_i_1)

Xo = bind_cols(X, select(data_i_1, x,y)) %>%
  arrange(x,y) %>%
  select(-c(x,y)) %>%
  as.matrix()


trt_i_1 = data_i_1$trt
blq_i_1 = data_i_1$blq
W_i_1 = f_W(select(data_i_1, c(x,y)))

calc_i_1 = sapply(seq(nsim), function(i){
  y_i_1 = data_i_1[ , paste0('sim', i)]
  f_calc(Xo, Z, W_i_1, y_i_1, trt_i_1, blq_i_1, ntot)
})

calc_i_1 = as.data.frame(t(calc_i_1))
ggplot(data_i_1)+
  aes(x,y, fill=sim1)+
  geom_tile()

ggplot(calc_i_1)+
  aes(OT)+
  geom_histogram()

# SIMULACION
nc = detectCores()
cl = makeCluster(nc-1)
clusterEvalQ(cl, {
  library(dplyr)
  library(ape)
  library(matrixcalc)
  library(Matrix)
  library(gstat)
})
clusterExport(cl, c('f_calc', 'f_W', 'f_calc_OT', 'X', 'Z', 'nsim', 'ntot'))

alfas = parLapply(cl, data_ord_tot, function(data_i) {
  trt_i = data_i$trt
  blq_i = data_i$blq
  W_i = f_W(select(data_i, c(x,y)))
  Xo = data_i %>% 
    select(x,y) %>% 
    bind_cols(X) %>%
    arrange(y,x) %>%
    select(-c(x,y)) %>%
    as.matrix()
  
  calc_i = sapply(seq(nsim), function(i){
    y_i = data_i[ , paste0('sim', i)]
    f_calc(Xo, Z, W_i, y_i, trt_i, blq_i, ntot)
  })
  sal = as.data.frame(t(calc_i))
  
  return(sal)
})
stopCluster(cl)

# # # # # # # # # # # # # # # # # # # # # # # # # # # 
names(alfas) = paste0('ord', c(0, num_ord))
alfas_df = plyr::ldply(alfas)

alfas_df = alfas_df %>% 
  mutate(name = factor(.id, names(alfas), ordered = T))

# nf = paste0('trt',ntrt,'blq',nblq,'rep',nrep,'sim',nsim)
# saveRDS(data_ord_tot, paste0('sim_',nf,'_052023.rds'))
# saveRDS(alfas_df, paste0('alfas_df_',nf,'_052023.rds'))
# # # # # # # # # # # # # # # # # # # # # # # # # # # 

int_sel = c(0, 12, 24, 36, 48, 60)

plot(alfas_df$media_As, alfas_df$media_Al, pch=16, cex=.1)
plot(alfas_df$media_alfas, alfas_df$media_alfal, pch=16, cex=.1)

alfas_df %>% 
  filter(name %in% paste0('ord', int_sel)) %>% 
  select(name, media_alfa_, media_alfas, media_alfal) %>% 
  tidyr::pivot_longer(-name,
                      names_to = 'var',
                      values_to = 'value') %>% 
  ggplot()+
  aes(value, fill=name)+
  geom_density(alpha=0.5)+
  facet_wrap(~var)

alfas_df %>% 
  filter(name %in% paste0('ord', int_sel)) %>% 
  select(name, var_alfa_, var_alfas, var_alfal) %>% 
  tidyr::pivot_longer(-name,
                      names_to = 'var',
                      values_to = 'value') %>% 
  ggplot()+
  aes(-log(value), fill=name)+
  geom_density(alpha=0.5)+
  facet_grid(name~var)


hist(alfas_df$media_As)
hist(alfas_df$media_Al)
hist(alfas_df$media_B)
range(alfas_df$var_B)


alfas_df %>% 
  filter(name %in% paste0('ord', int_sel)) %>% 
  ggplot()+
  # aes(sqrt(var_alfa))+
  aes(media_alfa)+
  geom_histogram()+
  facet_wrap(~name)+
  theme(axis.text.x = element_text(angle=90))

alfas_df %>% 
  filter(name %in% paste0('ord', int_sel)) %>% 
  ggplot()+
  aes(var_alfa, media_alfa)+
  geom_point()+
  geom_vline(aes(xintercept = mean(var_alfa)))+
  facet_wrap(~name)+
  theme(axis.text.x = element_text(angle=90))



alfas_df %>% 
  select(name, var_alfa) %>% 
  mutate(n_sim = rep(1:nsim, 31)) %>% 
  pivot_wider(names_from = name, values_from = var_alfa)


range(alfas_df$media_alfa)
hist(alfas_df$media_alfa)

ggplot(alfas_df)+
  aes(pmoran, pmoran_aov)+
  geom_point(size=0.1)+
  geom_hline(yintercept = 0.05, color='red')+
  geom_vline(xintercept = 0.05, color='red')

tbl = table(alfas_df$pmoran_aov < 0.05,
            alfas_df$pmoran < 0.05)
100 * (tbl[1,2] + tbl[2,1]) / sum(tbl)

alfas_df %>% 
  group_by(.id, y = pmoran<0.05, r = pmoran_aov<0.05) %>% 
  summarise(n = n()) %>% 
  filter(.id=='ord30')
