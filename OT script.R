# Load required libraries
library(ggplot2)
library(dplyr)
library(ape)
library(matrixcalc)
library(Matrix)
library(gstat)

# Set working directory
setwd("G:/My Drive/PC/UNAL/D/articulos/rezago_espacial/git")

# Load custom functions from source files
source('data_function.R')            # Load f_data() and f_data_ord() functions
source('WXZ_matriz_function.R')      # Load f_W(), f_X(), f_Z() functions
source('calculate_function.R')       # Load f_calc() and f_calc_OT() functions

# Define parameters
ntrt = 3
nblq = 2
nrep = 10
ntot = ntrt * nblq * nrep
nsim = 1000

# Set random seed for reproducibility
set.seed(123)

# Generate data using f_data() function
data = f_data(ntrt, nblq, nrep, nsim)

# Create a ggplot2 plot (plt1) based on the generated data
plt1 = ggplot(data) +
  aes(x, y, fill = sim1, label = paste(trt, blq, rep, sep = '.')) +
  geom_tile() +
  geom_text() +
  scale_fill_gradientn(colours = terrain.colors(10, rev = TRUE))

# Define a sequence of numbers for data ordering
num_ord = seq(2, nrow(data), 2)

# Create ordered data using f_data_ord() function
data_ord = lapply(num_ord, function(nr) {
  f_data_ord(data, nr)
})

# Print the length and dimensions of data_ord[[1]]
length(data_ord)
dim(data_ord[[1]])

# Combine ordered data with the original data
data_ord_tot = c(list(data), data_ord)

# Check the class and length of data_ord_tot
class(data_ord_tot)
length(data_ord_tot)

# Create the design matrix X
X = f_X(ntrt, nblq, nrep)

# Rank the design matrix and extract the rank
rnk = rankMatrix(X)[1]

# Create the reparametrization matrix Z
Z = f_Z(ntrt, nblq, rnk)

# --------------------------- #
# SIMULACION
alfas = lapply(seq(length(data_ord_tot)), function(i) {
  cat(i, '| ')
  
  # Extract data for the current iteration
  data_i = data_ord_tot[[i]]
  trt_i = data_i$trt
  blq_i = data_i$blq
  
  # Calculate the weight matrix 'W_i' using the 'f_W' function
  W_i = f_W(select(data_i, c(x, y)))
  
  # Prepare the design matrix 'Xo' by combining 'X' with 'x' and 'y' columns
  Xo = data_i %>% 
    select(x, y) %>% 
    cbind(X) %>%
    arrange(y, x) %>%
    select(-c(x, y)) %>%
    as.matrix()
  
  calc_i = sapply(seq(nsim), function(i){
    y_i = data_i[, paste0('sim', i)]
    f_calc(Xo, Z, W_i, y_i, trt_i, blq_i, ntot)
  })
  
  sal = as.data.frame(t(calc_i))
  
  return(sal)
})

# --------------------------- #
# Assign names to the 'alfas' list using 'paste0' function
names(alfas) = paste0('ord', c(0, num_ord))

# Convert the list of data frames 'alfas' into a single data frame 'alfas_df' using 'plyr::ldply'
alfas_df = plyr::ldply(alfas)

# Add a new column 'name' to 'alfas_df' containing factors based on '.id' (index) and names(alfas) with ordered levels
alfas_df = alfas_df %>% 
  mutate(name = factor(.id, names(alfas), ordered = TRUE))

# nf = paste0('trt',ntrt,'blq',nblq,'rep',nrep,'sim',nsim)
# saveRDS(data_ord_tot, paste0('sim_',nf,'.rds'))
# saveRDS(alfas_df, paste0('alfas_df_',nf,'.rds'))
