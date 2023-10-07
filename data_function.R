# --------------------------- #
# Function to generate simulated data
f_data <- function(ntrt, nblq, nrep, nsim) {
  # Calculate the total number of observations
  ntot = ntrt * nblq * nrep
  
  # Generate a grid of x and y values for spatial locations
  dataxy = expand.grid(x = seq(ntrt * nblq), y = seq(nrep))
  
  # Randomly sample a fraction of spatial locations
  dataxy = sample_frac(dataxy)
  
  # Create a data frame with treatment, block, and replication factors
  data0 = data.frame(
    trt = gl(ntrt, ntot / ntrt, ntot, paste0('T', 1:ntrt)),
    blq = gl(nblq, nrep, ntot, paste0('B', 1:nblq)),
    rep = seq(nrep)
  )
  
  # Create a dummy geostatistical model for spatial predictions
  g.dummy = gstat(formula=z~1, locations=~x+y, dummy=T, beta=2.0,
                  model=vgm(psill=0.2, range=5, model='Gau'), nmax=20)
  
  # Predict values using the dummy geostatistical model
  rto = predict(g.dummy, newdata=dataxy, nsim=nsim)
  
  # Combine the generated data with predictions and arrange it
  data_exp = data0 %>% 
    sample_frac() %>% 
    bind_cols(rto) %>% 
    arrange(trt, blq, rep)
  
  return(data_exp)
}


# --------------------------- #
# Function to reorder data based on a random sample of rows
f_data_ord <- function(data, nr){
  # Separate the data into coordinates (x and y) and other columns
  data0 = select(data, -c(x,y))
  dataxy = select(data, c(x,y))
  
  # Randomly sample rows and reorder them
  row_sel = sample(1:nrow(data), nr)
  row_ord = sample(row_sel)
  dataxy[row_sel, ] = dataxy[row_ord,]
  
  # Combine the reordered coordinates with other columns
  data_s = bind_cols(dataxy, data0)
  return(data_s)
}
