# --------------------------- #
# Function to generate a weight matrix based on pairwise distances
f_W <- function(data_xy) {
  # Calculate pairwise Euclidean distances between data points
  d = as.matrix(dist(data_xy))
  
  # Calculate the inverse of distances (reciprocal)
  di = 1 / d
  
  # Handle infinite distances (division by zero)
  di[is.infinite(di)] = 0
  
  # Initialize the weight matrix W0 with the inverse distances
  W0 = di
  
  # Standardize the weight matrix to generate weights that sum to 1
  # W = W0 / apply(W0, 1, sum)
  W = W0 / sum(W0)
  
  return(W)
}


# --------------------------- #
# Function to create a design matrix for a mixed-effects model
f_X <- function(ntrt, nblq, nrep) {
  # Create an empty matrix filled with zeros
  m_zero = matrix(0, nrep * ntrt * nblq, 1 + ntrt + nblq)
  
  # Set the first column of the matrix to all ones
  m_zero[, 1] = 1
  
  # Calculate the total number of replicates per treatment
  nrep_trt = nrep * nblq
  
  # Loop through treatments and mark corresponding columns with ones
  for (i in seq(ntrt)) {
    pos_trt = 1:nrep_trt + ((i - 1) * nrep_trt)
    m_zero[pos_trt, i + 1] = 1
  }
  
  # Calculate the total number of replicates per block
  nrep_blq = nrep
  
  # Loop through blocks and mark corresponding columns with ones
  for (i in seq(nblq)) {
    pos_blq = 1:nrep_blq + (i - 1) * nrep_blq + rep(0:(ntrt - 1), each = nrep_blq) * nrep_trt
    m_zero[pos_blq, i + ntrt + 1] = 1
  }
  
  return(m_zero)
}

# --------------------------- #
# Function to create a design matrix for random effects
f_Z <- function(ntrt, nblq, rnk) {
  # Create a matrix of ones for the first row
  f1 = matrix(1, 1, rnk)
  
  # Generate random values for treatment effects
  val = seq(-ntrt * nblq, ntrt * nblq)
  f_trt = matrix(sample(val, rnk * (ntrt - 1), replace = FALSE), ntrt - 1, byrow = TRUE)
  f_trt = rbind(f_trt, -colSums(f_trt))
  
  # Generate random values for block effects
  f_blq = matrix(sample(val, rnk * (nblq - 1), replace = FALSE), nblq - 1, byrow = TRUE)
  f_blq = rbind(f_blq, -colSums(f_blq))
  
  # Combine treatment and block effects into the design matrix mz
  mz = rbind(f1, f_trt, f_blq)
  
  # Check if the rank of mz matches the desired rank (rnk)
  if (rankMatrix(mz)[1] != rnk) {
    # If not, recursively call the function to generate a new design matrix
    mz = f_Z(ntrt, nblq, rnk)
  } else {
    # If the rank matches the desired rank, return the design matrix
    return(mz)
  }
}

