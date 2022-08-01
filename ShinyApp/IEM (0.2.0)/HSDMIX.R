
# Preamble ----
# Homogeneous Surface Diffusion Model (HSDM) for ion exchange (IX)
#
# Levi M. Haupert (haupert.levi@epa.gov) and David G. Wahman (wahman.david@epa.gov)
# USEPA, 2022
#
# For an overview of ion exchange column modeling, consult:
#   Slater, M.J., 2013. Principles of ion exchange technology. Butterworth-Heinemann.
# and
# Helfferich, F. G. (1995). Ion exchange. Courier Corporation.
# 
# For details on the numerical method of solution (Orthogonal Collocation), consult:
#   Crittenden, J. C., Hutzler, N. J., Geyer, D. G., Oravitz, J. L., & Friedman, G. (1986). 
# Transport of organic compounds with saturated groundwater flow: Model development and 
# parameter sensitivity. Water Resources Research, 22(3), 271-284.
# 
# Assumptions:
#   Constant selectivity.
# Plug flow.
# Fickian diffusion.
# Gel type resins
#
# Ds, kL for each ion. Entries for the presaturant (reference ion) are ignored.
#
# NOTES:
#
#  Currently, ion order in "ions" and "Cin" tabs of input file must be the same.
#  The presaturant ion (which is also the reference ion) must be listed first.
#  The reference ion must be monovalent.
#
#
# TODO: automatic approximation of kL based on temperature, molar volume, flow velocity.
#       Not sure if better to do here, or in companion input generator app.
#       We do want the option to specify it (for optimization).
#
# TODO: Plotting output still needs lots of work.
#
# TODO: Performance optimization:
#       Jacobian is denser than if one Ds, kL. Might not be much room for optimization there.
#       Loops need vectorized for performance.
#       Function overhead makes sapply() and mapply() even slower than loops, though.
#
#
# CAUTION: This preliminary code has not been tested for quality or correctness.
#
# DISCLAIMER:
# The United States Environmental Protection Agency (EPA) code is provided on an "as is" basis and 
# the user assumes responsibility for its use. EPA has relinquished control of the information and 
# no longer has responsibility to protect the integrity , confidentiality, or availability of the information. 
# Any reference to specific commercial products, processes, or services by service mark, trademark, 
# manufacturer, or otherwise, does not constitute or imply their endorsement, recomendation or favoring by EPA. 
# The EPA seal and logo shall not be used in any manner to imply endorsement of any commercial product or activity 
# by EPA or the United States Government.


# Libraries ----

library(deSolve)
library(orthopolynom)
library(dplyr)
library(readxl)
library(ggplot2)

# Constants ----
S_PER_HR <- 60 * 60 # seconds per hour


# Inputs ----
nt_report = 201 # number of reporting steps
in_file_name <- "inp-file.xlsx"

# Load input file ----

params <- read_excel(in_file_name, sheet = "params")
ions <- read_excel(in_file_name, sheet = "ions")
Cin <- read_excel(in_file_name, sheet = "Cin")

# NOTE: All inputfile processing (loading, editing, saving) will be handled outisde of solve function.
#       Same for plotting and other output processing

# Collocation functions ----
rad_colloc <- function(N){
  # For a grid of N collocation points.
  # Calculate B (madrix operator for 1-D radial Laplacian for a symmetric sphere)
  # and W (vector Gauss-Radau quadrature weights)
  # Ref: Villadsen, J., & Michelsen, M. L. (1978)
  
  # calculate number of interior collocation points symmetric around x = 0
  N_int <- N - 1
  
  # setup roots
  # get list of recurrence relations for the Jacobi polynomial (0, 1)
  # "p" is on the interval of -1 to 1
  # "g" is on the interval of 0 to 1 (i.e., shifted)
  # 1,1 is shifted legendre from python with 0
  # 2.5, 1.5 is spherical symmetry
  # 2.0, 1.0 is cylinder symmetry
  # 1.5, 0.5 is slab symmetry
  p_list <- jacobi.g.recurrences(N_int, 2.5, 1.5)
  
  # using the recurrence relations, construct monic orthogonal polynomials
  m.r <- monic.polynomial.recurrences(p_list)
  
  # returns roots of the monic orthogonal polynomials
  # take square root as the problem is symmetrical and roots are taken as x^2 terms
  # terms at zero and 1
  roots_non_sym <- c(rev(polynomial.roots(m.r)[[N]]), 1)
  
  # create a data.frame to store values
  derivatives <- data.frame(
    roots = roots_non_sym,
    p_1 = rep(0, N),
    p_2 = rep(0, N),
    p_3 = rep(0, N)
  )
  
  # set initial values
  p_1 <- c(1, rep(0, N-1))
  p_2 <- rep(0, N)
  p_3 <- rep(0, N)
  
  for (i in 1:N) {
    
    # set roots of interest
    x_i <- derivatives$roots[i]
    
    # set other roots to use
    j_values <- derivatives$roots[!derivatives$roots %in% x_i]
    
    # get deltas
    delta <- x_i - j_values
    
    for (j in 1:N_int) {
      
      # calculate derivatives for each j (i.e., other roots)
      p_1[j+1] <- delta[j] * p_1[j]
      p_2[j+1] <- delta[j] * p_2[j] + 2 * p_1[j]
      p_3[j+1] <- delta[j] * p_3[j] + 3 * p_2[j]
      
    }
    
    derivatives$p_1[i] <- p_1[N]
    derivatives$p_2[i] <- p_2[N]
    derivatives$p_3[i] <- p_3[N]
    
  }
  
  # define zero matrices
  Ar <- matrix(data = 0, N, N)
  Ar_sym <- matrix(data = 0, N, N)
  Br <- matrix(data = 0, N, N)
  Br_sym <- matrix(data = 0, N, N)
  
  # define A matrix values
  for (j in 1:N) {
    
    for (i in 1:N) {
      
      if(i == j) {
        Ar[i, j] <- 1 / 2 * derivatives$p_2[i] / derivatives$p_1[i]
      } else {
        Ar[i, j] <- 1 / (derivatives$roots[i] - derivatives$roots[j]) * derivatives$p_1[i] / derivatives$p_1[j]
      }
      
      # get symmertic equivalent
      Ar_sym[i, j] <- 2 * sqrt(derivatives$roots[i]) * Ar[i, j]
    }
  }
  
  # define B matrix values
  for (j in 1:N) {
    
    for (i in 1:N) {
      
      if(i == j) {
        Br[i, j] <- 1 / 3 * derivatives$p_3[i] / derivatives$p_1[i]
      } else {
        Br[i, j] <- 2 * Ar[i, j] * (Ar[i, i] - 1 / (derivatives$roots[i] - derivatives$roots[j]))
      }
      
      # get symmertic equivalent
      Br_sym[i, j] <- 4 * derivatives$roots[i] * Br[i, j] + 2 * 3 * Ar[i, j]
    }
  }
  
  # add roots for the symmetric case
  derivatives$roots_sym <- derivatives$roots^(1/2)
  
  # Manuscript formula (adjusted)
  a_weight <- 2
  derivatives$w_i_prime <- 1/(derivatives$roots * derivatives$p_1^2)
  derivatives$W_i_manu <- 1 / (a_weight + 1) * derivatives$w_i_prime * 1 / sum(derivatives$w_i_prime)
  
  B <- Br_sym
  W <- derivatives$W_i_manu
  
  return(list(B, W))
}

ax_colloc <- function(NZ) {
  NZ_int <- NZ - 2 # number of interior points.
  p_list = jacobi.g.recurrences(NZ_int, 1.0, 1.0)  # Shifted Legendre Poly
  m.r <-monic.polynomial.recurrences(p_list)
  roots_Z <- c(0, rev(polynomial.roots(m.r)[[NZ-1]]), 1)
  
  # create a data.frame to store values
  derivatives <- data.frame(
    roots = roots_Z,
    p_1 = rep(0, NZ),
    p_2 = rep(0, NZ),
    p_3 = rep(0, NZ)
  )
  
  # set initial values
  p_1 <- c(1, rep(0, NZ-1))
  p_2 <- rep(0, NZ)
  p_3 <- rep(0, NZ)
  
  for (i in 1:NZ) {
    
    # set roots of interest
    x_i <- derivatives$roots[i]
    
    # set other roots to use
    j_values <- derivatives$roots[!derivatives$roots %in% x_i]
    
    # get deltas
    delta <- x_i - j_values
    
    for (j in 1:(NZ-1)) {
      
      # calculate derivatives for each j (i.e., other roots)
      p_1[j+1] <- delta[j] * p_1[j]
      p_2[j+1] <- delta[j] * p_2[j] + 2 * p_1[j]
      p_3[j+1] <- delta[j] * p_3[j] + 3 * p_2[j]
      
    }
    
    derivatives$p_1[i] <- p_1[NZ]
    derivatives$p_2[i] <- p_2[NZ]
    derivatives$p_3[i] <- p_3[NZ]
    
  }
  
  # define zero matrices
  AZ <- matrix(data = 0, NZ, NZ)
  
  
  # define AZ matrix values
  for (j in 1:NZ) {
    
    for (i in 1:NZ) {
      
      if(i == j) {
        AZ[i, j] <- 1 / 2 * derivatives$p_2[i] / derivatives$p_1[i]
      } else {
        AZ[i, j] <- 1 / (derivatives$roots[i] - derivatives$roots[j]) * derivatives$p_1[i] / derivatives$p_1[j]
      }
    }
  }
  
  return(AZ)
  
}


# Solve function for Shiny App ----
HSDMIX_solve <- function (params, ions, Cin, nt_report){
  NR <- filter(params, name == "nr")$value # numer of grid points along bead radius
  NZ <- filter(params, name == "nz")$value # number of grid points along column axis.
  
  Q <- filter(params, name == "Q")$value # meq/L in resin beads
  L <- filter(params, name == "L")$value # bed depth (cm)
  v <- filter(params, name == "v")$value # superficial flow velocity (cm/s)
  EBED <- filter(params, name == "EBED")$value # bed porosity
  rb <- filter(params, name == "rb")$value # bead radius (cm)
  
  # Ion info
  # Presaturant ion (reference ion A) listed first
  ion_names <- ions$name
  KxA <- ions$KxA
  valence <- ions$valence 
  
  # mass transport paramters
  kL <- ions$kL # film transfer (cm/s)
  Ds <- ions$Ds # surface diffusion (sq. cm/s)
  
  # XXX: Obviously, we will want to load influent concentrations in a more R-idiomatic way. 
  # This is basically Fortran77 :/.
  C_in_t <- data.matrix(Cin)
  
  # Derived parameters ----
  Nt_interp <- dim(C_in_t)[1]
  NION <- length(ion_names)
  LIQUID <- NR + 1 # mnemonic device
  
  C_in_t[, 1] <- C_in_t[, 1] * S_PER_HR # convert time specification from hours to seconds
  
  
  t_max = C_in_t[Nt_interp, 1]
  times <- seq(0.0, t_max*0.99, length.out = nt_report) # seconds
  # times is just a bit short of hours_max to avoid problems with the interpolator.
  
  # XXX: Unfortunately, I can't find  whether deSolve has any way to provide the the timesteps the integrator actually takes
  # so we have to manually define the time scales for the inorganic ions and/or the longer eluting compounds.
  # This is super annoying for troubleshooting BDF or Radau computations 
  # and really inefficient+inconvenient for stiff problems in general.
  
  C_in_0 <- C_in_t[1, 2:(NION+1)] # initial influent concentration (meq/L)
  CT <- sum(C_in_0) # total charge equivalent concentration in feed
  EBCT <- L/v # empty bed contact time.
  tc <- 1.0 # characteristic time # vestigial?
  NEQ <- (NR+1) * NION * NZ
  grid_dims = c((NR+1), NION, NZ)
  
  dv_ions <- valence == 2
  mv_ions <- valence == 1
  mv_ions[1] <- FALSE # exclude presaturant (refrence ion)
  
  # Interpolating functions ----
  # for tracking C_in during integration.
  interp_list <- vector(mode = "list", length = NION)
  for (ii in 1:NION){
    interp_list[[ii]] <- approxfun(C_in_t[ , 1], y = C_in_t[ , ii+1])
  }
  
  # Initialize grid ----
  # Liquid phase is index (NR+1)
  x0 <- array(0.0, grid_dims)
  x0[LIQUID, , 1] <- C_in_0 # set inlet concentrations
  x0[LIQUID, 1, 2:NZ] <- CT  # Rest of liquid in column is full of presaturant
  x0[1:NR, 1, ] <- Q # resin intially loaded with presaturant
  dim(x0) <- c(NEQ)
  
  # collocation ----
  colloc <- rad_colloc(NR)
  BR <- colloc[[1]]  # 1-d radial Laplacian
  WR <- colloc[[2]]  # Gauss-Radau quadrature weights
  AZ <- ax_colloc(NZ) # 1st derivative along Z
  
  
  # Derivative function ----
  diffun <- function(t, x, parms){
    
    dim(x) <- grid_dims
    C <- x[LIQUID, , ]
    q <- x[1:NR, , ]
    qs <- x[NR, , ]
    
    CT_test <- colSums(C)
    
    # update influent concentrations
    for (ii in 1:NION){
      C[ii, 1] <- interp_list[[ii]](t)
    }
    
    # advection collocation intermediate step
    AZ_C <- array(0.0, c(NION, NZ))
    for (ii in 1:NION) {
      AZ_C[ii, ] <- AZ%*%C[ii, ]
    }
    
    
    dx_dt <- array(0.0, grid_dims)
    
    C_star <- array(0.0, c(NION, NZ))  
    if (2 %in% valence){
      # divalent isotherm
      for (ii in 2:NZ){
        cc <- -CT_test[ii] 
        bb <- 1 + (1/qs[1, ii]) * sum(qs[mv_ions, ii]/KxA[mv_ions])
        aa <- (1/qs[1,ii]**2) * qs[dv_ions, ii] / KxA[dv_ions]
        denom <- -bb - sqrt(bb**2 - 4 * aa * cc)
        C_star[1, ii] <- 2 * cc / denom
      }
      
      for (ii in 2:NION){
        C_star[ii, 2:NZ] <- qs[ii, 2:NZ]/KxA[ii]*(C_star[1, 2:NZ]/qs[1, 2:NZ])**valence[ii]
      }
      
      
    } else {
      # monovalent isotherm
      sum_terms <- array(0.0, c(NZ))
      
      for (ii in 2:NZ) {
        sum_terms[ii] <- sum(q[NR, ,ii] / KxA) / CT_test[ii]
      }
      
      for (ii in 2:NION) {
        C_star[ii, 2:NZ] <- q[NR, ii, 2:NZ] / KxA[ii] / sum_terms[2:NZ]
      }
    }
    
    
    J <- array(0.0, c(NION, NZ))  
    for (ii in 2:NION) {
      J[ii , 2:NZ] <- -kL[ii] * (C[ii , 2:NZ] - C_star[ii , 2:NZ])
    }  
    # surface flux calculation
    J[1, 2:NZ] <- - colSums(J[2:NION, 2:NZ]) # Implicitly calculate reference ion
    
    Jas <- 3 / rb * J
    
    dx_dt[LIQUID, , 2:NZ] <- (- v / L * AZ_C[ ,2:NZ] + (1 - EBED) * Jas[ ,2:NZ]) / EBED * tc 
    
    
    # internal diffusion (XXX: loops computationally slow)
    BR_q <- array(0.0, c(NR, NION, NZ))
    
    for (ii in 1:NION){
      for (jj in 2:NZ){
        BR_q[ , ii, jj] <- BR%*%q[ , ii, jj]
      }
    }
    
    dq_dt <- array(0.0, c(NR, NION, NZ))
    for (ii in 2:NION){
      dq_dt[ , ii, ] <- Ds[ii] * tc / rb**2 * BR_q[ , ii, ]
    }
    
    #  dq_dt[ , 1, 2:NZ] <- -rowSums(dq_dt[ , 2:NION, 2:NZ]) # Implicitly calculate reference ion
    # XXX: Why doesn't the above line work? It's not mathematically equivalent to the loop below?
    for (ii in 1:(NR-1)){
      dq_dt[ii, 1, 2:NZ] <- -colSums(dq_dt[ii, 2:NION, 2:NZ])
    }
    
    surf_term <- array(0.0, c(NION, NZ))
    for (ii in 1:NION){
      for (jj in 2:NZ){
        surf_term[ii, jj] <- WR[1:(NR-1)]%*%dq_dt[1:(NR-1), ii, jj]
      }
    }
    
    dx_dt[NR, , 2:NZ] <- (-tc / rb * J[ , 2:NZ] - surf_term[ , 2:NZ])/WR[NR]
    dx_dt[1:(NR-1), , 2:NZ] <- dq_dt[1:(NR-1), , 2:NZ]
    
    list(dx_dt) # return derivatives
  }
  
  # Integration ----
  out <- ode(y = x0, times = times, func = diffun, parms = NULL, method = "bdf")
  # XXX: is there something we can do with diagnose(out) ?
  
  t_out = out[ , 1]/60/60 # hours
  x_out = out[ , 2:(NEQ+1)]
  dim(x_out) <- c(nt_report, (NR+1), NION, NZ)
  
  # Check charge blances at outlet at end of simulation XXX: Maybe move inside of HSDMIX?
  stopifnot(all.equal(sum(x_out[nt_report, NR, , NZ]), Q))
  stopifnot(all.equal(sum(x_out[nt_report, (NR-1), , NZ]), Q))
  #stopifnot(all.equal(sum(x_out[nt_report, LIQUID, , NZ]), CT)) # XXX: TODO: tricky for timevarying infl.
  
  return(list(t_out, x_out)) # TODO: Name these and also provide success/fail info
}



# Run Simulation ----
out <- HSDMIX_solve(params, ions, Cin, nt_report)

# find outlet indices
outlet_id <- dim(out[[2]])[4]
liquid_id <- dim(out[[2]])[2]

# Plot reference ion outlet ----
# TODO: code for picking which curves to plot?
dat = data.frame(hours = out[[1]], conc = out[[2]][, liquid_id, 1, outlet_id])

ggplot(dat, aes(x = hours, y = conc)) +
  geom_point(shape = 1)