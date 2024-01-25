library(readxl)
library(shiny)
library(shinythemes)
library(deSolve)
library(orthopolynom)
library(plotly)
library(shinyjs)
library(tidyr)
library(DataEditR)
library(colorBlindness)
#library(xlsx)
library(ggplot2)


#------------------------------------------------------------------------------#
#~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*#

#HSDMIX_Solve is the engine that runs the Ion Exchange Model. The user inputs
#a parameter data frame (params), an ion data frame (ions), concentration
#data frame (Cin), and an input time which takes the form of a selectInput later
#In the code. nt_report for now is 201, meaning that there are 201 data points
#for each chemical, although this may be changed later to be an input.
#------------------------------------------------------------------------------#




#~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*#
#------------------------------------------------------------------------------#
#unit conversions
#------------------------------------------------------------------------------#
## length
m2cm<-100                               #meters to centimeters
mm2cm<-0.1                              #millimeters to centimeters
cm2cm<-1                                #centimeters to centimeters (for consistency)
in2cm<-2.54                             #inches to centimeters
ft2cm<-12 * in2cm                       #centimeters to feet
## time
sec2sec<-1
min2sec<-60
S_PER_HR <- 60 * 60                     # seconds per hour  #used by HSDMIX function
hour2sec<-60 * min2sec
day2sec<-24 * hour2sec
month2sec<-30 * day2sec                 #assumes 30 day month
year2sec<-365.25 * day2sec
## velocity
mpmin2cmps<-m2cm/min2sec                #meters per minute to centimeters per second
ftpmin2cmps<-ft2cm/min2sec              #feet per minute to centimeters per second
mph2cmps<-m2cm/hour2sec                 #meters per hour to centimeters per second
mmin2cms<-m2cm/min2sec
ftmin2cms<-ft2cm/min2sec
gal2ft3<-0.133680555556
gpmpft2cmps<-gal2ft3 * ft2cm / min2sec  #gallons per minute per foot squared
ft2ps2cm2ps<-(ft2cm)^2                  #feet squared per second to centimeters squared per second
m2ps2cm2ps<-(m2cm)^2                    #meters per second squared to centimeters per second squared
in2ps2cm2ps<-(in2cm)^2                  #inches per second squared to centimeters per second squared
ft2pm2cm2ps<-(ft2cm)^2 / (min2sec)      #feet per minute squared to centimeters per second squared
m2min2cm2s<-(m2cm^2) / (min2sec) 
## volume
gal2ml<-3785.411784
mgd2mlps<-1e6 * gal2ml/day2sec          #mgd to ml/sec
l2ml <- 1000.
k1<-10^(-6.42)
k2<-10^(-10.43)

#~~~~~~~~~~~~~~~~~~~ end unit conversions

#------------------------------------------------------------------------------#
#conversion dictionaries
#------------------------------------------------------------------------------#
##set up dictionaries   ### IF new values are added to drop-downs, must also be added here
length_conv <- c("m"=m2cm, "cm"=cm2cm, "mm"=mm2cm, "in"=in2cm, "ft"=ft2cm)
velocity_conv <- c("cm/s"=cm2cm, "m/s"=m2cm, "m/min"=mpmin2cmps, "m/h"=mph2cmps,
                   "m/hr"=mph2cmps, "in/s"=in2cm, "ft/s"=ft2cm, "ft/min"=ftpmin2cmps,
                   "gpm/ft^2"=gpmpft2cmps)
volumetric_conv <- c("cm^3/s"=cm2cm, "m^3/s"=m2cm^3, "ft^3/s"=ft2cm^3, 
                     "mL/s"=cm2cm, "L/min"=l2ml/min2sec, "mL/min"=1/min2sec,
                     "gpm"=gal2ml/min2sec, "mgd"=mgd2mlps)
time_conv <- c("Hours"=hour2sec, "Days"=day2sec, "Months"=month2sec, "Years"=year2sec,
               "hr"=hour2sec, "day"=day2sec, "month"=month2sec, "year"=year2sec)

kL_conv <- c("ft/s"=ft2cm, "m/s"=m2cm, "cm/s"=cm2cm, "in/s"=in2cm, 
             "m/min"=mpmin2cmps, "ft/min"=ftpmin2cmps, "m/h"=mph2cmps,
             "m/hr"=mph2cmps)
ds_conv <- c("ft^2/s"=ft2ps2cm2ps, "m^2/s"=m2ps2cm2ps, "cm^2/s"=cm2cm,
             "in^2/s"=in2ps2cm2ps)

mass_conv <- c("meq"=1, "meq/L"=1, "mg"=1, "ug"=1e-3, "ng"=1e-6, "mg/L"=1, "ug/L"=1e-3, "ng/L"=1e-6) ### changed

lengthvector<-c("cm", "m", "mm", "in", "ft")
velocityvector<-c("cm/s", "m/s", "m/min", "m/h", "in/s","ft/s","ft/min", "gpm/ft^2")
timevector <- c("hr","day")
flowratevector<-c("cm^3/s", "m^3/s", "ft^3/s", "mL/s", "L/min", "mL/min", "gpm", "mgd")
diametervector<-c("cm", "m", "mm", "in", "ft")
modelvector<-c("HSDM", "PSDM")


#------------------------------------------------------------------------------#
#HSDMIX Function
#------------------------------------------------------------------------------#

# Inputs ----
nt_report = 201 # number of reporting steps

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
HSDMIX_solve <- function (params, ions, Cin, inputtime, nt_report){
  
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
  
  C_in_t[, 1] <- C_in_t[, 1] * inputtime # convert time specification from hours to seconds
  
  
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

#------------------------------------------------------------------------------#
#PSDMIX Function
#------------------------------------------------------------------------------#

PSDMIX_solve <- function (params, ions, Cin, inputtime, nt_report){
  NR <- filter(params, name == "nr")$value # numer of grid points along bead radius
  NZ <- filter(params, name == "nz")$value # number of grid points along column axis.
  
  Q <- filter(params, name == "Q")$value # meq/L in resin beads
  L <- filter(params, name == "L")$value # bed depth (cm)
  v <- filter(params, name == "v")$value # superficial flow velocity (cm/s)
  EBED <- filter(params, name == "EBED")$value # bed porosity
  EPOR <- filter(params, name == "EPOR")$value # pellet porosity
  rb <- filter(params, name == "rb")$value # bead radius (cm)
  
  # Ion info
  # Presaturant ion (reference ion A) listed first
  ion_names <- ions$name
  KxA <- ions$KxA
  valence <- ions$valence 
  
  # mass transport paramters
  kL <- ions$kL # film transfer (cm/s)
  Ds <- ions$Ds # surface diffusion (sq. cm/s)
  Dp <- ions$Dp # pore diffusion (sq. cm/s)
  
  # XXX: Obviously, we will want to load influent concentrations in a more R-idiomatic way. 
  # This is basically Fortran77 :/.
  C_in_t <- data.matrix(Cin)
  
  # Derived parameters ----
  Nt_interp <- dim(C_in_t)[1]
  NION <- length(ion_names)
  LIQUID <- NR + 1 # mnemonic device
  
  C_in_t[, 1] <- C_in_t[, 1] * inputtime # convert time specification from hours to seconds
  
  
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
    Y <- x[1:NR, , ]
    q <- Y / (1 - EPOR)
    
    
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
    
    Cpore <- array(0.0, c(NR, NION, NZ))
    
    
    if (2 %in% valence){
      # divalent isotherm
      for (jj in 1:NR){
        for (ii in 2:NZ){
          cc <- -CT_test[ii] 
          bb <- 1 + (1/q[jj, 1, ii]) * sum(q[jj, mv_ions, ii]/KxA[mv_ions])
          aa <- (1/q[jj,1,ii]**2) * q[jj,dv_ions, ii] / KxA[dv_ions]
          denom <- -bb - sqrt(bb**2 - 4 * aa * cc)
          Cpore[jj,1, ii] <- 2 * cc / denom
        }
        
        for (ii in 2:NION){
          Cpore[jj, ii, 2:NZ] <- q[jj, ii, 2:NZ]/KxA[ii]*(Cpore[jj, 1, 2:NZ]/q[jj, 1, 2:NZ])**valence[ii]
        }
      }
      
      
    } else {
      # monovalent isotherm  
      sum_terms <- array(0.0, c(NZ))
      
      for (jj in 1:NR){
        for (ii in 2:NZ) {
          sum_terms[ii] <- sum(q[jj, ,ii] / KxA) / CT_test[ii]
        }
        
        for (ii in 2:NION) {
          Cpore[jj, ii, 2:NZ] <- q[jj, ii, 2:NZ] / KxA[ii] / sum_terms[2:NZ]
        }
      }
    }
    
    C_star <- Cpore[NR, , ]  
    
    
    J <- array(0.0, c(NION, NZ))  
    for (ii in 2:NION) {
      J[ii , 2:NZ] <- -kL[ii] * (C[ii , 2:NZ] - C_star[ii , 2:NZ])
    }  
    # surface flux calculation
    J[1, 2:NZ] <- - colSums(J[2:NION, 2:NZ]) # Implicitly calculate reference ion
    
    Jas <- 3 / rb * J
    
    dx_dt[LIQUID, , 2:NZ] <- (- v / L * AZ_C[ ,2:NZ] + (1 - EBED) * Jas[ ,2:NZ]) / EBED * tc 
    
    
    # internal diffusion (XXX: loops computationally slow)
    BR_Y <- array(0.0, c(NR, NION, NZ))
    BR_Cpore <- array(0.0, c(NR, NION, NZ))
    for (ii in 1:NION){
      for (jj in 2:NZ){
        BR_Y[ , ii, jj] <- BR%*%Y[ , ii, jj]
        BR_Cpore[ , ii, jj] <- BR%*%Cpore[ , ii, jj]
      }
    }
    
    
    dY_dt <- array(0.0, c(NR, NION, NZ))
    for (ii in 2:NION){
      dY_dt[ , ii, ] <- tc * (EPOR * (Dp[ii] - Ds[ii]) * BR_Cpore[ , ii, ] + Ds[ii] * BR_Y[ , ii, ]) / rb**2
    }
    
    for (ii in 1:(NR-1)){
      dY_dt[ii, 1, 2:NZ] <- -colSums(dY_dt[ii, 2:NION, 2:NZ]) # Implicitly calculate reference ion
    }
    
    surf_term <- array(0.0, c(NION, NZ))
    for (ii in 1:NION){
      for (jj in 2:NZ){
        surf_term[ii, jj] <- WR[1:(NR-1)]%*%dY_dt[1:(NR-1), ii, jj]
      }
    }
    
    dx_dt[NR, , 2:NZ] <- (-tc / rb * J[ , 2:NZ] - surf_term[ , 2:NZ])/WR[NR]
    dx_dt[1:(NR-1), , 2:NZ] <- dY_dt[1:(NR-1), , 2:NZ]
    
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

#~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*#


#------------------------------------------------------------------------------#
#FUNCTION LIBRARY
#------------------------------------------------------------------------------#


#~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*#
#------------------------------------------------------------------------------#
#process_files
# process_files is a function that reads in the input file, then splits each 
# page into 3 csv files and saves them separate. This may seem unnecessary at
# first, but the reason we do this is because the package DataEditR can only 
# read csv files and not xlsx files. We wanted to have a data frame that is as 
# easily editable as an excel page, and this is the only function that exists to
# do that
#------------------------------------------------------------------------------#

process_files <- function (file) {
  
  effluent<-data.frame(time=(0), CHLORIDE=(0))
  
  params<-read_xlsx(file, sheet="params")
  ions<-read_xlsx(file, sheet="ions")
  cin<-read_xlsx(file, sheet="Cin")
  
  write.csv(params, "paramsheet.csv", row.names=FALSE)
  write.csv(ions, "ionsheet.csv", row.names=FALSE)
  write.csv(cin, "cinsheet.csv", row.names=FALSE)
  
  #Checks for effluent data, if unavailable use empty dataset
  #It is very likely that the user will not have the data available and 
  #Has no intention of plotting it, which makes tryCatch useful ehre
  
  tryCatch({
    
    eff<-read_xlsx(file, sheet='effluent')
    write.csv(eff, "effluent.csv", row.names=FALSE)
    
  },
  
  warning=function(war){
    #pass
  },
  error=function(err){
    
    write.csv(effluent, "effluent.csv", row.names=FALSE)
    
  })
  
}


#------------------------------------------------------------------------------#
#Create Plotly
#This function is used to create the plot for just the ions. 
#frame1 = Computed Data
#frame2 = Effluent Data
#frame3 = Influent Data
#------------------------------------------------------------------------------#

create_plotly<-function(frame1, frame2, frame3){
  
  
  #Create a subset of data that 
  counterionframe<-subset(frame1, name %in% c("CHLORIDE", "SULFATE", "NITRATE", "BICARBONATE"))
  counterioneff<-subset(frame2, name %in% c("CHLORIDE_effluent", "SULFATE_effluent", "NITRATE_effluent", "BICARBONATE_effluent"))
  counterioninfluent<-frame3
  
  
  # counterionfig2<-ggplot(counterionframe, aes(hours, conc, color=name))+
  #   geom_line(data=counterionframe)+geom_point(data=counterioneff)+geom_line(data=counterioninfluent)+geom_point(data=counterioninfluent)
  # 
  # counterionfig<-ggplotly(counterionfig2)
  
  #Using the curated data, plot
  counterionfig<-plot_ly(counterionframe, x=~hours, y=~conc, type='scatter', mode='lines', color=~name, colors=SteppedSequential5Steps)%>%
    add_trace(data=counterioneff, x=~hours, y=~conc, mode='markers')%>%
    add_trace(data=counterioninfluent, x=~hours, y=~conc, mode='lines+markers')
  
  return(counterionfig)
  
}



#Same thing as create_plotly but just for the ions
create_plotly2<-function(frame1, frame2, frame3){
  
  ionframe<-subset(frame1, !(name %in% c("CHLORIDE", "SULFATE", "NITRATE", "BICARBONATE")))
  ioneff<-subset(frame2, !(name %in% c("CHLORIDE_effluent", "SULFATE_effluent", "NITRATE_effluent", "BICARBONATE_effluent")))
  ioninfluent<-frame3
  
  ionfig<-plot_ly(ionframe, x=~hours, y=~conc, type='scatter', mode='lines', color=~name, colors=SteppedSequential5Steps)%>%
    add_trace(data=ioneff, x=~hours, y=~conc, mode='markers')%>%
    add_trace(data=ioninfluent, x=~hours, y=~conc, mode='lines+markers')
  
  return(ionfig)
  
}


#------------------------------------------------------------------------------#
#Bed Volume
#------------------------------------------------------------------------------#

get_bv_in_sec <- function(input) {
  #get number of seconds per bv
  if (input$veloselect == 'Linear') {
    Vv = input$Vv*velocity_conv[input$VelocityUnits]
  } else {
    Vv = input$Fv * volumetric_conv[input$FlowrateUnits]/(pi/4 * ((input$Dv * length_conv[input$DiameterUnits])**2))
  }
  
  # print(input)
  
  ## divide converted length by velocity to get BV in seconds
  return(input$Lv*length_conv[input$LengthUnits]/Vv)
  
}



#The Cin data that gets plotted may or may not be in the right units
#So this function corrects that based on a column in the ions tab 
#All units in this tab are in meq
cin_correct<-function(ions, cins){
  
  corr_cin <- cins
  
  ## we now know that compounds are in both lists, assuming error == 0
  
  for (item in 1:nrow(ions)) {
    ## convert mass units
    mass_mult <- 1.
    mass_units <- ions[item, "conc_units"]  # convenience variable
    if (mass_units != 'meq') {
      mass_mult <- mass_conv[mass_units] / (ions[item, "mw"]) * ions[item, "valence"]  ## TODO: check this math
    }
    
    #should multiply a column by conversion factor
    compound <- ions[item, "name"]  # convenience variable 
    
    corr_cin[, compound] <- cins[, compound] * mass_mult
  }
  
  
  return(corr_cin)
  
  
  
}



effluent_data_processor<-function(ion, effluent){
  if(nrow(effluent)>1){                                      #If effluent data is not empty
    
    mydata<-mass_converter_mgl(ion, effluent)                              #convert to mgl
    colnames(mydata)<-paste(colnames(mydata), "effluent", sep="_")#Distinguish the names from the simulated data
    
    timevec<-data.frame(time=c(mydata[,1]))
    concframe<-gather(mydata[,2:ncol(mydata)])                    #Gather into shape that is easy to convert and plot
    effframe<-cbind(timevec, concframe)
    
    colnames(effframe)<-c("time", "name", "conc")
    
    return(effframe)
    
  }
  
  else{
    
    effframe<-data.frame(time=NA, name=NA, conc=NA)
    
    return(effframe)
    
  }
}


influent_chemical_renamer<-function(influent, influent_hours){
  cindata<-influent[,2:ncol(influent_hours)]
  time<-influent[,1]
  
  colnames(cindata)<-paste(colnames(cindata), "influent", sep="_")
  alldat<-cbind(time, cindata)
  
  return(alldat)
}


influent_organizer<-function(influent, influent_hours){
  
  cindat_organized<-tidyr::gather(influent[2:ncol(influent_hours)])
  cin_time<-influent[,1]
  cin_prepped<-cbind(cin_time, cindat_organized)
  
  colnames(cin_prepped)<-c("hours", "name", "conc")
  
  return(cin_prepped)
  
}


#------------------------------------------------------------------------------#
#HSDMIX Prep
#This function makes sure that the appropriate data frames are created
#and that they have the converted values 
#------------------------------------------------------------------------------#

HSDMIX_prep <- function (input, iondata, concdata, nt_report) {
  ## prepare paramdataframe for use by HSDMIX_solve
  if (input$veloselect == 'Linear') {
    Vv = input$Vv*velocity_conv[input$VelocityUnits]
  } else {
    Vv = input$Fv * volumetric_conv[input$FlowrateUnits]/(pi/4 * ((input$Dv * length_conv[input$DiameterUnits])**2))
  }
  
  paramdataframe <- data.frame(
    name=c("Q", "EBED", "L", "v", "rb", "kL", "Ds", "nr", "nz", "time"),
    value=c(input$Qv,
            input$EBEDv,
            input$Lv*length_conv[input$LengthUnits],
            Vv,
            input$rbv*length_conv[input$rbunits], NA, NA,
            input$nrv,
            input$nzv, 1),
    units=c("meq/L", NA, "cm", "cm/s", "cm", NA, NA, NA, NA, input$timeunits)
  )
  
  
  ## check that ions and concdata match
  error <- 0
  for (item in 1:nrow(iondata)) {
    ## checking ions are in concentration 
    if (!(iondata[item, 'name'] %in% colnames(concdata))) {
      print(paste0(iondata[item, 'name'], " not found in Concentration Data Columns"))
      error <- error + 1
    }
  }
  
  for (item in colnames(concdata)) {
    ## checking concentration ions in ion list
    if (item != "time") {
      if (!(item %in% iondata[, 'name'])) {
        print(paste0(item, " not found in Ion Data"))
        error <- error + 1
      }
    }
  } 
  
  print(paste0("Number of errors: ", error))
  
  ## replicate inputs, will change in next section
  corr_ions <- iondata
  corr_cin <- concdata
  
  ## we now know that compounds are in both lists, assuming error == 0
  if (error == 0) {
    corr_cin <- cin_correct(iondata, concdata)
    for (item in 1:nrow(iondata)) {
      
      ### TODO: Need to check the mass transfer unit conversions
      ## convert kL to cm/s
      corr_ions[item, 'kL'] <- iondata[item, 'kL'] * kL_conv[iondata[item, 'kL_units']]
      
      ## convert Ds to cm/s^2
      corr_ions[item, 'Ds'] <- iondata[item, 'Ds'] * ds_conv[iondata[item, 'Ds_units']]
      
    }
  }
  
  
  timeconverter <- time_conv[input$timeunits2]  ### TODO: Is this really necessary, or doing what we think it is doing? 
  
  
  if (error == 0) {
    return (HSDMIX_solve(paramdataframe, corr_ions, corr_cin, timeconverter, nt_report))
  } else {
    return (error)
  }
}



PSDMIX_prep<-function(input, iondata, concdata, nt_report){
  # prepare paramdataframe for use by HSDMIX_solve
  if (input$veloselect == 'Linear') {
    Vv = input$Vv*velocity_conv[input$VelocityUnits]
  } else {
    Vv = input$Fv * volumetric_conv[input$FlowrateUnits]/(pi/4 * ((input$Dv * length_conv[input$DiameterUnits])**2))
  }
  
  
  paramdataframe<-data.frame(name=c("Q", "EBED", "EPOR", "L", "v", "rb", "kL", "Ds", "nr", "nz", "time"),
                             value=c(input$Qv,
                                     input$EBEDv,
                                     input$EPORv,
                                     input$Lv*length_conv[input$LengthUnits],
                                     Vv,
                                     input$rbv*length_conv[input$rbunits], NA, NA,
                                     input$nrv,
                                     input$nzv, 1),
                             units=c("meq/L", NA, "cm", "cm/s", "cm", NA, NA, NA, NA, NA, input$timeunits)
  )
  
  
  # check that ions and concdata match
  error <- 0
  for (item in 1:nrow(iondata)) {
    ## checking ions are in concentration
    if (!(iondata[item, 'name'] %in% colnames(concdata))) {
      print(paste0(iondata[item, 'name'], " not found in Concentration Data Columns"))
      error <- error + 1
    }
  }
  
  for (item in colnames(concdata)) {
    ## checking concentration ions in ion list
    if (item != "time") {
      if (!(item %in% iondata[, 'name'])) {
        print(paste0(item, " not found in Ion Data"))
        error <- error + 1
      }
    }
  }
  
  print(paste0("Number of errors: ", error))
  
  ## replicate inputs, will change in next section
  corr_ions <- iondata
  corr_cin <- concdata
  
  ## we now know that compounds are in both lists, assuming error == 0
  if (error == 0) {
    corr_cin <- cin_correct(iondata, concdata)
    for (item in 1:nrow(iondata)) {
      
      ### TODO: Need to check the mass transfer unit conversions
      ## convert kL to cm/s
      corr_ions[item, 'kL'] <- iondata[item, 'kL'] * kL_conv[iondata[item, 'kL_units']]
      
      ## convert Ds to cm/s^2
      corr_ions[item, 'Ds'] <- iondata[item, 'Ds'] * ds_conv[iondata[item, 'Ds_units']]
      
      # Dp units are the same as Ds
      corr_ions[item, 'Dp'] <- iondata[item, 'Dp'] * ds_conv[iondata[item, 'Dp_units']]
      
    }
  }
  
  
  timeconverter <- time_conv[input$timeunits2]  ### TODO: Is this really necessary, or doing what we think it is doing?
  
  
  if (error == 0) {
    return (PSDMIX_solve(paramdataframe, corr_ions, corr_cin, timeconverter, nt_report))
  } else {
    return (error)
  }
}



HSDMIX_in_hours_mgl<-function(HSDMIXoutput, ions, time){
  
  allchems_meq<-HSDMIXoutput
  iondata<-ions
  timedata<-time
  
  massvector_meq<-c(iondata$mw/iondata$valence)
  
  correctedchems<-mapply('*', allchems_meq, massvector_meq)
  correctedchemsframe<-data.frame(correctedchems)
  allchem<-cbind(time, correctedchemsframe)
  allchemgathered<-tidyr::gather(correctedchemsframe)
  nameandconcs<-data.frame(name=allchemgathered[,1],
                           conc=allchemgathered[,2])
  
  allchemicalsmgl<-cbind(timedata, nameandconcs)
  
  return(allchemicalsmgl)
  
}



HSDMIX_cc0<-function(HSDMIXoutput, c0values){
  
  allchemicalscc0<-mapply('/', HSDMIXoutput, c0values)
  dataframecc0<-data.frame(allchemicalscc0)
  organized<-tidyr::gather(dataframecc0)
  cc0frame<-data.frame(name=organized[,1],
                       conc=organized[,2])
  
  return(cc0frame)
  
}


effluent_cc0<-function(effluent, ions, c0values){
  
  if(nrow(effluent)>1){
    effdat<-cin_correct(ions, effluent)
    cc0dat<-c0values
    
    time<-effdat[,1]
    subseteffdat<-effdat[,2:ncol(effdat)]
    
    effcc0<-mapply('/', subseteffdat, cc0dat)
    effcc02<-data.frame(effcc0)
    effcc03<-tidyr::gather(effcc02)
    effcc04<-data.frame(name=effcc03[,1],
                        conc=effcc03[,2])
    return(effcc04) ## use return() for clarity
    
  }
  
  else{
    effcc04<-data.frame(hours=c(NA), name=c(NA), conc=c(NA))
    return(effcc04)
  }
  
}


#------------------------------------------------------------------------------#
#cc0 Conversion function
#------------------------------------------------------------------------------#


cc0_conv_meq <- function (iondata, concdata) {
  error <- 0
  for (item in 1:nrow(iondata)) {
    ## checking ions are in concentration 
    if (!(iondata[item, 'name'] %in% colnames(concdata))) {
      print(paste0(iondata[item, 'name'], " not found in Concentration Data Columns"))
      error <- error + 1
    }
  }
  
  for (item in colnames(concdata)) {
    ## checking concentration ions in ion list
    if (item != "time") {
      if (!(item %in% iondata[, 'name'])) {
        print(paste0(item, " not found in Ion Data"))
        error <- error + 1
      }
    }
  } 
  if (error == 0) {
    cc0 <- c()
    meq_conv <- c()
    for (item in 1:nrow(iondata)) {
      if (iondata[item, "conc_units"] == "meq" || iondata[item, "conc_units"] == "meq/L") {
        meq_conv[iondata[item, "name"]] <- 1
      }
      else{
        meq_conv[iondata[item, "name"]] <- mass_conv[iondata[item, "conc_units"]] * iondata[item,"valence"] / iondata[item, "mw"]
      }
    }
    for (item in colnames(concdata)) {
      ## checking concentration ions in ion list
      if (item != "time") {
        cc0[item] <- concdata[1, item] * meq_conv[item]
      }
    }
    #print(cc0)
    return(cc0)
  } else{
    return(concdata[2,2:ncol(concdata)])
  }
  
} 


mass_converter_mgl <- function (iondata, concs) {
  corr_cin <- concs
  for (item in 1:nrow(iondata)){
    mass_mult <- 1.
    mass_units <- iondata[item, "conc_units"]  # convenience variable
    if (mass_units == 'meq' || mass_units == 'meq/L') {
      mass_mult <- (iondata[item, "mw"])/iondata[item, "valence"] ## TODO: check this math
    }
    else { 
      mass_mult <- mass_conv[iondata[item,'conc_units']]
      
    }
    
    #should multiply a column by conversion factor
    compound <- iondata[item, "name"]  # convenience variable 
    
    corr_cin[, compound] <- concs[, compound] * mass_mult
  }
  
  #print(corr_cin)
  return(corr_cin)
  
}




#~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*#





#~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*#


wd <- getwd()
process_files(paste0("config.xlsx"))
#process_effluent(paste0("config.xlsx"))


#~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*#
#------------------------------------------------------------------------------#
#BEGIN UI#
# The UI section deals with all the visuals that the user will interact with
# This section is split up into 3 tabs, Input, Output, and About
# The Input section has 2 tabs, parameters and ions, these have default values
# So that the code can be ran without inputting anything, but the user can upload
# an excel file that will overwrite this data, and the user can edit it further
# within these tabs
# The Output section displays the plots where the x and y-axis can be converted
# into other units
# The about section contains information about the tool for further details
# The side bar panel which is present throughout the app is used for minor 
# adjustments as well as the axis conversions and file inputs
#------------------------------------------------------------------------------#

ui <- fluidPage(
  
  useShinyjs(),
  # 
  # tags$html(class = "no-js", lang="en"),
  # tags$head(
  #   tags$link(rel="stylesheet", media="all", href="https://www.epa.gov/themes/epa_theme/css/styles.css?r6lsex")),
  
  # Site Header
  #This HTML code is a speical theme that makes it look like an EPA app
  
  # HTML(
  #   '
  #   <div>
  #     <div class="js-view-dom-id-epa-alerts--public">
  #       <noscript>
  #         <div class="usa-site-alert usa-site-alert--info">
  #           <div class="usa-alert">
  #             <div class="usa-alert__body">
  #               <div class="usa-alert__text">
  #                 <p>JavaScript appears to be disabled on this computer. Please <a href="/alerts">click here to see any active alerts</a>.</p>
  #               </div>
  #             </div>
  #           </div>
  #         </div>
  #       </noscript>
  #     </div>
  #   </div>
  #   <header class="l-header">
  #     <div class="usa-overlay"></div>
  #     <div class="l-constrain">
  #       <div class="l-header__navbar">
  #         <div class="l-header__branding">
  #           <a class="site-logo" href="/" aria-label="Home" title="Home" rel="home">
  #             <span class="site-logo__image">
  #               <svg class="site-logo__svg" viewBox="0 0 1061 147" aria-hidden="true" xmlns="http://www.w3.org/2000/svg">
  #                 <path d="M112.8 53.5C108 72.1 89.9 86.8 69.9 86.8c-20.1 0-38-14.7-42.9-33.4h.2s9.8 10.3-.2 0c3.1 3.1 6.2 4.4 10.7 4.4s7.7-1.3 10.7-4.4c3.1 3.1 6.3 4.5 10.9 4.4 4.5 0 7.6-1.3 10.7-4.4 3.1 3.1 6.2 4.4 10.7 4.4 4.5 0 7.7-1.3 10.7-4.4 3.1 3.1 6.3 4.5 10.9 4.4 4.3 0 7.4-1.2 10.5-4.3zM113.2 43.5c0-24-19.4-43.5-43.3-43.5-24 0-43.5 19.5-43.5 43.5h39.1c-4.8-1.8-8.1-6.3-8.1-11.6 0-7 5.7-12.5 12.5-12.5 7 0 12.7 5.5 12.7 12.5 0 5.2-3.1 9.6-7.6 11.6h38.2zM72.6 139.3c.7-36.9 29.7-68.8 66.9-70 0 37.2-30 68-66.9 70zM67.1 139.3c-.7-36.9-29.7-68.8-67.1-70 0 37.2 30.2 68 67.1 70zM240 3.1h-87.9v133.1H240v-20.4h-60.3v-36H240v-21h-60.3v-35H240V3.1zM272.8 58.8h27.1c9.1 0 15.2-8.6 15.1-17.7-.1-9-6.1-17.3-15.1-17.3h-25.3v112.4h-27.8V3.1h62.3c20.2 0 35 17.8 35.2 38 .2 20.4-14.8 38.7-35.2 38.7h-36.3v-21zM315.9 136.2h29.7l12.9-35h54.2l-8.1-21.9h-38.4l18.9-50.7 39.2 107.6H454L400.9 3.1h-33.7l-51.3 133.1zM473.3.8v22.4c0 1.9.2 3.3.5 4.3s.7 1.7 1 2.2c1.2 1.4 2.5 2.4 3.9 2.9 1.5.5 2.8.7 4.1.7 2.4 0 4.2-.4 5.5-1.3 1.3-.8 2.2-1.8 2.8-2.9.6-1.1.9-2.3 1-3.4.1-1.1.1-2 .1-2.6V.8h4.7v24c0 .7-.1 1.5-.4 2.4-.3 1.8-1.2 3.6-2.5 5.4-1.8 2.1-3.8 3.5-6 4.2-2.2.6-4 .9-5.3.9-1.8 0-3.8-.3-6.2-1.1-2.4-.8-4.5-2.3-6.2-4.7-.5-.8-1-1.8-1.4-3.2-.4-1.3-.6-3.3-.6-5.9V.8h5zM507.5 14.5v-2.9l4.6.1-.1 4.1c.2-.3.4-.7.8-1.2.3-.5.8-.9 1.4-1.4.6-.5 1.4-.9 2.3-1.3.9-.3 2.1-.5 3.4-.4.6 0 1.4.1 2.4.3.9.2 1.9.6 2.9 1.2s1.8 1.5 2.4 2.6c.6 1.2.9 2.8.9 4.7l-.4 17-4.6-.1.4-16c0-.9 0-1.7-.2-2.4-.1-.7-.5-1.3-1.1-1.9-1.2-1.2-2.6-1.8-4.3-1.8-1.7 0-3.1.5-4.4 1.7-1.3 1.2-2 3.1-2.1 5.7l-.3 14.5-4.5-.1.5-22.4zM537.2.9h5.5V6h-5.5V.9m.5 10.9h4.6v25.1h-4.6V11.8zM547.8 11.7h4.3V6.4l4.5-1.5v6.8h5.4v3.4h-5.4v15.1c0 .3 0 .6.1 1 0 .4.1.7.4 1.1.2.4.5.6 1 .8.4.3 1 .4 1.8.4 1 0 1.7-.1 2.2-.2V37c-.9.2-2.1.3-3.8.3-2.1 0-3.6-.4-4.6-1.2-1-.8-1.5-2.2-1.5-4.2V15.1h-4.3v-3.4zM570.9 25.2c-.1 2.6.5 4.8 1.7 6.5 1.1 1.7 2.9 2.6 5.3 2.6 1.5 0 2.8-.4 3.9-1.3 1-.8 1.6-2.2 1.8-4h4.6c0 .6-.2 1.4-.4 2.3-.3 1-.8 2-1.7 3-.2.3-.6.6-1 1-.5.4-1 .7-1.7 1.1-.7.4-1.5.6-2.4.8-.9.3-2 .4-3.3.4-7.6-.2-11.3-4.5-11.3-12.9 0-2.5.3-4.8 1-6.8s2-3.7 3.8-5.1c1.2-.8 2.4-1.3 3.7-1.6 1.3-.2 2.2-.3 3-.3 2.7 0 4.8.6 6.3 1.6s2.5 2.3 3.1 3.9c.6 1.5 1 3.1 1.1 4.6.1 1.6.1 2.9 0 4h-17.5m12.9-3v-1.1c0-.4 0-.8-.1-1.2-.1-.9-.4-1.7-.8-2.5s-1-1.5-1.8-2c-.9-.5-2-.8-3.4-.8-.8 0-1.5.1-2.3.3-.8.2-1.5.7-2.2 1.3-.7.6-1.2 1.3-1.6 2.3-.4 1-.7 2.2-.8 3.6h13zM612.9.9h4.6V33c0 1 .1 2.3.2 4h-4.6l-.1-4c-.2.3-.4.7-.7 1.2-.3.5-.8 1-1.4 1.5-1 .7-2 1.2-3.1 1.4l-1.5.3c-.5.1-.9.1-1.4.1-.4 0-.8 0-1.3-.1s-1.1-.2-1.7-.3c-1.1-.3-2.3-.9-3.4-1.8s-2.1-2.2-2.9-3.8c-.8-1.7-1.2-3.9-1.2-6.6.1-4.8 1.2-8.3 3.4-10.5 2.1-2.1 4.7-3.2 7.6-3.2 1.3 0 2.4.2 3.4.5.9.3 1.6.7 2.2 1.2.6.4 1 .9 1.3 1.4.3.5.6.8.7 1.1V.9m0 23.1c0-1.9-.2-3.3-.5-4.4-.4-1.1-.8-2-1.4-2.6-.5-.7-1.2-1.3-2-1.8-.9-.5-2-.7-3.3-.7-1.7 0-2.9.5-3.8 1.3-.9.8-1.6 1.9-2 3.1-.4 1.2-.7 2.3-.7 3.4-.1 1.1-.2 1.9-.1 2.4 0 1.1.1 2.2.3 3.4.2 1.1.5 2.2 1 3.1.5 1 1.2 1.7 2 2.3.9.6 2 .9 3.3.9 1.8 0 3.2-.5 4.2-1.4 1-.8 1.7-1.8 2.1-3 .4-1.2.7-2.4.8-3.4.1-1.4.1-2.1.1-2.6zM643.9 26.4c0 .6.1 1.3.3 2.1.1.8.5 1.6 1 2.3.5.8 1.4 1.4 2.5 1.9s2.7.8 4.7.8c1.8 0 3.3-.3 4.4-.8 1.1-.5 1.9-1.1 2.5-1.8.6-.7 1-1.5 1.1-2.2.1-.7.2-1.2.2-1.7 0-1-.2-1.9-.5-2.6-.4-.6-.9-1.2-1.6-1.6-1.4-.8-3.4-1.4-5.9-2-4.9-1.1-8.1-2.2-9.5-3.2-1.4-1-2.3-2.2-2.9-3.5-.6-1.2-.8-2.4-.8-3.6.1-3.7 1.5-6.4 4.2-8.1 2.6-1.7 5.7-2.5 9.1-2.5 1.3 0 2.9.2 4.8.5 1.9.4 3.6 1.4 5 3 .5.5.9 1.1 1.2 1.7.3.5.5 1.1.6 1.6.2 1.1.3 2.1.3 2.9h-5c-.2-2.2-1-3.7-2.4-4.5-1.5-.7-3.1-1.1-4.9-1.1-5.1.1-7.7 2-7.8 5.8 0 1.5.5 2.7 1.6 3.5 1 .8 2.6 1.4 4.7 1.9 4 1 6.7 1.8 8.1 2.2.8.2 1.4.5 1.8.7.5.2 1 .5 1.4.9.8.5 1.4 1.1 1.9 1.8s.8 1.4 1.1 2.1c.3 1.4.5 2.5.5 3.4 0 3.3-1.2 6-3.5 8-2.3 2.1-5.8 3.2-10.3 3.3-1.4 0-3.2-.3-5.4-.8-1-.3-2-.7-3-1.2-.9-.5-1.8-1.2-2.5-2.1-.9-1.4-1.5-2.7-1.7-4.1-.3-1.3-.4-2.4-.3-3.2h5zM670 11.7h4.3V6.4l4.5-1.5v6.8h5.4v3.4h-5.4v15.1c0 .3 0 .6.1 1 0 .4.1.7.4 1.1.2.4.5.6 1 .8.4.3 1 .4 1.8.4 1 0 1.7-.1 2.2-.2V37c-.9.2-2.1.3-3.8.3-2.1 0-3.6-.4-4.6-1.2-1-.8-1.5-2.2-1.5-4.2V15.1H670v-3.4zM705.3 36.9c-.3-1.2-.5-2.5-.4-3.7-.5 1-1.1 1.8-1.7 2.4-.7.6-1.4 1.1-2 1.4-1.4.5-2.7.8-3.7.8-2.8 0-4.9-.8-6.4-2.2-1.5-1.4-2.2-3.1-2.2-5.2 0-1 .2-2.3.8-3.7.6-1.4 1.7-2.6 3.5-3.7 1.4-.7 2.9-1.2 4.5-1.5 1.6-.1 2.9-.2 3.9-.2s2.1 0 3.3.1c.1-2.9-.2-4.8-.9-5.6-.5-.6-1.1-1.1-1.9-1.3-.8-.2-1.6-.4-2.3-.4-1.1 0-2 .2-2.6.5-.7.3-1.2.7-1.5 1.2-.3.5-.5.9-.6 1.4-.1.5-.2.9-.2 1.2h-4.6c.1-.7.2-1.4.4-2.3.2-.8.6-1.6 1.3-2.5.5-.6 1-1 1.7-1.3.6-.3 1.3-.6 2-.8 1.5-.4 2.8-.6 4.2-.6 1.8 0 3.6.3 5.2.9 1.6.6 2.8 1.6 3.4 2.9.4.7.6 1.4.7 2 .1.6.1 1.2.1 1.8l-.2 12c0 1 .1 3.1.4 6.3h-4.2m-.5-12.1c-.7-.1-1.6-.1-2.6-.1h-2.1c-1 .1-2 .3-3 .6s-1.9.8-2.6 1.5c-.8.7-1.2 1.7-1.2 3 0 .4.1.8.2 1.3s.4 1 .8 1.5.9.8 1.6 1.1c.7.3 1.5.5 2.5.5 2.3 0 4.1-.9 5.2-2.7.5-.8.8-1.7 1-2.7.1-.9.2-2.2.2-4zM714.5 11.7h4.3V6.4l4.5-1.5v6.8h5.4v3.4h-5.4v15.1c0 .3 0 .6.1 1 0 .4.1.7.4 1.1.2.4.5.6 1 .8.4.3 1 .4 1.8.4 1 0 1.7-.1 2.2-.2V37c-.9.2-2.1.3-3.8.3-2.1 0-3.6-.4-4.6-1.2-1-.8-1.5-2.2-1.5-4.2V15.1h-4.3v-3.4zM737.6 25.2c-.1 2.6.5 4.8 1.7 6.5 1.1 1.7 2.9 2.6 5.3 2.6 1.5 0 2.8-.4 3.9-1.3 1-.8 1.6-2.2 1.8-4h4.6c0 .6-.2 1.4-.4 2.3-.3 1-.8 2-1.7 3-.2.3-.6.6-1 1-.5.4-1 .7-1.7 1.1-.7.4-1.5.6-2.4.8-.9.3-2 .4-3.3.4-7.6-.2-11.3-4.5-11.3-12.9 0-2.5.3-4.8 1-6.8s2-3.7 3.8-5.1c1.2-.8 2.4-1.3 3.7-1.6 1.3-.2 2.2-.3 3-.3 2.7 0 4.8.6 6.3 1.6s2.5 2.3 3.1 3.9c.6 1.5 1 3.1 1.1 4.6.1 1.6.1 2.9 0 4h-17.5m12.9-3v-1.1c0-.4 0-.8-.1-1.2-.1-.9-.4-1.7-.8-2.5s-1-1.5-1.8-2c-.9-.5-2-.8-3.4-.8-.8 0-1.5.1-2.3.3-.8.2-1.5.7-2.2 1.3-.7.6-1.2 1.3-1.6 2.3-.4 1-.7 2.2-.8 3.6h13zM765.3 29.5c0 .5.1 1 .2 1.4.1.5.4 1 .8 1.5s.9.8 1.6 1.1c.7.3 1.6.5 2.7.5 1 0 1.8-.1 2.5-.3.7-.2 1.3-.6 1.7-1.2.5-.7.8-1.5.8-2.4 0-1.2-.4-2-1.3-2.5s-2.2-.9-4.1-1.2c-1.3-.3-2.4-.6-3.6-1-1.1-.3-2.1-.8-3-1.3-.9-.5-1.5-1.2-2-2.1-.5-.8-.8-1.9-.8-3.2 0-2.4.9-4.2 2.6-5.6 1.7-1.3 4-2 6.8-2.1 1.6 0 3.3.3 5 .8 1.7.6 2.9 1.6 3.7 3.1.4 1.4.6 2.6.6 3.7h-4.6c0-1.8-.6-3-1.7-3.5-1.1-.4-2.1-.6-3.1-.6h-1c-.5 0-1.1.2-1.7.4-.6.2-1.1.5-1.5 1.1-.5.5-.7 1.2-.7 2.1 0 1.1.5 1.9 1.3 2.3.7.4 1.5.7 2.1.9 3.3.7 5.6 1.3 6.9 1.8 1.3.4 2.2 1 2.8 1.7.7.7 1.1 1.4 1.4 2.2.3.8.4 1.6.4 2.5 0 1.4-.3 2.7-.9 3.8-.6 1.1-1.4 2-2.4 2.6-1.1.6-2.2 1-3.4 1.3-1.2.3-2.5.4-3.8.4-2.5 0-4.7-.6-6.6-1.8-1.8-1.2-2.8-3.3-2.9-6.3h5.2zM467.7 50.8h21.9V55h-17.1v11.3h16.3v4.2h-16.3v12.1H490v4.3h-22.3zM499 64.7l-.1-2.9h4.6v4.1c.2-.3.4-.8.7-1.2.3-.5.8-1 1.3-1.5.6-.5 1.4-1 2.3-1.3.9-.3 2-.5 3.4-.5.6 0 1.4.1 2.4.2.9.2 1.9.5 2.9 1.1 1 .6 1.8 1.4 2.5 2.5.6 1.2 1 2.7 1 4.7V87h-4.6V71c0-.9-.1-1.7-.2-2.4-.2-.7-.5-1.3-1.1-1.9-1.2-1.1-2.6-1.7-4.3-1.7-1.7 0-3.1.6-4.3 1.8-1.3 1.2-2 3.1-2 5.7V87H499V64.7zM524.6 61.8h5.1l7.7 19.9 7.6-19.9h5l-10.6 25.1h-4.6zM555.7 50.9h5.5V56h-5.5v-5.1m.5 10.9h4.6v25.1h-4.6V61.8zM570.3 67c0-1.8-.1-3.5-.3-5.1h4.6l.1 4.9c.5-1.8 1.4-3 2.5-3.7 1.1-.7 2.2-1.2 3.3-1.3 1.4-.2 2.4-.2 3.1-.1v4.6c-.2-.1-.5-.2-.9-.2h-1.3c-1.3 0-2.4.2-3.3.5-.9.4-1.5.9-2 1.6-.9 1.4-1.4 3.2-1.3 5.4v13.3h-4.6V67zM587.6 74.7c0-1.6.2-3.2.6-4.8.4-1.6 1.1-3 2-4.4 1-1.3 2.2-2.4 3.8-3.2 1.6-.8 3.6-1.2 5.9-1.2 2.4 0 4.5.4 6.1 1.3 1.5.9 2.7 2 3.6 3.3.9 1.3 1.5 2.8 1.8 4.3.2.8.3 1.5.4 2.2v2.2c0 3.7-1 6.9-3 9.5-2 2.6-5.1 4-9.3 4-4-.1-7-1.4-9-3.9-1.9-2.5-2.9-5.6-2.9-9.3m4.8-.3c0 2.7.6 5 1.8 6.9 1.2 2 3 3 5.6 3.1.9 0 1.8-.2 2.7-.5.8-.3 1.6-.9 2.3-1.7.7-.8 1.3-1.9 1.8-3.2.4-1.3.6-2.9.6-4.7-.1-6.4-2.5-9.6-7.1-9.6-.7 0-1.5.1-2.4.3-.8.3-1.7.8-2.5 1.6-.8.7-1.4 1.7-1.9 3-.6 1.1-.9 2.8-.9 4.8zM620.2 64.7l-.1-2.9h4.6v4.1c.2-.3.4-.8.7-1.2.3-.5.8-1 1.3-1.5.6-.5 1.4-1 2.3-1.3.9-.3 2-.5 3.4-.5.6 0 1.4.1 2.4.2.9.2 1.9.5 2.9 1.1 1 .6 1.8 1.4 2.5 2.5.6 1.2 1 2.7 1 4.7V87h-4.6V71c0-.9-.1-1.7-.2-2.4-.2-.7-.5-1.3-1.1-1.9-1.2-1.1-2.6-1.7-4.3-1.7-1.7 0-3.1.6-4.3 1.8-1.3 1.2-2 3.1-2 5.7V87h-4.6V64.7zM650 65.1l-.1-3.3h4.6v3.6c1.2-1.9 2.6-3.2 4.1-3.7 1.5-.4 2.7-.6 3.8-.6 1.4 0 2.6.2 3.6.5.9.3 1.7.7 2.3 1.1 1.1 1 1.9 2 2.3 3.1.2-.4.5-.8 1-1.3.4-.5.9-1 1.5-1.6.6-.5 1.5-.9 2.5-1.3 1-.3 2.2-.5 3.5-.5.9 0 1.9.1 3 .3 1 .2 2 .7 3 1.3 1 .6 1.7 1.5 2.3 2.7.6 1.2.9 2.7.9 4.6v16.9h-4.6V70.7c0-1.1-.1-2-.2-2.5-.1-.6-.3-1-.6-1.3-.4-.6-1-1.2-1.8-1.6-.8-.4-1.8-.6-3.1-.6-1.5 0-2.7.4-3.6 1-.4.3-.8.5-1.1.9l-.8.8c-.5.8-.8 1.8-1 2.8-.1 1.1-.2 2-.1 2.6v14.1h-4.6V70.2c0-1.6-.5-2.9-1.4-4-.9-1-2.3-1.5-4.2-1.5-1.6 0-2.9.4-3.8 1.1-.9.7-1.5 1.2-1.8 1.7-.5.7-.8 1.5-.9 2.5-.1.9-.2 1.8-.2 2.6v14.3H650V65.1zM700.5 75.2c-.1 2.6.5 4.8 1.7 6.5 1.1 1.7 2.9 2.6 5.3 2.6 1.5 0 2.8-.4 3.9-1.3 1-.8 1.6-2.2 1.8-4h4.6c0 .6-.2 1.4-.4 2.3-.3 1-.8 2-1.7 3-.2.3-.6.6-1 1-.5.4-1 .7-1.7 1.1-.7.4-1.5.6-2.4.8-.9.3-2 .4-3.3.4-7.6-.2-11.3-4.5-11.3-12.9 0-2.5.3-4.8 1-6.8s2-3.7 3.8-5.1c1.2-.8 2.4-1.3 3.7-1.6 1.3-.2 2.2-.3 3-.3 2.7 0 4.8.6 6.3 1.6s2.5 2.3 3.1 3.9c.6 1.5 1 3.1 1.1 4.6.1 1.6.1 2.9 0 4h-17.5m12.8-3v-1.1c0-.4 0-.8-.1-1.2-.1-.9-.4-1.7-.8-2.5s-1-1.5-1.8-2c-.9-.5-2-.8-3.4-.8-.8 0-1.5.1-2.3.3-.8.2-1.5.7-2.2 1.3-.7.6-1.2 1.3-1.6 2.3-.4 1-.7 2.2-.8 3.6h13zM725.7 64.7l-.1-2.9h4.6v4.1c.2-.3.4-.8.7-1.2.3-.5.8-1 1.3-1.5.6-.5 1.4-1 2.3-1.3.9-.3 2-.5 3.4-.5.6 0 1.4.1 2.4.2.9.2 1.9.5 2.9 1.1 1 .6 1.8 1.4 2.5 2.5.6 1.2 1 2.7 1 4.7V87h-4.6V71c0-.9-.1-1.7-.2-2.4-.2-.7-.5-1.3-1.1-1.9-1.2-1.1-2.6-1.7-4.3-1.7-1.7 0-3.1.6-4.3 1.8-1.3 1.2-2 3.1-2 5.7V87h-4.6V64.7zM752.3 61.7h4.3v-5.2l4.5-1.5v6.8h5.4v3.4h-5.4v15.1c0 .3 0 .6.1 1 0 .4.1.7.4 1.1.2.4.5.6 1 .8.4.3 1 .4 1.8.4 1 0 1.7-.1 2.2-.2V87c-.9.2-2.1.3-3.8.3-2.1 0-3.6-.4-4.6-1.2-1-.8-1.5-2.2-1.5-4.2V65.1h-4.3v-3.4zM787.6 86.9c-.3-1.2-.5-2.5-.4-3.7-.5 1-1.1 1.8-1.7 2.4-.7.6-1.4 1.1-2 1.4-1.4.5-2.7.8-3.7.8-2.8 0-4.9-.8-6.4-2.2-1.5-1.4-2.2-3.1-2.2-5.2 0-1 .2-2.3.8-3.7.6-1.4 1.7-2.6 3.5-3.7 1.4-.7 2.9-1.2 4.5-1.5 1.6-.1 2.9-.2 3.9-.2s2.1 0 3.3.1c.1-2.9-.2-4.8-.9-5.6-.5-.6-1.1-1.1-1.9-1.3-.8-.2-1.6-.4-2.3-.4-1.1 0-2 .2-2.6.5-.7.3-1.2.7-1.5 1.2-.3.5-.5.9-.6 1.4-.1.5-.2.9-.2 1.2h-4.6c.1-.7.2-1.4.4-2.3.2-.8.6-1.6 1.3-2.5.5-.6 1-1 1.7-1.3.6-.3 1.3-.6 2-.8 1.5-.4 2.8-.6 4.2-.6 1.8 0 3.6.3 5.2.9 1.6.6 2.8 1.6 3.4 2.9.4.7.6 1.4.7 2 .1.6.1 1.2.1 1.8l-.2 12c0 1 .1 3.1.4 6.3h-4.2m-.5-12.1c-.7-.1-1.6-.1-2.6-.1h-2.1c-1 .1-2 .3-3 .6s-1.9.8-2.6 1.5c-.8.7-1.2 1.7-1.2 3 0 .4.1.8.2 1.3s.4 1 .8 1.5.9.8 1.6 1.1c.7.3 1.5.5 2.5.5 2.3 0 4.1-.9 5.2-2.7.5-.8.8-1.7 1-2.7.1-.9.2-2.2.2-4zM800.7 50.9h4.6V87h-4.6zM828.4 50.8h11.7c2.1 0 3.9.1 5.5.4.8.2 1.5.4 2.2.9.7.4 1.3.9 1.8 1.6 1.7 1.9 2.6 4.2 2.6 7 0 2.7-.9 5.1-2.8 7.1-.8.9-2 1.7-3.6 2.2-1.6.6-3.9.9-6.9.9h-5.7V87h-4.8V50.8m4.8 15.9h5.8c.8 0 1.7-.1 2.6-.2.9-.1 1.8-.3 2.6-.7.8-.4 1.5-1 2-1.9.5-.8.8-2 .8-3.4s-.2-2.5-.7-3.3c-.5-.8-1.1-1.3-1.9-1.7-1.6-.5-3.1-.8-4.5-.7h-6.8v11.9zM858.1 67c0-1.8-.1-3.5-.3-5.1h4.6l.1 4.9c.5-1.8 1.4-3 2.5-3.7 1.1-.7 2.2-1.2 3.3-1.3 1.4-.2 2.4-.2 3.1-.1v4.6c-.2-.1-.5-.2-.9-.2h-1.3c-1.3 0-2.4.2-3.3.5-.9.4-1.5.9-2 1.6-.9 1.4-1.4 3.2-1.3 5.4v13.3H858V67zM875.5 74.7c0-1.6.2-3.2.6-4.8.4-1.6 1.1-3 2-4.4 1-1.3 2.2-2.4 3.8-3.2 1.6-.8 3.6-1.2 5.9-1.2 2.4 0 4.5.4 6.1 1.3 1.5.9 2.7 2 3.6 3.3.9 1.3 1.5 2.8 1.8 4.3.2.8.3 1.5.4 2.2v2.2c0 3.7-1 6.9-3 9.5-2 2.6-5.1 4-9.3 4-4-.1-7-1.4-9-3.9-1.9-2.5-2.9-5.6-2.9-9.3m4.8-.3c0 2.7.6 5 1.8 6.9 1.2 2 3 3 5.6 3.1.9 0 1.8-.2 2.7-.5.8-.3 1.6-.9 2.3-1.7.7-.8 1.3-1.9 1.8-3.2.4-1.3.6-2.9.6-4.7-.1-6.4-2.5-9.6-7.1-9.6-.7 0-1.5.1-2.4.3-.8.3-1.7.8-2.5 1.6-.8.7-1.4 1.7-1.9 3-.7 1.1-.9 2.8-.9 4.8zM904.1 61.7h4.3v-5.2l4.5-1.5v6.8h5.4v3.4h-5.4v15.1c0 .3 0 .6.1 1 0 .4.1.7.4 1.1.2.4.5.6 1 .8.4.3 1 .4 1.8.4 1 0 1.7-.1 2.2-.2V87c-.9.2-2.1.3-3.8.3-2.1 0-3.6-.4-4.6-1.2-1-.8-1.5-2.2-1.5-4.2V65.1h-4.3v-3.4zM927.2 75.2c-.1 2.6.5 4.8 1.7 6.5 1.1 1.7 2.9 2.6 5.3 2.6 1.5 0 2.8-.4 3.9-1.3 1-.8 1.6-2.2 1.8-4h4.6c0 .6-.2 1.4-.4 2.3-.3 1-.8 2-1.7 3-.2.3-.6.6-1 1-.5.4-1 .7-1.7 1.1-.7.4-1.5.6-2.4.8-.9.3-2 .4-3.3.4-7.6-.2-11.3-4.5-11.3-12.9 0-2.5.3-4.8 1-6.8s2-3.7 3.8-5.1c1.2-.8 2.4-1.3 3.7-1.6 1.3-.2 2.2-.3 3-.3 2.7 0 4.8.6 6.3 1.6s2.5 2.3 3.1 3.9c.6 1.5 1 3.1 1.1 4.6.1 1.6.1 2.9 0 4h-17.5m12.9-3v-1.1c0-.4 0-.8-.1-1.2-.1-.9-.4-1.7-.8-2.5s-1-1.5-1.8-2c-.9-.5-2-.8-3.4-.8-.8 0-1.5.1-2.3.3-.8.2-1.5.7-2.2 1.3-.7.6-1.2 1.3-1.6 2.3-.4 1-.7 2.2-.8 3.6h13zM966.1 69.8c0-.3 0-.8-.1-1.4-.1-.6-.3-1.1-.6-1.8-.2-.6-.7-1.2-1.4-1.6-.7-.4-1.6-.6-2.7-.6-1.5 0-2.7.4-3.5 1.2-.9.8-1.5 1.7-1.9 2.8-.4 1.1-.6 2.2-.7 3.2-.1 1.1-.2 1.8-.1 2.4 0 1.3.1 2.5.3 3.7.2 1.2.5 2.3.9 3.3.8 2 2.4 3 4.8 3.1 1.9 0 3.3-.7 4.1-1.9.8-1.1 1.2-2.3 1.2-3.6h4.6c-.2 2.5-1.1 4.6-2.7 6.3-1.7 1.8-4.1 2.7-7.1 2.7-.9 0-2.1-.2-3.6-.6-.7-.2-1.4-.6-2.2-1-.8-.4-1.5-1-2.2-1.7-.7-.9-1.4-2.1-2-3.6-.6-1.5-.9-3.5-.9-6.1 0-2.6.4-4.8 1.1-6.6.7-1.7 1.6-3.1 2.7-4.2 1.1-1 2.3-1.8 3.6-2.2 1.3-.4 2.5-.6 3.7-.6h1.6c.6.1 1.3.2 1.9.4.7.2 1.4.5 2.1 1 .7.4 1.3 1 1.8 1.7.9 1.1 1.4 2.1 1.7 3.1.2 1 .3 1.8.3 2.6h-4.7zM973.6 61.7h4.3v-5.2l4.5-1.5v6.8h5.4v3.4h-5.4v15.1c0 .3 0 .6.1 1 0 .4.1.7.4 1.1.2.4.5.6 1 .8.4.3 1 .4 1.8.4 1 0 1.7-.1 2.2-.2V87c-.9.2-2.1.3-3.8.3-2.1 0-3.6-.4-4.6-1.2-1-.8-1.5-2.2-1.5-4.2V65.1h-4.3v-3.4zM993.5 50.9h5.5V56h-5.5v-5.1m.5 10.9h4.6v25.1H994V61.8zM1006.1 74.7c0-1.6.2-3.2.6-4.8.4-1.6 1.1-3 2-4.4 1-1.3 2.2-2.4 3.8-3.2 1.6-.8 3.6-1.2 5.9-1.2 2.4 0 4.5.4 6.1 1.3 1.5.9 2.7 2 3.6 3.3.9 1.3 1.5 2.8 1.8 4.3.2.8.3 1.5.4 2.2v2.2c0 3.7-1 6.9-3 9.5-2 2.6-5.1 4-9.3 4-4-.1-7-1.4-9-3.9-1.9-2.5-2.9-5.6-2.9-9.3m4.7-.3c0 2.7.6 5 1.8 6.9 1.2 2 3 3 5.6 3.1.9 0 1.8-.2 2.7-.5.8-.3 1.6-.9 2.3-1.7.7-.8 1.3-1.9 1.8-3.2.4-1.3.6-2.9.6-4.7-.1-6.4-2.5-9.6-7.1-9.6-.7 0-1.5.1-2.4.3-.8.3-1.7.8-2.5 1.6-.8.7-1.4 1.7-1.9 3-.6 1.1-.9 2.8-.9 4.8zM1038.6 64.7l-.1-2.9h4.6v4.1c.2-.3.4-.8.7-1.2.3-.5.8-1 1.3-1.5.6-.5 1.4-1 2.3-1.3.9-.3 2-.5 3.4-.5.6 0 1.4.1 2.4.2.9.2 1.9.5 2.9 1.1 1 .6 1.8 1.4 2.5 2.5.6 1.2 1 2.7 1 4.7V87h-4.6V71c0-.9-.1-1.7-.2-2.4-.2-.7-.5-1.3-1.1-1.9-1.2-1.1-2.6-1.7-4.3-1.7-1.7 0-3.1.6-4.3 1.8-1.3 1.2-2 3.1-2 5.7V87h-4.6V64.7zM479.1 100.8h5.2l14.1 36.1h-5.3l-3.8-9.4h-16.2l-3.8 9.4h-5l14.8-36.1m-4.4 22.7H488l-6.5-17.8-6.8 17.8zM508.7 138.8c.1.7.2 1.4.4 1.9.2.6.5 1.1.9 1.6.8.9 2.3 1.4 4.4 1.5 1.6 0 2.8-.3 3.7-.9.9-.6 1.5-1.4 1.9-2.4.4-1.1.6-2.3.7-3.7.1-1.4.1-2.9.1-4.6-.5.9-1.1 1.7-1.8 2.3-.7.6-1.5 1-2.3 1.3-1.7.4-3 .6-3.9.6-1.2 0-2.4-.2-3.8-.6-1.4-.4-2.6-1.2-3.7-2.5-1-1.3-1.7-2.8-2.1-4.4-.4-1.6-.6-3.2-.6-4.8 0-4.3 1.1-7.4 3.2-9.5 2-2.1 4.6-3.1 7.6-3.1 1.3 0 2.3.1 3.2.4.9.3 1.6.6 2.1 1 .6.4 1.1.8 1.5 1.2l.9 1.2v-3.4h4.4l-.1 4.5v15.7c0 2.9-.1 5.2-.2 6.7-.2 1.6-.5 2.8-1 3.7-1.1 1.9-2.6 3.2-4.6 3.7-1.9.6-3.8.8-5.6.8-2.4 0-4.3-.3-5.6-.8-1.4-.5-2.4-1.2-3-2-.6-.8-1-1.7-1.2-2.7-.2-.9-.3-1.8-.4-2.7h4.9m5.3-5.8c1.4 0 2.5-.2 3.3-.7.8-.5 1.5-1.1 2-1.8.5-.6.9-1.4 1.2-2.5.3-1 .4-2.6.4-4.8 0-1.6-.2-2.9-.4-3.9-.3-1-.8-1.8-1.4-2.4-1.3-1.4-3-2.2-5.2-2.2-1.4 0-2.5.3-3.4 1-.9.7-1.6 1.5-2 2.4-.4 1-.7 2-.9 3-.2 1-.2 2-.2 2.8 0 1 .1 1.9.3 2.9.2 1.1.5 2.1 1 3 .5.9 1.2 1.6 2 2.2.8.7 1.9 1 3.3 1zM537.6 125.2c-.1 2.6.5 4.8 1.7 6.5 1.1 1.7 2.9 2.6 5.3 2.6 1.5 0 2.8-.4 3.9-1.3 1-.8 1.6-2.2 1.8-4h4.6c0 .6-.2 1.4-.4 2.3-.3 1-.8 2-1.7 3-.2.3-.6.6-1 1-.5.4-1 .7-1.7 1.1-.7.4-1.5.6-2.4.8-.9.3-2 .4-3.3.4-7.6-.2-11.3-4.5-11.3-12.9 0-2.5.3-4.8 1-6.8s2-3.7 3.8-5.1c1.2-.8 2.4-1.3 3.7-1.6 1.3-.2 2.2-.3 3-.3 2.7 0 4.8.6 6.3 1.6s2.5 2.3 3.1 3.9c.6 1.5 1 3.1 1.1 4.6.1 1.6.1 2.9 0 4h-17.5m12.9-3v-1.1c0-.4 0-.8-.1-1.2-.1-.9-.4-1.7-.8-2.5s-1-1.5-1.8-2.1c-.9-.5-2-.8-3.4-.8-.8 0-1.5.1-2.3.3-.8.2-1.5.7-2.2 1.3-.7.6-1.2 1.3-1.6 2.3-.4 1-.7 2.2-.8 3.7h13zM562.9 114.7l-.1-2.9h4.6v4.1c.2-.3.4-.8.7-1.2.3-.5.8-1 1.3-1.5.6-.5 1.4-1 2.3-1.3.9-.3 2-.5 3.4-.5.6 0 1.4.1 2.4.2.9.2 1.9.5 2.9 1.1 1 .6 1.8 1.4 2.5 2.5.6 1.2 1 2.7 1 4.7V137h-4.6v-16c0-.9-.1-1.7-.2-2.4-.2-.7-.5-1.3-1.1-1.9-1.2-1.1-2.6-1.7-4.3-1.7-1.7 0-3.1.6-4.3 1.8-1.3 1.2-2 3.1-2 5.7V137h-4.6v-22.3zM607 119.8c0-.3 0-.8-.1-1.4-.1-.6-.3-1.1-.6-1.8-.2-.6-.7-1.2-1.4-1.6-.7-.4-1.6-.6-2.7-.6-1.5 0-2.7.4-3.5 1.2-.9.8-1.5 1.7-1.9 2.8-.4 1.1-.6 2.2-.7 3.2-.1 1.1-.2 1.8-.1 2.4 0 1.3.1 2.5.3 3.7.2 1.2.5 2.3.9 3.3.8 2 2.4 3 4.8 3.1 1.9 0 3.3-.7 4.1-1.9.8-1.1 1.2-2.3 1.2-3.6h4.6c-.2 2.5-1.1 4.6-2.7 6.3-1.7 1.8-4.1 2.7-7.1 2.7-.9 0-2.1-.2-3.6-.6-.7-.2-1.4-.6-2.2-1-.8-.4-1.5-1-2.2-1.7-.7-.9-1.4-2.1-2-3.6-.6-1.5-.9-3.5-.9-6.1 0-2.6.4-4.8 1.1-6.6.7-1.7 1.6-3.1 2.7-4.2 1.1-1 2.3-1.8 3.6-2.2 1.3-.4 2.5-.6 3.7-.6h1.6c.6.1 1.3.2 1.9.4.7.2 1.4.5 2.1 1 .7.4 1.3 1 1.8 1.7.9 1.1 1.4 2.1 1.7 3.1.2 1 .3 1.8.3 2.6H607zM629.1 137.1l-3.4 9.3H621l3.8-9.6-10.3-25h5.2l7.6 19.8 7.7-19.8h5z"/>
  #               </svg>
  #             </span>
  #           </a>
  #           <button class="usa-menu-btn usa-button l-header__menu-button">Menu</button>
  #         </div>
  # 
  #     <div class="l-header__nav">
  #       <nav class="usa-nav usa-nav--epa" role="navigation" aria-label="EPA header navigation">
  #         <div class="usa-nav__inner">
  #           <button class="usa-nav__close" aria-label="Close">
  #             <svg class="icon icon--nav-close" aria-hidden="true" role="img">
  #               <title>Primary navigation</title>
  #               <use xlink:href="https://www.epa.gov/themes/epa_theme/images/sprite.artifact.svg#close"></use>
  #             </svg> </button>
  #           <div class="usa-nav__menu">
  #              <ul class="menu menu--main">
  #               <li class="menu__item"><a href="https://github.com/David-Cola/David-Cola.github.io" class="menu__link">Ion Exchange Model</a></li>
  # 
  #             </ul>
  #           </div>
  #         </div>
  #       </nav>
  #     </div>
  #   </header>
  #   <main id="main" class="main" role="main" tabindex="-1">'
  # ),
  
  tags$head(
    tags$style(HTML('.navbar-default .navbar-nav > li > a:hover, .navbar-default .navbar-nav > li > a:focus {
      color: #000; /*Sets the text hover color on navbar*/
    }

    .navbar-default .navbar-nav > .active > a, .navbar-default .navbar-nav > .active >
      a:hover, .navbar-default .navbar-nav > .active > a:focus {
        color: white; /*BACKGROUND color for active*/
          background-color: #0e6cb6;
      }

    .navbar-default {
      background-color: #0e6cb6;
        border-color: #030033;
    }



    .navbar-default .navbar-nav > li > a {
      color: white; /*Change active text color here*/
    }'))),
  
  tags$style(HTML("
    .tabbable > .nav > li > a                  {background-color: #D3D3D3;  color:black}
  ")),
  
  
  
  
  
  
  navbarPage(
    
    div(style ="
          margin-top:-0.4em",
        tags$img(src = "EPA_logo4.png", alt="EPA logo", height="35", width="128"),
        ),
    
             
             #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*#             
             #------------------------------------------------------------------------------#
             #INPUT SECTION#
             #------------------------------------------------------------------------------#            
             #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*#     
             
             tabPanel("Input", 
                      sidebarLayout(
                        sidebarPanel(
                          selectInput("model", "Model Selection", c("HSDM", "PSDM")),
                          fileInput("file1", "Choose .xlsx File", accept = ".xlsx"),
                          textOutput("selectedfile"),
                          textOutput("reject"),
                          textOutput("OutputConcentration"),
                          sliderInput("nrv", "Radial Collocation Points",3, 18, 7),
                          sliderInput("nzv", "Axial Collocation Points", 3, 18, 13),
                          
                          br(),
                          
                          actionButton("run_button", "Run Analysis", icon=icon("play")),
                          textOutput("ionadded"),
                          textOutput("concentrationadded"),
                          textOutput("analysisran")
                        ),
                        
                        
                        
                        mainPanel(
                          tabsetPanel(
                            #------------------------------------------------------------------------------#
                            #COLUMN PARAMETERS#
                            #------------------------------------------------------------------------------#                                 
                            tabPanel("Column Parameters",
                                     
                                     br(),
                                     
                                     #--------------------------Resin Characteristics-------------------------------#  
                    
                                     
                                    

                              
                                     br(),
                                     
                                     
                                  
                                    
                                     fluidRow(
                                       column(3,HTML(paste0("<h4>","<strong>", "Resin Characteristics", "</strong>", "</h4>"))),
                                       column(1,),
                                       column(2,shinyWidgets::autonumericInput(
                                         inputId = "Qv",
                                         label="Resin Capacity",
                                         value = 1400,
                                         decimalPlaces = 2,
                                         digitGroupSeparator = ",",
                                         decimalCharacter = ".")),
                                       column(3, selectInput("qunits", "Resin Capacity Units", c("meq/L")))),
                                         

                                     fluidRow(
                                       column(3, ),
                                       column(1,),
                                       column(2, shinyWidgets::autonumericInput(
                                         inputId = "rbv",
                                         label="Bead Radius",
                                         value = 0.03375,
                                         decimalPlaces = 5,
                                         digitGroupSeparator = ",",
                                         decimalCharacter = ".")),
                                       column(3, selectInput("rbunits", "Bead Radius Units", c("cm", "m", "mm", "in", "ft")))),
                                        

                                     fluidRow(
                                       column(3, ),
                                       column(1,),
                                       column(2, shinyWidgets::autonumericInput(
                                         inputId = "EBEDv",
                                         label="Bed Porosity",
                                         value = 0.35,
                                         currencySymbolPlacement = "p",
                                         decimalPlaces = 3,
                                         digitGroupSeparator = ",",
                                         decimalCharacter = "."
                                       ))),
                                      

                                     fluidRow(
                                       column(3, ),
                                       column(1,),
                                       column(2, shinyWidgets::autonumericInput(
                                         inputId = "EPORv",
                                         label="Bead Porosity",
                                         value = 0.2,
                                         currencySymbolPlacement = "p",
                                         decimalPlaces = 3,
                                         digitGroupSeparator = ",",
                                         decimalCharacter = "."
                                       ))),
                                     
                                     #------------------------------------------------------------------------------#                                       
                                     
                                     hr(),
                                     #Parameters Row 2#
                                     #--------------------------Column Specifications-------------------------------#                                     
                                     fluidRow(
                                       column(3,
                                           
                                              HTML(paste0("<h4>","<strong>", "Column Specifications", "</strong>", "</h4>")),
                                        
                                              #This radio button toggles between Linear and volumetric flowrate
                                              br(),
                                              radioButtons("veloselect", "", c("Linear", "Volumetric"))),
                                       
                                       column(1,),
                                       
                                       column(2, #offset=1,
                                                         
                                              shinyWidgets::autonumericInput(
                                                inputId = "Lv",
                                                label="Length",
                                                value = 14.765, 
                                                decimalPlaces = 3,
                                                digitGroupSeparator = ",",
                                                decimalCharacter = "."),
                                              

                                              shinyWidgets::autonumericInput(
                                                inputId = "Vv",
                                                label="Velocity",
                                                value = 0.123, 
                                                decimalPlaces = 3,
                                                digitGroupSeparator = ",",
                                                decimalCharacter = "."),
                                              
                                              shinyWidgets::autonumericInput(
                                                inputId = "Dv",
                                                label="Diameter",
                                                value = 4, 
                                                decimalPlaces = 3,
                                                digitGroupSeparator = ",",
                                                decimalCharacter = "."),
                                              
                                         
                                              shinyWidgets::autonumericInput(
                                                inputId = "Fv",
                                                label="Flow Rate",
                                                value = 1.546, 
                                                decimalPlaces = 5,
                                                digitGroupSeparator = ",",
                                                decimalCharacter = ".")),
                                       column(3,
                                              
                                              selectInput("LengthUnits", "Length Units", c("cm", "m", "mm", "in", "ft")),
                                              div(style ="
                                               margin-top:-0.33em", 
                                              selectInput("VelocityUnits", "Velocity Units", c("cm/s", "m/s", "m/min", "m/h", "in/s","ft/s","ft/min", "gpm/ft^2"))),
                                              div(style ="
                                               margin-top:-0.33em",          
                                       selectInput("DiameterUnits","Diameter Units",c("cm", "m", "in", "ft"))),
                                       div(style ="
                                               margin-top:-0.33em", 
                                       selectInput("FlowrateUnits","Flow Rate Units",c("cm^3/s", "m^3/s", "ft^3/s", "mL/s", "L/min", "mL/min", "gpm", "mgd"))))),
                                              
                                            
                                       # column(3,
                                       #        selectInput("LengthUnits", "", c("cm", "m", "mm", "in", "ft")),
                                       #        div(style="margin-top:-0.5em",
                                       #            selectInput("VelocityUnits", "", c("cm/s", "m/s", "m/min", "m/h", "in/s","ft/s","ft/min", "gpm/ft^2"))),
                                       #        
                                       #        div(style ="
                                       #        margin-top:-0.5em", 
                                       #            selectInput("DiameterUnits","",c("cm", "m", "in", "ft")),
                                       #            selectInput("FlowrateUnits","",c("cm^3/s", "m^3/s", "ft^3/s", "mL/s", "L/min", "mL/min", "gpm", "mgd"))))),
                                     #------------------------------------------------------------------------------#                                     
                                     
                                     hr(),
                                     #Parameters Row 4#
                                     
                                     fluidRow(
                                       column(3,
                                              HTML(paste0("<h4>","<strong>", "Concentration Time", "</strong>", "</h4>"))),
                                       column(3,),
                                       column(3, 
                                                  selectInput("timeunits2", "Time Units", c("hr", "day")))),
                                           
                                     
                                     #------------------------------------------------------------------------------#
                                     #END COLUMN PARAMETERS#
                                     #------------------------------------------------------------------------------#                                       
                                     
                                     #------------------------------------------------------------------------------#
                                     #IONS TAB#
                                     #------------------------------------------------------------------------------#                                     
                                     
                            ),
                            tabPanel("Ions",
                                     br(),
                                     h4("Ion List"),
                                     dataEditUI("edit-1"),
                                     br(), br(),
                                     h4("Concentration Points"),
                                     dataEditUI("edit-2"),
                                     br(), br(),
                                     h4("Effluent Data"),
                                     dataEditUI("edit-3")
                            ),
                            
                            #------------------------------------------------------------------------------#
                            #END IONS TAB#
                            #------------------------------------------------------------------------------#
                            
                            #------------------------------------------------------------------------------#
                            #ALKALINITY TAB#
                            #------------------------------------------------------------------------------#
                            
                            tabPanel("Alkalinity",
                                     br(),
                                     h4("Bicarbonate Concentration of Alkalinity"),
                                     textOutput("AlkConv"),
                                     br(),
                                     fluidRow(
                                       column(4,
                                              numericInput("alkvalue", "Alkalinity Value", 5),
                                              numericInput("pH", "pH", 7)),
                                       column(4, offset=1,
                                              selectInput("alkunits", "Concentration Units", c("meq", "mg/L")),
                                              
                                              div(style ="
                                              margin-top:2em",
                                                  h5("Bicarbonate Concentration (meq)")),
                                              div(style ="
                                              margin-top:-1em",
                                                  textOutput("bicarbcin"))),
                                       column(4, offset=1,
                                              div(style ="
                                              margin-top:2em",
                                                  h5("Bicarbonate Concentration (mg/L)")),
                                              div(style ="
                                              margin-top:-1em",
                                                  textOutput("bicarbcinmgl"))),
                                       br()
                                       
                                     )#fluid row
                            ),#tab panel
                            #------------------------------------------------------------------------------#
                            #END ALKALINITY TAB#
                            #------------------------------------------------------------------------------#
                          )
                        )#mainPanel
                      ),#sidebarLayout
                      
             ),#navbarPage
             
             #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*#             
             #------------------------------------------------------------------------------#
             #END INPUT SECTION#
             #------------------------------------------------------------------------------#            
             #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*# 
             
             #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*#             
             #------------------------------------------------------------------------------#
             #OUTPUT SECTION#
             #------------------------------------------------------------------------------#            
             #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*#   
             
             tabPanel("Output",
                      
                      sidebarLayout(
                        sidebarPanel(
                          selectInput("OCunits", "Output Concentration Units", c("mg/L", "ug/L", "ng/L", "c/c0")),
                          selectInput("timeunits","Output Time Units",c("Days", "Bed Volumes (x1000)", "Hours", "Months", "Years")),
                          
                          checkboxInput("computeddata", "Computed Data", TRUE),
                          checkboxInput("effluentdata", "Effluent Data", FALSE),
                          checkboxInput("influentdata", "Influent Data", FALSE),
                          
                          downloadButton("save_button", "Save Data")
                        ),
                        
                        mainPanel(
                          
                          shinycssloaders::withSpinner(
                            plotlyOutput("Plot")),#Counterions
                          br(),
                          textOutput("CounterIonPlot"),
                          br(),
                          plotlyOutput("ExtraChemicals"),
                          br(),
                          textOutput("IonPlot")))), #Ions
             
             #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*#             
             #------------------------------------------------------------------------------#
             #END OUTPUT SECTION#
             #------------------------------------------------------------------------------#            
             #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*#
             
             #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*#             
             #------------------------------------------------------------------------------#
             #ABOUT SECTION#
             #------------------------------------------------------------------------------#            
             #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*#
             
             tabPanel("About",
                      h5("Ion Exchange Model"),
                      textOutput("about"),
                      br(),
                      tags$a(href="https://github.com/USEPA/Water_Treatment_Models/", "Read more about the Ion Exchange Model"),
                      br(), br(),
                      #textOutput("how2use"),
                      h5("There are two ways to start this model:"),
                      textOutput("how2use2"),
                      br(),
                      textOutput("how2use3"),
                      textOutput("how2use4"),
                      br(),
                      #textOutput("how2use5"),
                      h5("Developed By"),
                      textOutput("how2use6"),
                      textOutput("how2use7"),
                      textOutput("how2use8"))
  )
)


#~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*#
#------------------------------------------------------------------------------#
#~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*#
#------------------------------------------------------------------------------#
#BEGIN SERVER#
#------------------------------------------------------------------------------#
#~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*#
#------------------------------------------------------------------------------#
#~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*#

server <- function(input, output, session) {
  
  #------------------------------------------------------------------------------#
  #STATIC DISPLAY TEXTS#
  #------------------------------------------------------------------------------#
  
  #--------------------------Resin Characteristics-------------------------------#  
  output$RC<-renderText("Resin Characteristics")
  output$Q<-renderText("Resin Capacity")
  output$rb<-renderText("Bead Radius")
  output$EBED<-renderText("Bed Porosity")
  output$name<-renderText("Name")
  #------------------------------------------------------------------------------#
  
  #The reason these are multiple if else statements is because it gives more
  #control of where the items will appear in the UI
  
  output$EPOR<-renderUI({
    renderText("Bead Porosity")
  })
  
  # output$EPORv<-renderUI({
  #     numericInput("EPORvalue", "", 0.2)
  # })
  
  observe({
    toggleState("EPORv", condition=input$model!="HSDM")
  })
  
  
  
  
  #------------------------------------------------------------------------------#
  
  output$CS<-renderText("Column Specifications")
  output$MC<-renderText("Material Characteristics")
  output$CS3<-renderText("Solver Related")
  #--------------------------Column Specifications-------------------------------#  
  output$Length<-renderText("Length")
  output$Velocity<-renderText("Velocity")
  output$Diameter<-renderText("Diameter")
  output$Flowrate<-renderText("Flow Rate")
  #------------------------------------------------------------------------------#  
  output$SR<-renderText("Solver Related")
  output$nr<-renderText("Radial Collocation Points")
  output$nz<-renderText("Axial Collocation Points")
  
  output$conctime<-renderText("Time")
  output$TS<-renderText("Time Step")
  
  output$ChemicalNames<-renderText("Chemical Names")
  output$mw<-renderText("mw")
  output$KxA<-renderText("KxA")
  output$Valence<-renderText("Valence")
  
  output$Name2<-renderText("Name")
  output$InitialTime<-renderText("Inital")
  output$FinalTime<-renderText("Final")
  
  output$OC<-renderText("Units")
  
  output$IonList<-renderText("Ion List")
  output$ConcentrationList<-renderText("Concentration Points")
  
  output$about<-renderText("The Ion Exchange Model is a tool used to model a strong-base anion exchange unit operation in a drinking water treatment plant. This model relies on selectivity coefficient parameters and other information about the anion exchange resin and predicts the breakthrough behavior for unit operation design.")
  
  output$AlkConv<-renderText("Bicarbonate is the common chemical used to measure alkalinity in this model, however, the user may have the pH of their water without the Bicarbonate specifications. If this is the case then the user can use this calculator to take their pH measurement and find the corresponding Bicarbonate concentrations.  ")
  output$bicarbion<-renderTable(bicarbion)
  
  #output$selectedfile<-renderText(input$file1)
  
  
  output$how2use2<-renderText("1) Use an Excel file to describe parameters of water treatment unit operation (examples provided). One can upload such file by clicking 'Browse' in the top left corner of the Input page.")
  output$how2use3<-renderText("2) Start with the data that is provided in the user interface and manipulate the data from there. Once the parameters have been decided ions can be added, either in the xlsx file or on the ions tab, as well as concentration points. When the user is satisfied with their settings, click 'run analysis' to begin the computation. Simulation time can take a few seconds to minutes depending on how many ions are added.")
  output$how2use4<-renderText(" Once the parameters have been decided ions can be added, either in the xlsx file or on the ions tab, as well as concentration points. When the user is satisfied with their settings, click 'run analysis' to begin the computation. Simulation time can take a few seconds to minutes depending on how many ions are added.")
  output$how2use5<-renderText("Developed By")
  output$how2use6<-renderText("David Colantonio")
  output$how2use7<-renderText("Levi Haupert")
  output$how2use8<-renderText("Jonathan Burkhardt")
  
  output$CounterIonPlot<-eventReactive(input$run_button, {
    print("Above is a graph potraying the concentration of counter ion chemicals over time using the Ion Exchange technique.")
  })
  
  output$IonPlot<-eventReactive(input$run_button, {
    print("Above is a graph protraying the concentration of ion chemicals over time using the Ion Exchange tenchnique.")
  })
 
  #------------------------------------------------------------------------------#
  #INPUT FILE HANDLING#
  #------------------------------------------------------------------------------#
  
  #The function process_files was defined at the beginning, now it is being called
  #When a file is inputted into the UI
  observeEvent(input$file1, {
    file <- input$file1
    process_files(file$datapath)
  })
  
  #GUI rejects a file upload that is not an xlsx
  output$reject<-renderPrint({
    req(input$file1)
    if(input$file1$type != "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"){ stop("Please upload a .xlsx file")}
  })
  
  ### TODO: Should the reject be within the process_files? Or somehow before that.... so merging the above 2 items?
  
  #When a file is uploaded, the session is reloaded. We do this because there does
  #Not seem to be any other way to overwrite the DataEditR tables. See the process_file
  #Notes for further elaboration.
  observeEvent(input$file1, {
    session$reload()
  })
  
  
  
  #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*#             
  #------------------------------------------------------------------------------#
  #DATA PREP SECTION#
  #------------------------------------------------------------------------------#            
  #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*#  
  
  
  
  
  
  #------------------------------------------------------------------------------#
  #PARAMS DATA HANDLING#
  #------------------------------------------------------------------------------#  
  
  #When the param sheet is read in, make it a reactiveVal so that the data can
  #Be saved between refreshes and edited
  paramsheet<-reactiveVal(read.csv("paramsheet.csv"))
  
  
  ##Flow rate V. Linear Velocity
  ##Some water treatment users may want to use a linear velocity and some may want to use a flow rate
  ##Given that we can have one option, both options, or neither
  
  
  test_df<-data.frame(C=c('v','flrt','diam'))
  flags<-reactive({test_df$C %in% paramsheet()$name}) ##flags are in order [1] velocity [2] flowrate and [3] diameter
  
  velocity<-reactiveVal()
  velocityvector2<-reactiveVal()
  velocityvector3<-reactiveVal()
  
  flowrate<-reactiveVal()
  flowrate2<-reactiveVal()
  flowrate3<-reactiveVal()
  
  diameter<-reactiveVal()
  diameter2<-reactiveVal()
  diameter3<-reactiveVal()
  
  
  observe({if (flags()[1]){
    # velocity read in
    velocity(filter(paramsheet(), name=='v')$value)
    updateNumericInput(session, "Vv", value=velocity())
    
    velocityvector2(c(filter(paramsheet(), name=='v')$units, velocityvector))
    velocityvector3<-unique(velocityvector2())
    
    updateSelectInput(session, "VelocityUnits", choices=velocityvector3())
    
    ##add toggle of velocity selector
    updateRadioButtons(session, "veloselect", selected="Linear")
    
  }
    else if(flags()[2] & flags()[3]){
      
      flowrate(filter(paramsheet(), name=='flrt')$value)
      diameter(filter(paramsheet(), name=='diam')$value)
      
      updateNumericInput(session, "Fv", value=flowrate())
      updateNumericInput(session, "Dv", value=diameter())
      
      
      flowrate2(c(filter(paramsheet(), name=='flrt')$units, flowratevector))
      flowrate3(unique(flowrate2()))
      
      diameter2(c(filter(paramsheet(), name=='diam')$units, diametervector))
      diameter3(unique(diameter2()))
      
      updateSelectInput(session, "FlowrateUnits", choices=flowrate3())
      updateSelectInput(session, "DiameterUnits", choices=diameter3())                                     
      
      updateRadioButtons(session, "veloselect", selected="Volumetric")
    }
    else{
      print("Warning: No flow data provided, defaults used")
    }
    
  })
  
  
  
  #Take the data from the file that the user uploaded and overwrite the default frame
  capacity<-reactive({filter(paramsheet(), name=="Q")$value})
  eebed<-reactive({filter(paramsheet(), name=="EBED")$value})
  length2<-reactive({filter(paramsheet(), name=="L")$value})
  beadradius<-reactive({filter(paramsheet(), name=="rb")$value})
  film<-reactive({filter(paramsheet(), name=="kL")$value})
  diffuse<-reactive({filter(paramsheet(), name=="Ds")$value})
  radial<-reactive({filter(paramsheet(), name=="nr")$value})
  axial<-reactive({filter(paramsheet(), name=="nz")$value})
  time<-reactive({filter(paramsheet(), name=="time")$value})
  
  observe({
    updateNumericInput(session, "Qv", value=format(capacity(), digits=4, scientific=FALSE))
    updateNumericInput(session, "EBEDv", value=format(eebed(), digit=4, scientific=FALSE))
    updateNumericInput(session, "Lv", value=length2())
    updateNumericInput(session, "rbv", value=prettyNum(beadradius(), digits=4, scientific=FALSE))
    updateNumericInput(session, "kLv", value=film())
    updateNumericInput(session, "Dsv", value=diffuse())
    updateNumericInput(session, "nrv", value=radial())
    updateNumericInput(session, "nzv", value=axial())
  })
  
  
  
  #This is the to update the selectInput to the unit that you have in the excel
  #file. There doesn't seem to be a direct way to update to the selectInput to the one in
  #the xlsx file, but, you can update the options as a whole. So, I create a vector
  #of all the unit options, then take the unit you have in your file and combine them
  #into one vector. The value in your excel file becomes the first value in the combined vector. 
  #Then from there
  #I run it through a unique() function which gets rid of the duplicate in your vector
  #Then, I can update the selectInput vector where it has reorganized your  initial
  #vector to start
  #with the users input.
  
  
  modelvec<-reactive({
    
    model<-c(filter(paramsheet(), name=="model")$value, modelvector)
    updatedmodel<-unique(model)
    
    updatedmodel
    
  })
  
  
  lengthvector2<-reactive({c(filter(paramsheet(), name=="L")$units, lengthvector)})
  lengthvector3<-reactive({unique(lengthvector2())})
  
  rbvector<-reactive({c(filter(paramsheet(), name=='rb')$units, lengthvector)})
  rbvector2<-reactive({unique(rbvector())})
  
  timevector2<-reactive(c(filter(paramsheet(), name=='time')$units, timevector))
  timevector3<-reactive({unique(timevector2())})
  
  observe({
    updateSelectInput(session, "rbunits", choices=rbvector2())
    updateSelectInput(session, "LengthUnits", choices=lengthvector3())
    updateSelectInput(session, "timeunits2", choices=timevector3())
    updateSelectInput(session, "model", choices=modelvec())
  })
  
  observe({
    toggleState("Vv", condition=input$veloselect!="Volumetric")
    toggleState("Fv", condition=input$veloselect!="Linear")
    toggleState("Dv", condition=input$veloselect!="Linear")
  })
  
  velocityvar<-reactiveVal()
  
  #------------------------------------------------------------------------------#
  #IBICARBONATE TO ALKALINITY CONVERTER#
  #------------------------------------------------------------------------------#  
  
  
  bicarbconverted<-reactiveVal()
  bicarbmeq2mgl<-50.045001
  
  h_plus<-reactiveVal()
  observe({h_plus(10^-input$pH)})
  
  calcium_carb_alpha<-reactive({k1*h_plus()/(h_plus()**2 + k1*h_plus()+k1*k2)})
  
  observe({
    if(input$alkunits=='meq'){
      bicarbconverted(calcium_carb_alpha()*input$alkvalue)
    }
    else{
      bicarbconverted(calcium_carb_alpha()*input$alkvalue/bicarbmeq2mgl) #mw/valence -> mw
    }
  })
  
  
  output$bicarbcin<-renderText(bicarbconverted())
  output$bicarbcinmgl<-renderText(bicarbconverted()*bicarbmeq2mgl)
  
  
  
  #------------------------------------------------------------------------------#
  #IONS TAB DATA HANDLING#
  #------------------------------------------------------------------------------#  
  
  iondat<- dataEditServer("edit-1", data = "ionsheet.csv")
  dataOutputServer("output-1", data = iondat)
  
  
  #------------------------------------------------------------------------------#
  #CIN TAB DATA HANDLING#
  #The prep here seems like it can mostly be done mostly in one or two functions
  #instead of the multiple functions that are used here, but, some of the 
  #independent steps are used throughout the code. So this makes it easier
  #To break it up
  #------------------------------------------------------------------------------#   
  
  cindat<-dataEditServer("edit-2",read_args=list(colClasses=c("numeric")),data="cinsheet.csv") ## read_args should make all columns numeric, which seems to address the "initial read in as integer issues"
  dataOutputServer("output-2", data = cindat)
  
  #Convert the cin data time to hours if it is not already
  cindat_hours<-reactive({
    cindata<-cindat()
    time<-cindata['time']
    if(input$timeunits2=='hr'){
      newtime<-time
    }
    else{
      newtime<-time*24
    }
    newcin<-cbind(newtime, cindata[,2:ncol(cindata)])
    return(newcin)
  })
  
  
  #convert cindat to meq if it is not already
  cin_hours_meq<-eventReactive(input$run_button,{cin_correct(iondat(), cindat_hours())})
  cin_hours_mgl<-eventReactive(input$run_button,{mass_converter_mgl(iondat(), cindat_hours())})
  
  #When the file is first uploaded the influent data and simulated data both have
  #the same names, so to differentiate them I rename them to chemical_influent 
  #with this function
  cin_hours_meq_renamer<-eventReactive(input$run_button, {influent_chemical_renamer(cin_hours_meq(), cindat_hours())})
  cin_hours_mgl_renamer<-eventReactive(input$run_button, {influent_chemical_renamer(cin_hours_mgl(), cindat_hours())})
  
  #The cin tab is now in the correct units and named appropriately. Gather then 
  #brings them to a shape that makes it easy to convert and easy to plot
  #Using the tidyr::gather function with the hours still attached gives a bad 
  #result, so the time is Seperated out and then reattached in 
  #cindat_meq_hours_preprepped
  
  cin_meq_hours_prep<-eventReactive(input$run_button,{influent_organizer(cin_hours_meq_renamer(), cindat_hours())})
  cin_mgl_hours_prep<-eventReactive(input$run_button,{influent_organizer(cin_hours_mgl_renamer(), cindat_hours())})
  
  
  
  
  #------------------------------------------------------------------------------#
  #EFFLUENT TAB HANDLING#
  #------------------------------------------------------------------------------# 
  
  effluentdat<-dataEditServer("edit-3",read_args=list(colClasses=c("numeric")), data="effluent.csv")
  dataOutputServer("output-1", data=effluentdat)
  
  #Put effluent data into plot data format
  effdata<-reactive({effluent_data_processor(iondat(), effluentdat())})
  
  #Make effluent time match input time
  effdat_hours<-reactive({
    effluentdata<-effdata()
    time<-effluentdata['time']
    if(input$timeunits2=='hr'){
      newtime<-time
    }
    else{
      newtime<-time*24
    }
    neweff<-cbind(newtime, effluentdata[,2:ncol(effluentdata)])
    return(neweff)
  })
  
  
  
  
  
  
  #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*#             
  #------------------------------------------------------------------------------#
  #END DATA PREP SECTION#
  #------------------------------------------------------------------------------#            
  #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*#  
  
  
  
  
  
  
  #------------------------------------------------------------------------------#
  #CALLING THE HSDMIX FUNCTION#
  #------------------------------------------------------------------------------#
  
  #HSDMIX values is stored in this reactiveVal "out"
  out<-reactiveVal()
  
  observeEvent(input$run_button, {
    if(input$model=="HSDM"){
      out(HSDMIX_prep(input, iondat(), cindat(), nt_report))
    }
    else{
      out(PSDMIX_prep(input, iondat(), cindat(), nt_report))
    }
  })
  
  
  
  # find outlet indices
  outlet_id <- reactive({dim(out()[[2]])[4]})
  liquid_id <- reactive({dim(out()[[2]])[2]})
  
  
  #------------------------------------------------------------------------------#
  #                 IEX CONCENTRATION OUTPUT DATAFRAME
  #------------------------------------------------------------------------------#
  ### TODO: add better error handling HSDMIX_prep can now return an 'error' value which is an integer, or the full data
  #### only want to proceed if it isn't an error state
  
  timeframe<-reactive({data.frame(hours=out()[[1]])})
  allchemicalconcs<-list()
  
  
  #HSDMIX outputs a list, so this takes the list and binds them into a dataframe
  allchemicals_hours_meq<-eventReactive(input$run_button, {for (x in 1:nrow(iondat())){
    conc<-out()[[2]][, liquid_id(), x, outlet_id()]
    allchemicalconcs[[x]]<-conc
  }
    allconcdf<-data.frame(allchemicalconcs)
    colnames(allconcdf)<-iondat()$name
    allconcdf
  })
  
  allchemicals_hours_mgl<-reactive({HSDMIX_in_hours_mgl(allchemicals_hours_meq(), iondat(), timeframe())})
  
  
  #------------------------------------------------------------------------------#
  #                   END IEX CONCENTRATION OUTPUTDATAFRAME#
  #------------------------------------------------------------------------------#
  #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*#
  
  #------------------------------------------------------------------------------#
  #GENERATE C/C0 DATAFRAMES FROM MG/L CONCENTRATION DATAFRAMES#
  
  #The goal here is: For each column in allchemicals, divide every element
  #in each column by the first value in that column.
  #------------------------------------------------------------------------------#
  ## need to create meq/L c0 for math to work out
  cc0vector_meq<-reactive({cc0_conv_meq(iondat(), cindat())})
  
  
  computedcc0<-reactive({HSDMIX_cc0(allchemicals_hours_meq(), cc0vector_meq())})
  effluentcc0<-reactive({effluent_cc0(effluentdat(), iondat(), cc0vector_meq())})
  influentcc0<-reactive({HSDMIX_cc0(cin_hours_meq_renamer()[,2:ncol(cin_hours_meq_renamer())], cc0vector_meq())})
  
  
  
  #------------------------------------------------------------------------------#
  #END IEX CONCENTRATION OUTPUTDATAFRAME#
  #------------------------------------------------------------------------------#
  #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*#
  
  
  
  
  
  
  #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*#
  #------------------------------------------------------------------------------#
  #CONVERTING OUTPUT DATAFRAMES#
  #------------------------------------------------------------------------------#
  
  
  outputcounterions<-reactiveValues()
  outputions<-reactiveValues()
  outputeffluent<-reactiveValues()
  outputinfluent<-reactiveValues()
  
  
  ion_list<-c("CHLORIDE", "SULFATE", "BICARBONATE", "NITRATE")
  ion_flag<-reactive({ion_list %in% colnames(cindat())})
  number_of_ions<-reactive({length(ion_flag()[ion_flag()==TRUE])})
  
  counterIon_loc<-reactive({ number_of_ions() * nt_report })
  addIon_loc<-reactive({ counterIon_loc() + 1 })
  
  
  counteriondata<-reactive({allchemicals_hours_mgl()[0:counterIon_loc(),]})
  iondata<-reactive({allchemicals_hours_mgl()[addIon_loc():nrow(allchemicals_hours_mgl()),]})
  
  counteriondatacc0<-reactive({computedcc0()[0:counterIon_loc(),]})
  iondatacc0<-reactive({computedcc0()[addIon_loc():nrow(allchemicals_hours_mgl()),]})
  
  outputcounterions$name<-reactive({counteriondata()$name})
  outputions$name<-reactive({iondata()$name})
  
  
  
  
  
  observe({
    ## convert time units for graphing
    # calculating kBV
    if (input$timeunits == "Bed Volumes (x1000)") {
      bv_conv <- get_bv_in_sec(input)
      outputcounterions$time <- counteriondata()$hours / (bv_conv / hour2sec) / 1e3
      outputions$time <- iondata()$hours / (bv_conv / hour2sec) / 1e3
      
      outputeffluent$time<- effdat_hours()$time/ (bv_conv / hour2sec) / 1e3
      outputinfluent$hours<-cin_mgl_hours_prep()$hours  / (bv_conv / hour2sec) / 1e3  ## should this be $time?
      
    } else {
      outputcounterions$time <- counteriondata()$hours / (time_conv[input$timeunits] / hour2sec)
      outputions$time <- iondata()$hours / (time_conv[input$timeunits] / hour2sec)
      
      outputeffluent$time<- effdat_hours()$time/ (time_conv[input$timeunits] / hour2sec)
      outputinfluent$hours<-cin_mgl_hours_prep()$hours/ (time_conv[input$timeunits] / hour2sec) ## should this be $time?
    }
  })
  
  
  
  
  
  
  observe({
    ### convert y-axis/mass units for graphing
    if(input$OCunits=="c/c0"){
      ## just replicates the returned data
      outputcounterions$conc <-  counteriondatacc0()$conc
      outputions$conc<-iondatacc0()$conc
      
      outputeffluent$conc<- effluentcc0()$conc
      outputinfluent$conc <- influentcc0()$conc#influentcc04()$conc
    } else {
      outputcounterions$conc <- counteriondata()$conc / mass_conv[input$OCunits]
      outputions$conc <- iondata()$conc / mass_conv[input$OCunits]
      
      outputeffluent$conc <- effdat_hours()$conc/mass_conv[input$OCunits]
      outputinfluent$conc <- cin_mgl_hours_prep()$conc/mass_conv[input$OCunits]
    }
    
  })
  
  
  
  
  
  
  
  
  ### graph data
  
  counterion_data_processed<-reactive({
    if(input$computeddata==TRUE){
      plot_data <- counteriondata()
      plot_data$conc <- outputcounterions$conc
      plot_data$hours <- outputcounterions$time
      plot_data
    }
    else{
      plot_data <- data.frame(hours=c(NA), name=c(NA), conc=c(NA))
      plot_data
    }
  })
  
  ion_data_processed<-reactive({
    if(input$computeddata==TRUE){
      plot_data2 <- iondata()
      plot_data2$conc <- outputions$conc
      plot_data2$hours <- outputions$time
      plot_data2
    }
    else{
      plot_data2 <- data.frame(hours=c(NA), name=c(NA), conc=c(NA))
      plot_data2
    }
  })
  
  
  
  effluent_processed<-reactive({
    if(input$effluentdata==TRUE){
      plot_data3<-effdat_hours()
      plot_data3$conc<-outputeffluent$conc
      plot_data3$hours<- outputeffluent$time
      plot_data3
    }
    else{
      plot_data3 <- data.frame(hours=c(NA), name=c(NA), conc=c(NA))
      plot_data3
    }
  })
  
  
  
  
  influent_processed<-reactive({
    if(input$influentdata==TRUE){
      plot_data4<-cin_mgl_hours_prep()
      plot_data4$conc<-outputinfluent$conc
      plot_data4$hours<-outputinfluent$hours
      plot_data4
    }
    else{
      plot_data4 <- data.frame(hours=c(NA), name=c(NA), conc=c(NA))
      plot_data4
    }
  })
  
  
  
  cindat_converter_counter<-reactive({
    dat<-influent_processed()
    newdat<-subset(dat, name=="SULFATE_influent" | name=="CHLORIDE_influent"
                   | name=="NITRATE_influent" | name=="BICARBONATE_influent")
    newdat
  })
  
  cindat_converter_ion<-reactive({
    dat<-influent_processed()
    newdat<-subset(dat, name!="SULFATE_influent" & name!="CHLORIDE_influent"
                   & name!="NITRATE_influent" & name!="BICARBONATE_influent")
    newdat
  })
  
  
  
  fig<-reactive({create_plotly(counterion_data_processed(), effluent_processed(), cindat_converter_counter())})
  counterionfigure<-reactive({fig()%>%layout(title="Concentration over Time", showlegend=TRUE,
                                             legend=list(orientation='h', y=1),
                                             xaxis=list(title=input$timeunits, gridcolor = 'ffff'),
                                             yaxis=list(title=paste0("Concentration (",input$OCunits,")"), showexponent='all',
                                                        exponentformat='e', gridcolor = 'ffff'))})
  
  bonusfig<-reactive({create_plotly2(ion_data_processed(), effluent_processed(), cindat_converter_ion())})
  ionfigure<-reactive({bonusfig()%>%layout(title="Concentration over Time", showlegend=TRUE,
                                           legend=list(orientation='h', y=1),
                                           xaxis=list(title=input$timeunits, gridcolor = 'ffff'),
                                           yaxis=list(title=paste0("Concentration (",input$OCunits,")"), showexponent='all',
                                                      exponentformat='e', gridcolor = 'ffff'))})
  
  
  
  output$Plot<-renderPlotly(
    counterionfigure())
  
  output$ExtraChemicals <- renderPlotly(
    ionfigure())
  
  
  
  
  outputOptions(output, "Plot", suspendWhenHidden = FALSE)
  outputOptions(output, "ExtraChemicals", suspendWhenHidden = FALSE)
  
  paramdf<-reactive({data.frame(name=c("Q", "EBED", "L", "v", "rb", "kL", "Ds", "nr", "nz", "time"),
                                value=c(input$Qv, input$EBEDv, input$Lv, input$Vv, input$rbv, NA, NA, input$nrv, input$nzv, input$timeunits2),
                                units=c(input$qunits, NA, input$LengthUnits, input$VelocityUnits, input$rbunits, NA, "cm2/s", NA, NA, input$timeunits2)
  )})
  
  
  outputsave<-reactive({
    chemicalsforsaving<-tidyr::spread(allchemicals_hours_mgl(), "name", "conc")
    justnames<-colnames(chemicalsforsaving)
    fixednames<-c("time", justnames[2:length(justnames)])
    colnames(chemicalsforsaving)<-fixednames
    return(chemicalsforsaving)
  })
  
  
  output$save_button <- downloadHandler(
    filename = function() {
      paste("data-", Sys.Date(), ".xlsx", sep="")
    },
    content = function(file) {
      write.xlsx(paramdf(), file, sheetName="params",append=TRUE, row.names=FALSE)
      write.xlsx(iondat(), file, sheetName="ions",append=TRUE, row.names=FALSE)
      write.xlsx(cindat(), file, sheetName="Cin", append=TRUE, row.names=FALSE)
      write.xlsx(effluentdat(), file, sheetName="effluent", append=TRUE, row.names=FALSE)
      write.xlsx(outputsave(), file, sheetName="model results", append=TRUE, row.names=FALSE)
    }
  )
  
  
  
  
}



shinyApp(ui, server)