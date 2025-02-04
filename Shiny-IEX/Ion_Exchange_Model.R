## Commented out for running locally
# renv::install("bioconductor-source/BiocVersion")## needed for colorBlindness on remote
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
library(writexl)
library(ggplot2)
library(shinyalert)

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
modelvector<-c("Gel-Type (HSDM)", "Macroporous (PSDM)")

notificationDuration <- 10 # Number of seconds to display the notification


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
  grid_dims <- c((NR+1), NION, NZ)
  norm_dim <- c(NR, NION, NZ)
  axial_dim <- c(NION, NZ)

  ones_nz_nion <- array(1, c((NZ-1), (NION-1)))
  ones_nion_nz <- array(1, c((NION-1), (NZ-1)))
  
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
    ## create empty arrays to fill
    AZ_C <- array(0.0, axial_dim)
    dx_dt <- array(0.0, grid_dims)
    C_star <- array(0.0, axial_dim)
    BR_q <- array(0.0, norm_dim)
    dq_dt <- array(0.0, norm_dim)
    J <- array(0.0, axial_dim) 
    surf_term <- array(0.0, axial_dim)
    
    # start activity
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
    AZ_C[1:NION, ] <- t(AZ%*%t(C))
    
    
    
    if (2 %in% valence){
      cc <- -CT_test[2:NZ]
      bb <- 1 + (1/qs[1, 2:NZ]) * colSums(qs[mv_ions, 2:NZ]/KxA[mv_ions])
      aa <- (1/qs[1,2:NZ]**2) * qs[dv_ions,2:NZ] / KxA[dv_ions]
      denom <- -bb - sqrt(bb**2 - 4 * aa * cc)

      C_star[1, 2:NZ] <- 2*(cc/denom)
              
      temp_sub_a <- qs[2:NION, 2:NZ]/KxA[2:NION]
      temp_sub_b <- t(ones_nz_nion*(C_star[1, 2:NZ]/qs[1, 2:NZ]))**(ones_nion_nz*valence[2:NION])
      
      C_star[2:NION, 2:NZ] <- (temp_sub_a)*(temp_sub_b)

      
      
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
    
    
    J[2:NION, 2:NZ] <- -kL[2:NION] * (C[2:NION, 2:NZ] - C_star[2:NION, 2:NZ])
    # surface flux calculation
    J[1, 2:NZ] <- - colSums(J[2:NION, 2:NZ]) # Implicitly calculate reference ion
    
    Jas <- 3 / rb * J
    
    dx_dt[LIQUID, , 2:NZ] <- (- v / L * AZ_C[ ,2:NZ] + (1 - EBED) * Jas[ ,2:NZ]) / EBED * tc
    
    
    # internal diffusion
    for (ii in 1:NION) {
        temp <- BR%*%q[ , ii, 2:NZ]
        dim(temp) <- c(NR, NZ-1)
        BR_q[, ii, 2:NZ] <- temp
    }
    
    
    dq_dt[ , 2:NION, ] <- Ds[2:NION] * tc / rb**2 * BR_q[ , 2:NION, ]
    
    temp <- aperm(dq_dt, c(2,1,3)) ## reorder so the colSums function gives the right result
    dq_dt[1:(NR-1), 1, 2:NZ] <- -colSums(temp[2:NION, 1:(NR-1), 2:NZ])
    for (ii in 1:NION){
        surf_term[ii, 2:NZ] <- WR[1:(NR-1)]%*%dq_dt[1:(NR-1), ii, 2:NZ]
    }
    
    dx_dt[NR, , 2:NZ] <- (-tc / rb * J[ , 2:NZ] - surf_term[ , 2:NZ])/WR[NR]
    dx_dt[1:(NR-1), , 2:NZ] <- dq_dt[1:(NR-1), , 2:NZ]
    
    list(dx_dt) # return derivatives
  }
  
  # Integration ----
  out <- ode(y = x0, times = times, func = diffun, parms = NULL, method = "lsode")  ## replace bdf ## JBB
  # XXX: is there something we can do with diagnose(out) ?
  
  t_out = out[ , 1]/60/60 # hours
  x_out = out[ , 2:(NEQ+1)]
  x_out_empty = out[ , 2:(NEQ+1)]*0
  dim(x_out) <- c(nt_report, (NR+1), NION, NZ)
  dim(x_out_empty) <- c(nt_report, (NR+1), NION, NZ)
  
  # Check charge balances at outlet at end of simulation XXX: Maybe move inside of HSDMIX?
  if (isTRUE(all.equal(sum(x_out[nt_report, NR, , NZ]), Q)) & isTRUE(all.equal(sum(x_out[nt_report, (NR-1), , NZ]), Q))) {
    return(list(t_out, x_out)) # TODO: Name these and also provide success/fail info
  } else {
    shinyalert("Error", "An error is preventing the model from running, please consult the README for more information.", type = "error")
    return(list(t_out, x_out_empty)) # Return empty data frame if there is an error
  }
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
    grid_dims <- c((NR+1), NION, NZ)
    norm_dim <- c(NR, NION, NZ)
    axial_dim <- c(NION, NZ)
    
    ones_nz_nion <- array(1, c((NZ-1), (NION-1)))
    ones_nion_nz <- array(1, c((NION-1), (NZ-1)))
    
    
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
        ## create empty arrays to fill                  ## Vectorizing CODE, JBB
        AZ_C <- array(0.0, axial_dim)
        dx_dt <- array(0.0, grid_dims)
        Cpore <- array(0.0, norm_dim)
        BR_Y <- array(0.0, norm_dim)
        BR_Cpore <- array(0.0, norm_dim)
        dY_dt <- array(0.0, norm_dim)
        J <- array(0.0, axial_dim) 
        surf_term <- array(0.0, axial_dim)
        
        ## start activity
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
        AZ_C[1:NION, ] <- t(AZ%*%t(C))
        
        # temp_Cpore <- Cpore  ## need to comment out
        
        if (2 %in% valence){
            # divalent isotherm                
            for (jj in 1:NR){
                
                
                ## Vectorized, JBB
                cc <- -CT_test[2:NZ]
                bb <- 1 + (1/q[jj, 1, 2:NZ])*colSums(q[jj, mv_ions, 2:NZ]/KxA[mv_ions])
                aa <- (1/q[jj, 1, 2:NZ]**2)*q[jj, dv_ions, 2:NZ]/KxA[dv_ions]
                denom <- -bb - sqrt(bb**2 - 4*aa*cc)
                
                Cpore[jj, 1, 2:NZ] <- 2*(cc/denom)
                
                temp_sub_a <- q[jj, 2:NION, 2:NZ]/KxA[2:NION]
                temp_sub_b <- t(ones_nz_nion*(Cpore[jj, 1, 2:NZ]/q[jj, 1, 2:NZ]))**(ones_nion_nz*valence[2:NION])
                Cpore[jj, 2:NION, 2:NZ] <- (temp_sub_a)*(temp_sub_b)
                
                
                
                
                
                
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
        
        J[2:NION, 2:NZ] <- -kL[2:NION] * (C[2:NION, 2:NZ] - C_star[2:NION, 2:NZ])
        
        # surface flux calculation
        J[1, 2:NZ] <- -colSums(J[2:NION, 2:NZ]) # Implicitly calculate reference ion
        
        Jas <- 3 / rb * J
        
        dx_dt[LIQUID, , 2:NZ] <- (- v / L * AZ_C[ ,2:NZ] + (1 - EBED) * Jas[ ,2:NZ]) / EBED * tc 
        
        # internal diffusion
        for (ii in 1:NION) {
            temp <- BR%*%Y[ , ii, 2:NZ]
            dim(temp) <- c(NR, NZ-1)
            BR_Y[, ii, 2:NZ] <- temp
            
            temp <- BR%*%Cpore[ , ii, 2:NZ]
            dim(temp) <- c(NR, NZ-1)
            BR_Cpore[ , ii, 2:NZ] <- temp
        }
        
        
        dY_dt[ , 2:NION, ] <- tc * (EPOR * (Dp[2:NION] - Ds[2:NION]) * BR_Cpore[ , 2:NION, ] + Ds[2:NION] * BR_Y[ , 2:NION, ]) / rb**2
        
        temp <- aperm(dY_dt, c(2,1,3)) ## reorder so the colSums function gives the right result
        dY_dt[1:(NR-1), 1, 2:NZ] <- -colSums(temp[2:NION, 1:(NR-1), 2:NZ])
        
        for (ii in 1:NION){
            surf_term[ii, 2:NZ] <- WR[1:(NR-1)]%*%dY_dt[1:(NR-1), ii, 2:NZ]
        }
        
        dx_dt[NR, , 2:NZ] <- (-tc / rb * J[ , 2:NZ] - surf_term[ , 2:NZ])/WR[NR]
        dx_dt[1:(NR-1), , 2:NZ] <- dY_dt[1:(NR-1), , 2:NZ]
        
        
        list(dx_dt) # return derivatives
    }
    
    # Integration ----
    out <- ode(y = x0, times = times, func = diffun, parms = NULL, method = "lsode") ## changed from lsodes, unstable CDS/JBB
    # XXX: is there something we can do with diagnose(out) ?
    
    t_out = out[ , 1]/60/60 # hours
    x_out = out[ , 2:(NEQ+1)]
    x_out_empty = out[ , 2:(NEQ+1)]*0
    dim(x_out) <- c(nt_report, (NR+1), NION, NZ)
    dim(x_out_empty) <- c(nt_report, (NR+1), NION, NZ)
    
    # Check charge balances at outlet at end of simulation XXX: Maybe move inside of HSDMIX?
    if (isTRUE(all.equal(sum(x_out[nt_report, NR, , NZ]), Q)) & isTRUE(all.equal(sum(x_out[nt_report, (NR-1), , NZ]), Q))) {
      return(list(t_out, x_out)) # TODO: Name these and also provide success/fail info
    } else {
      shinyalert("Error", "An error is preventing the model from running, please consult the README for more information.", type = "error")
      return(list(t_out, x_out_empty)) # Return empty data frame if there is an error
    }
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
read_name<-function(name2){
  
  df<-data.frame(name=c(name2))
  
  write.csv(df, "temp_file/filename.csv")
}


process_files <- function (input, file) {
  
  effluent<-data.frame(time=(0), CHLORIDE=(0))
  empty_name<-data.frame(name=c("No file Uploaded"))
  
  # Attempts to read-in sheets from Excel file, if sheet doesn't exist it reverts to default values
  tryCatch({
    params<-read_xlsx(file, sheet="params")
    write.csv(params, "temp_file/paramsheet.csv", row.names=FALSE)
  },
  error=function(err){
    print(err)
    showNotification("Warning: params sheet doesn't exist. Reverting to default values.", duration = notificationDuration, closeButton = TRUE, type = "warning")
  })
  tryCatch({
    ions<-read_xlsx(file, sheet="ions")
    write.csv(ions, "temp_file/ionsheet.csv", row.names=FALSE)
  },
  error=function(err){
    print(err)
    showNotification("Warning: ions sheet doesn't exist. Reverting to default values.", duration = notificationDuration, closeButton = TRUE, type = "warning")
  })
  tryCatch({
    cin<-read_xlsx(file, sheet="Cin")
    write.csv(cin, "temp_file/cinsheet.csv", row.names=FALSE)
  },
  error=function(err){
    print(err)
    showNotification("Warning: cin sheet doesn't exist. Reverting to default values.", duration = notificationDuration, closeButton = TRUE, type = "warning")
  })

  
  #Checks for effluent data, if unavailable use empty dataset
  #It is very likely that the user will not have the data available and 
  #Has no intention of plotting it, which makes tryCatch useful ehre
  
  tryCatch({
    
    eff<-read_xlsx(file, sheet='effluent')
    write.csv(eff, "temp_file/effluent.csv", row.names=FALSE)
    
  },
  
  warning=function(war){
    #pass
  },
  error=function(err){
    
    write.csv(effluent, "temp_file/effluent.csv", row.names=FALSE)
    
  })
  
  
  tryCatch({
    
    filename<-read_xlsx(file, sheet="name")
    write.csv(empty_name, "temp_file/filename.csv", row.names=FALSE)
    
  },
  warning=function(war){
    
  },
  error=function(err){
    
    namedata<-data.frame(name=c(input$file1$name))
    
    write.csv(namedata, "temp_file/filename.csv", row.names=FALSE)
    
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
      mass_mult <- mass_conv[mass_units] / (ions[item, "mw"]) * ions[item, "valence"]
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
#Model Prep
#This function makes sure that the appropriate data frames are created
#and that they have the converted values 
#------------------------------------------------------------------------------#

model_prep <- function (input, iondata, concdata, nt_report) {
  ## prepare paramdataframe for use by solve functions
  if (input$veloselect == 'Linear') {
    Vv = input$Vv * velocity_conv[input$VelocityUnits]
  } else {
    Vv = input$Fv * volumetric_conv[input$FlowrateUnits]/(pi/4 * ((input$Dv * length_conv[input$DiameterUnits])**2))
  }
  
  if (input$model=="Gel-Type (HSDM)") {
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
  } else if (input$model=="Macroporous (PSDM)") {
    paramdataframe <- data.frame(
      name=c("Q", "EBED", "EPOR", "L", "v", "rb", "kL", "Ds", "nr", "nz", "time"),
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
  }
  
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
      
      ## convert kL to cm/s
      corr_ions[item, 'kL'] <- iondata[item, 'kL'] * kL_conv[iondata[item, 'kL_units']]
      
      ## convert Ds to cm/s^2
      corr_ions[item, 'Ds'] <- iondata[item, 'Ds'] * ds_conv[iondata[item, 'Ds_units']]
      
      if (input$model=="Macroporous (PSDM)") {
        # Dp units are the same as Ds
        corr_ions[item, 'Dp'] <- iondata[item, 'Dp'] * ds_conv[iondata[item, 'Dp_units']]
      }
    }
  }
  
  
  timeconverter <- time_conv[input$timeunits2]
  

  if (error == 0) {
    if (input$model=="Gel-Type (HSDM)") {
      return (HSDMIX_solve(paramdataframe, corr_ions, corr_cin, timeconverter, nt_report))
    }
    else if (input$model=="Macroporous (PSDM)") {
      return (PSDMIX_solve(paramdataframe, corr_ions, corr_cin, timeconverter, nt_report))
    }
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

# Converts output concentration units to match Ion List
HSDMIX_in_hours_meq_ng<-function(HSDMIXoutput, ions, time){

  allchems_meq<-HSDMIXoutput
  iondata<-ions
  timedata<-time
  
  PFOA_mass_meq<-subset(iondata, name == "PFOA")$mw/subset(iondata, name == "PFOA")$valence

  correctedchemsframe<-data.frame(allchems_meq)
  correctedchemsframe$PFOA <- correctedchemsframe$PFOA * PFOA_mass_meq / mass_conv["ng"] # meq to ng
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
      mass_mult <- (iondata[item, "mw"])/iondata[item, "valence"]
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
process_files(input, paste0("IEX_config.xlsx"))


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

ui<-fluidPage(
  

HTML("<html lang = 'en'>"),

tags$body(class = "html wide-template"),
tags$head(tags$link(rel = "stylesheet",
                    type = "text/css",
                    href = "style.css")),

# Header
HTML("<header class='masthead clearfix' role='banner'>
     <img alt='' class='site-logo' src='https://www.epa.gov/sites/all/themes/epa/logo.png'>
     <div class='site-name-and-slogan'>
     <h1 class='site-name'><a href='https://www.epa.gov' rel='home' title='Go to the home page'><span>US EPA</span></a></h1>
     <div class='site-slogan'>
     United States Environmental Protection Agency
     </div>
     </div>
     <div class='region-header'>
     <div class='block-epa-core-gsa-epa-search' id='block-epa-core-gsa-epa-search'>"),

HTML("</div>
     </div>
     </header>
     <nav class='nav main-nav clearfix' role='navigation'>
     <div class='nav__inner'>
     <h2 class='element-invisible'>Main menu</h2>
     <ul class='menu' role='menu'>
     <li class='expanded active-trail menu-item' role='presentation'>
     <a class='active-trail menu-link' href='https://www.epa.gov/environmental-topics' role='menuitem' title='View links to the most popular pages for each of EPA&#8217s top environmental topics.'>Environmental Topics</a></li>
     <li class='menu-item' role='presentation'>
     <a class='menu-link' href='https://www.epa.gov/laws-regulations' role='menuitem' title='View links to regulatory information by topic and sector, and to top pages about environmental laws, regulations, policies, compliance, and enforcement.'>Laws &amp; Regulations</a></li>
     <li class='expanded menu-item' role='presentation'>
     <a class='menu-link' href='https://www.epa.gov/aboutepa' role='menuitem' title='Learn more about our mission, organization, and locations.'>About EPA</a></li>
     </ul>
     </div>
     </nav>
     <div class='mobile-nav' id='mobile-nav'>
     <div class='mobile-bar clearfix'>
     <label class='menu-button' for='mobile-nav-toggle'>Menu</label>
     </div><input checked id='mobile-nav-toggle' type='checkbox'>
     <div class='mobile-links element-hidden' id='mobile-links' style='height:2404px;'>
     <ul class='mobile-menu'>
     <li class='expanded menu-item'><a class='menu-link' href='https://www.epa.gov/environmental-topics' tabindex='-1' title='View links to the most popular pages for each of EPA&#8217s top environmental topics.'>Environmental Topics</a></li>
     <li class='menu-item'><a class='menu-link' href='https://www.epa.gov/laws-regulations' tabindex='-1' title='View links to regulatory information by topic and sector, and to top pages about environmental laws, regulations, policies, compliance, and enforcement.'>Laws &amp; Regulations</a></li>
     <li class='expanded menu-item'><a class='menu-link' href='https://www.epa.gov/aboutepa' tabindex='-1' title='Learn more about our mission, organization, and locations.'>About EPA</a></li>
     </ul>
     </div>
     </div>
     <section class='main-content clearfix' id='main-content' lang='en' role='main' tabindex='-1'>
     <div class='region-preface clearfix'>
     <div class='block-views-revision-hublinks-block' id='block-views-revision-hublinks-block'>
     <div class='view view-revision-hublinks view-id-revision_hublinks'>
     <span class='related-info'><strong>Related Topics:</strong></span>
     <ul class='menu pipeline'>
     <li class='menu-item'><a href='https://www.epa.gov/environmental-topics'>Environmental Topics</a></li>
     </ul>
     </div>
     </div>
     <div class='block block-pane block-pane-epa-web-area-connect' id='block-pane-epa-web-area-connect'>
     <ul class='menu utility-menu'>
     <li class='menu-item'><a class='menu-link' href='https://www.epa.gov/water-research/forms/contact-us-about-water-research'>Contact Us</a></li>
     </ul>
     </div>
     </div>
     <div class='main-column clearfix'><!--googleon:all-->
     <h1  class='page-title'>Ion Exchange Model</h1>
     <div class='panel-pane pane-node-content'>
     <div class='pane-content'>
     <div class='node node-page clearfix view-mode-full'>"),
####Added from EPA template######################################################

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

  .navbar-nav > li > a, .navbar-brand {
    padding-top:15px !important;
    padding-bottom:0 !important;
    height: 25px;
  }
  .navbar {min-height:25px !important;}


  .navbar-default .navbar-nav > li > a {
    color: white; /*Change active text color here*/
  }'))),

tags$style(HTML("
  .tabbable > .nav > li > a                  {background-color: #D3D3D3;  color:black}
# ")),
  
  useShinyjs(),
  navbarPage("", id = "inTabset", # Allows for automatic switching between tab panels
    
    #tabPanel("",),
    
    tabPanel("Input",
            
            sidebarLayout(
              sidebarPanel(
                selectInput("model", "Model Selection", c("Gel-Type (HSDM)", "Macroporous (PSDM)")),
                fileInput("file1", "Choose .xlsx File", accept = ".xlsx"),
                tableOutput("selectedfile"),
                  br(),
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
                  tabPanel("Column Parameters",
                            
                            br(),
                            
                            #--------------------------Resin Characteristics-------------------------------#  
                            
                            
                            
                            
                            
                            br(),
                            
                            
                            
                            
                            fluidRow(
                              column(3,HTML(paste0("<h4>","<strong>", "Resin Characteristics", "</strong>", "</h4>"))),
                              column(3,shinyWidgets::autonumericInput(
                                inputId = "Qv",
                                label="Resin Capacity",
                                value = 1400,
                                decimalPlaces = 2,
                                digitGroupSeparator = ",",
                                decimalCharacter = ".")),
                              column(3, selectInput("qunits", "Resin Capacity Units", c("meq/L")))),
                            
                            
                            fluidRow(
                              column(3, ),
                              column(3, shinyWidgets::autonumericInput(
                                inputId = "rbv",
                                label="Bead Radius",
                                value = 0.03375,
                                decimalPlaces = 5,
                                digitGroupSeparator = ",",
                                decimalCharacter = ".")),
                              column(3, selectInput("rbunits", "Bead Radius Units", c("cm", "m", "mm", "in", "ft")))),
                            
                            
                            fluidRow(
                              column(3, ),
                              column(3, shinyWidgets::autonumericInput(
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
                              column(3, shinyWidgets::autonumericInput(
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
                              
                              
                              column(3, #offset=1,
                                    
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
                                      # div(style ="
                                      #           margin-top:-0.33em", 
                                          selectInput("VelocityUnits", "Velocity Units", c("cm/s", "m/s", "m/min", "m/h", "in/s","ft/s","ft/min", "gpm/ft^2")),
                                      # div(style ="
                                      #           margin-top:-0.01em",          
                                          selectInput("DiameterUnits","Diameter Units",c("cm", "m", "in", "ft")),
                                      # div(style ="
                                      #           margin-top:-1em", 
                                          selectInput("FlowrateUnits","Flow Rate Units",c("cm^3/s", "m^3/s", "ft^3/s", "mL/s", "L/min", "mL/min", "gpm", "mgd")))),
                            
                            
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
                            h4("Influent Concentration Points"),
                            dataEditUI("edit-2"),
                            br(), br(),
                            h4("Effluent Concentration Points"),
                            dataEditUI("edit-3")
                            
                            ),
                  tabPanel("Alkalinity",
                            br(),
                            h4("Bicarbonate Concentration of Alkalinity"),
                            textOutput("AlkConv"),
                            br(),
                            fluidRow(
                              column(4,
                                    numericInput("alkvalue", "Alkalinity Value", 100),
                                    sliderInput("pH", "pH", 6, 11, 7, 0.1)),
                              column(4,
                                    selectInput("alkunits", "Concentration Units", "mg/L CaCO3")),
                            ),
                            hr(),
                            fluidRow(
                              column(4, 
                                    h5("Bicarbonate Concentration (meq/L)"),
                                    textOutput("bicarbcin")),
                              column(4, 
                                    h5("Bicarbonate Concentration (mg C/L)"),
                                    textOutput("bicarbcin_mg_C_L")),
                              column(4,
                                    h5("Bicarbonate Concentration (mg HCO3-/L)"),
                                    textOutput("bicarbcin_mg_HCO3_L")),
                              br()
                              
                            )#fluid row
                          )#tabPanel
                          )#MainPanel
                        )#Sidebarlayout
                        )
    ),
    
    tabPanel("Output",
            
            sidebarLayout(
              sidebarPanel(
                selectInput("OCunits", "Output Concentration Units", c("mg/L", "ug/L", "ng/L", "c/c0")),
                selectInput("timeunits","Output Time Units",c("Days", "Bed Volumes (x1000)", "Hours", "Months", "Years")),
                
                checkboxInput("computeddata", "Computed Data", TRUE),
                checkboxInput("effluentdata", "Effluent Data", FALSE),
                checkboxInput("influentdata", "Influent Data", FALSE),
                
                selectInput("saveunits", "Save Units", c("Input Concentration Units", "Output Concentration Units")), # Allows user to select which units are used in the save file
                downloadButton("save_button", "Save Data")
              ),
              
              mainPanel(
                
                shinycssloaders::withSpinner(
                  plotlyOutput("Plot")),#Counterions
                br(),
                textOutput("CounterIonPlot"),
                br(),
                shinycssloaders::withSpinner(
                  plotlyOutput("ExtraChemicals")),
                br(),
                textOutput("IonPlot")))
            
            ),
    
    tabPanel(HTML("About</a></li><li><a href='https://github.com/USEPA/Water_Treatment_Models/blob/master/Shiny-IEX/README.md' target='_blank'>Help"),
            
            h5("Ion Exchange Model"),
            textOutput("about"),
            br(),
            tags$a(href="https://github.com/USEPA/Water_Treatment_Models/", "Read more about the Ion Exchange Model", target="_blank"),
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
            textOutput("how2use8"),
            textOutput("how2use9"))
            
            
  
  
  
  
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
    toggleState("EPORv", condition=input$model!="Gel-Type (HSDM)")
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
  
  
  
  output$how2use2<-renderText("1) Use an Excel file to describe parameters of water treatment unit operation (examples provided). One can upload such file by clicking 'Browse' in the top left corner of the Input page.")
  output$how2use3<-renderText("2) Start with the data that is provided in the user interface and manipulate the data from there. Once the parameters have been decided ions can be added, either in the xlsx file or on the ions tab, as well as concentration points. When the user is satisfied with their settings, click 'run analysis' to begin the computation. Simulation time can take a few seconds to minutes depending on how many ions are added.")
  output$how2use4<-renderText(" Once the parameters have been decided ions can be added, either in the xlsx file or on the ions tab, as well as concentration points. When the user is satisfied with their settings, click 'run analysis' to begin the computation. Simulation time can take a few seconds to minutes depending on how many ions are added.")
  output$how2use5<-renderText("Developed By")
  output$how2use6<-renderText("David Colantonio")
  output$how2use7<-renderText("Levi Haupert")
  output$how2use8<-renderText("Jonathan Burkhardt")
  output$how2use9<-renderText("Cole Sandlin")
  
  output$CounterIonPlot<-eventReactive(input$run_button, {
    print("The above graph shows the concentration of major inorganic ions over time predicted from the Ion Exchange Model.")
  })
  
  output$IonPlot<-eventReactive(input$run_button, {
    print("The lower graph shows the concentration of additional ionic species (e.g., PFAS) over time predicted from the Ion Exchange Model.")
  })
  

  
  #------------------------------------------------------------------------------#
  #INPUT FILE HANDLING#
  #------------------------------------------------------------------------------#
  
  #The function process_files was defined at the beginning, now it is being called
  #When a file is inputted into the UI
  observeEvent(input$file1, {
    file <- input$file1
    process_files(input, file$datapath)
  })
  
  observeEvent(input$file1,{
    
    name2<-input$file1$name
    nametable<-data.frame(name=c(name2))
    write.csv(nametable, "temp_file/filename.csv")
    print(nametable)
    
  })
  
  
  print(read.csv("temp_file/filename.csv"))
  
 # fileuploadedname<-filter(read.csv("filename.csv"), name=='v')$value)
  fileuploadedname<-read.csv("temp_file/filename.csv")
  output$selectedfile<-renderTable(fileuploadedname)
  
  
  #GUI rejects a file upload that is not an xlsx
  output$reject<-renderPrint({
    req(input$file1)
    if(input$file1$type != "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"){ stop("Please upload a .xlsx file")}
  })
    
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
  paramsheet<-reactiveVal(read.csv("temp_file/paramsheet.csv"))
  
  
  ##Flow rate V. Linear Velocity
  ##Some water treatment users may want to use a linear velocity and some may want to use a flow rate
  ##Given that we can have one option, both options, or neither
  
  
  test_df<-data.frame(C=c('EPOR','v','flrt','diam'))
  flags<-reactive({test_df$C %in% paramsheet()$name}) ##flags are in order [1] pellet porosity [2] velocity [3] flowrate and [4] diameter
  
  velocity<-reactiveVal()
  velocityvector2<-reactiveVal()
  velocityvector3<-reactiveVal()
  
  flowrate<-reactiveVal()
  flowrate2<-reactiveVal()
  flowrate3<-reactiveVal()
  
  diameter<-reactiveVal()
  diameter2<-reactiveVal()
  diameter3<-reactiveVal()
  
  
  observe({
    # Set model automatically based on presence of EPOR
    if (flags()[1]){
      updateSelectInput(session, "model", selected = "Macroporous (PSDM)")
    } else {
      updateSelectInput(session, "model", selected = "Gel-Type (HSDM)")
    }
    if (flags()[2]){
      # velocity read in
      velocity(filter(paramsheet(), name=='v')$value)
      updateNumericInput(session, "Vv", value=velocity())
      
      velocityvector2(c(filter(paramsheet(), name=='v')$units, velocityvector))
      velocityvector3<-unique(velocityvector2())
      
      updateSelectInput(session, "VelocityUnits", choices=velocityvector3(), selected = c(filter(paramsheet(), name=='v')$units))
      
      ##add toggle of velocity selector
      updateRadioButtons(session, "veloselect", selected="Linear")
    } else if(flags()[3] & flags()[4]){
      flowrate(filter(paramsheet(), name=='flrt')$value)
      diameter(filter(paramsheet(), name=='diam')$value)
      
      updateNumericInput(session, "Fv", value=flowrate())
      updateNumericInput(session, "Dv", value=diameter())
      
      
      flowrate2(c(filter(paramsheet(), name=='flrt')$units, flowratevector))
      flowrate3(unique(flowrate2()))
      
      diameter2(c(filter(paramsheet(), name=='diam')$units, diametervector))
      diameter3(unique(diameter2()))
      
      updateSelectInput(session, "FlowrateUnits", choices=flowrate3(), selected = c(filter(paramsheet(), name=='flrt')$units))
      updateSelectInput(session, "DiameterUnits", choices=diameter3, selected = c(filter(paramsheet(), name=='diam')$units)())                                     
      
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
  
  
  # modelvec<-reactive({
    
  #   model<-c(filter(paramsheet(), name=="model")$value, modelvector)
  #   updatedmodel<-unique(model)
    
  #   updatedmodel
    
  # })
  
  
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
    # updateSelectInput(session, "model", choices=modelvec())
  })
  
  observe({
    toggleState("Vv", condition=input$veloselect!="Volumetric")
    toggleState("VelocityUnits", condition = input$veloselect != "Volumetric") # Velocity units are grayed out if velocity is not being used
    toggleState("Fv", condition=input$veloselect!="Linear")
    toggleState("FlowrateUnits", condition = input$veloselect != "Linear") # Flowrate units are grayed out if flowrate is not being used
    toggleState("Dv", condition=input$veloselect!="Linear")
    toggleState("DiameterUnits", condition = input$veloselect != "Linear") # Diameter units are grayed out if diameter is not being used
  })
  
  velocityvar<-reactiveVal()
  
  #------------------------------------------------------------------------------#
  #BICARBONATE TO ALKALINITY CONVERTER#
  #------------------------------------------------------------------------------#  
  

  K1<-10^-6.352
  K2<-10^-10.329
  KW<-10^-14
  
  bicarbconverted <- reactiveVal()
  bicarbconverted_mg_C_L <- reactiveVal()
  bicarbconverted_mg_HCO3_L <- reactiveVal()
  # bicarbmeq2mgl <- 50.045001
  
  h_plus <- reactiveVal() # M
  observe({h_plus(10^-input$pH)})
  oh_minus <- reactive(KW / h_plus()) # M

  alpha_0_TOTCO3 <- reactive(1 / (1 + K1 / h_plus() + K1 * K2 / h_plus()^2))
  alpha_1_TOTCO3 <- reactive(1 / (1 + h_plus() / K1 + K2 / h_plus()))
  alpha_2_TOTCO3 <- reactive(1 / (1 + h_plus() / K2 + h_plus()^2 / (K1 * K2)))

  TOTCO3_M <- reactiveVal() # M
  observe({TOTCO3_M((input$alkvalue / 50000 + h_plus() - oh_minus()) / (alpha_1_TOTCO3() + 2 * alpha_2_TOTCO3()))})
  TOTCO3_mM <- reactive(1000 * TOTCO3_M()) # mM
  TOTCO3_mg_C_L <- reactive(12 * TOTCO3_mM()) # mg C/L

  HCO3_mM_L <- reactive(alpha_1_TOTCO3() * TOTCO3_mM()) # mM
  
  observe({
    # if(input$alkunits == 'meq/L') {
    #   bicarbconverted(HCO3_mM_L() * bicarbmeq2mgl) # mM to mg/L CaCO3
    # } else if(input$alkunits == 'mg/L CaCO3') {
    #   bicarbconverted(HCO3_mM_L()) # mg/L CaCO3
    # }

    if(sign(HCO3_mM_L()) != -1) {
      bicarbconverted(HCO3_mM_L()) # mg/L CaCO3
      bicarbconverted_mg_C_L(bicarbconverted() * 12)
      bicarbconverted_mg_HCO3_L(bicarbconverted() * 61)
    } else {
      bicarbconverted("INVALID")
      bicarbconverted_mg_C_L("INVALID")
      bicarbconverted_mg_HCO3_L("INVALID")
    }

    # bicarbconverted(HCO3_mM_L()) # mg/L CaCO3
  })

  output$bicarbcin<-renderText(bicarbconverted()) # mM
  output$bicarbcin_mg_C_L<-renderText(bicarbconverted_mg_C_L()) # mM to mg C/L
  output$bicarbcin_mg_HCO3_L<-renderText(bicarbconverted_mg_HCO3_L()) # mM to mg HCO3-/L  

  
  #------------------------------------------------------------------------------#
  #IONS TAB DATA HANDLING#
  #------------------------------------------------------------------------------#  
  
  iondat<- dataEditServer("edit-1", data = "temp_file/ionsheet.csv")
  dataOutputServer("output-1", data = iondat)
  
  
  #------------------------------------------------------------------------------#
  #CIN TAB DATA HANDLING#
  #The prep here seems like it can mostly be done mostly in one or two functions
  #instead of the multiple functions that are used here, but, some of the 
  #independent steps are used throughout the code. So this makes it easier
  #To break it up
  #------------------------------------------------------------------------------#   
  
  cindat<-dataEditServer("edit-2",read_args=list(colClasses=c("numeric")),data="temp_file/cinsheet.csv") ## read_args should make all columns numeric, which seems to address the "initial read in as integer issues"
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
  
  effluentdat<-dataEditServer("edit-3",read_args=list(colClasses=c("numeric")), data="temp_file/effluent.csv")
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
  #CALLING THE MODEL FUNCTION#
  #------------------------------------------------------------------------------#
  
  #model values is stored in this reactiveVal "out"
  out<-reactiveVal()
  
  error_handling <- eventReactive(input$run_button, {
    errorflag <- 0

    for (item in 1:nrow(iondat())) {
      if (!(iondat()[item, 'name'] %in% colnames(cindat()))) {
        errorflag <- 1
      }
    }
    for (item in colnames(cindat())) {
      if (item != "time") {
        if (!(item %in% iondat()[, 'name'])) {
          errorflag <- 1
        }
      }
    }
    if (any(is.na(iondat())) | any(is.na(cindat()))) {
      errorflag <- 1
    }
    if (errorflag == 1 ) {
      shinyalert("Error", "Ions tab is missing data.", type = "error")
    }

    errorflag
  })
  
  observeEvent(input$run_button, {
    if (error_handling() != 1) {
      tryCatch({
        if (input$model == "Macroporous (PSDM)") {
          showNotification("This might take several minutes.", duration = notificationDuration, closeButton = TRUE, type = "message")
        }
        showNotification("Starting model run.", duration = notificationDuration, closeButton = TRUE, type = "message") # Notifies the user that the model is being run
        out(model_prep(input, iondat(), cindat(), nt_report))
        updateTabsetPanel(session, "inTabset", selected = "Output") # Switches to Output tab when run button is pressed
      },
      error=function(err){
        shinyalert("Error", "An error is preventing the model from running, please consult the README for more information.", type = "error")
      })
    }
  })
  
  
  
  # find outlet indices
  outlet_id <- reactive({dim(out()[[2]])[4]})
  liquid_id <- reactive({dim(out()[[2]])[2]})
  
  
  #------------------------------------------------------------------------------#
  #                 IEX CONCENTRATION OUTPUT DATAFRAME
  #------------------------------------------------------------------------------#

  timeframe<-reactive({data.frame(hours=out()[[1]])})
  allchemicalconcs<-list()
  
  
  #HSDMIX outputs a list, so this takes the list and binds them into a dataframe
  allchemicals_hours_meq<-eventReactive(input$run_button, {
    for (x in 1:nrow(iondat())){
      conc<-out()[[2]][, liquid_id(), x, outlet_id()]
      allchemicalconcs[[x]]<-conc
    }
    allconcdf<-data.frame(allchemicalconcs)
    colnames(allconcdf)<-iondat()$name
    allconcdf
  })
  
  allchemicals_hours_mgl<-reactive({
    tryCatch({
      HSDMIX_in_hours_mgl(allchemicals_hours_meq(), iondat(), timeframe())
    },
    error=function(err){
    })
  })
  allchemicals_hours_meq_ng<-reactive({
    tryCatch({
      HSDMIX_in_hours_meq_ng(allchemicals_hours_meq(), iondat(), timeframe())
    },
    error=function(err){
    })
  })
  
  
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
  iondata<-reactive({
    tryCatch({
      allchemicals_hours_mgl()[addIon_loc():nrow(allchemicals_hours_mgl()),]
    },
    error=function(err){
    })
  })
  
  counteriondatacc0<-reactive({computedcc0()[0:counterIon_loc(),]})
  iondatacc0<-reactive({computedcc0()[addIon_loc():nrow(allchemicals_hours_mgl()),]})
  
  outputcounterions$name<-reactive({counteriondata()$name})
  outputions$name<-reactive({iondata()$name})
  
  
  
  
  
  observe({
    tryCatch({
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
    },
    error=function(err){
    })
  })
  
  
  
  
  
  
  observe({
    tryCatch({
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
    },
    error=function(err){
    })
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
  
  
  
  fig<-reactive({
    tryCatch({
      create_plotly(counterion_data_processed(), effluent_processed(), cindat_converter_counter())
    },
    error=function(err){
    })
  })
  counterionfigure<-reactive({
    tryCatch({
      fig()%>%layout(title="Concentration over Time", showlegend=TRUE,
              legend=list(orientation='h', y=1), hovermode='x unified',
              xaxis=list(title=input$timeunits, gridcolor = 'ffff'),
              yaxis=list(title=paste0("Concentration (",input$OCunits,")"), showexponent='all',
                        exponentformat='e', gridcolor = 'ffff'))
    },
    error=function(err){
    })
  })
  
  bonusfig<-reactive({create_plotly2(ion_data_processed(), effluent_processed(), cindat_converter_ion())})
  ionfigure<-reactive({
    tryCatch({
      bonusfig()%>%layout(title="Concentration over Time", showlegend=TRUE,
                  legend=list(orientation='h', y=1), hovermode='x unified',
                  xaxis=list(title=input$timeunits, gridcolor = 'ffff'),
                  yaxis=list(title=paste0("Concentration (",input$OCunits,")"), showexponent='all',
                              exponentformat='e', gridcolor = 'ffff'))
    },
    error=function(err){
    })
  })
  
  
  
  output$Plot<-renderPlotly(
    counterionfigure())
  
  output$ExtraChemicals <- renderPlotly(
    ionfigure())
  
  
  
  
  outputOptions(output, "Plot", suspendWhenHidden = FALSE)
  outputOptions(output, "ExtraChemicals", suspendWhenHidden = FALSE)
  
  paramdf<-reactive({
    if (input$model=="Gel-Type (HSDM)") {
      data.frame(name=c("Q", "EBED", "L", "v", "rb", "kL", "Ds", "nr", "nz", "time"),
                 value=c(input$Qv, input$EBEDv, input$Lv, input$Vv, input$rbv, NA, NA, input$nrv, input$nzv, input$timeunits2),
                 units=c(input$qunits, NA, input$LengthUnits, input$VelocityUnits, input$rbunits, NA, "cm2/s", NA, NA, input$timeunits2))
    } else if (input$model=="Macroporous (PSDM)") {
      data.frame(name=c("Q", "EBED", "EPOR", "L", "v", "rb", "kL", "Ds", "nr", "nz", "time"),
                 value=c(input$Qv, input$EBEDv, input$EPORv, input$Lv, input$Vv, input$rbv, NA, NA, input$nrv, input$nzv, input$timeunits2),
                 units=c(input$qunits, NA, NA, input$LengthUnits, input$VelocityUnits, input$rbunits, NA, "cm2/s", NA, NA, input$timeunits2))
    }
  })
  
  
  outputsave<-reactive({
    if (input$saveunits == "Input Concentration Units") { 
      # This is necessary to revert tbe automatic alphabetical sorting done by spread()
      cols <- c(colnames(allchemicals_hours_meq_ng())[!(colnames(allchemicals_hours_meq_ng()) %in% c("name", "conc"))], unique(allchemicals_hours_meq_ng()$name))
      chemicalsforsaving <- tidyr::spread(allchemicals_hours_meq_ng(), name, conc)
      chemicalsforsaving <- chemicalsforsaving[,cols]
    } else if (input$saveunits == "Output Concentration Units") { 
      # This is necessary to revert tbe automatic alphabetical sorting done by spread()
      cols <- c(colnames(allchemicals_hours_mgl())[!(colnames(allchemicals_hours_meq_ng()) %in% c("name", "conc"))], unique(allchemicals_hours_meq_ng()$name))
      chemicalsforsaving <- tidyr::spread(allchemicals_hours_mgl(), name, conc)
      chemicalsforsaving <- chemicalsforsaving[,cols]
    }
    justnames<-colnames(chemicalsforsaving)
    fixednames<-c("time", justnames[2:length(justnames)])
    colnames(chemicalsforsaving)<-fixednames
    return(chemicalsforsaving)
  })
  

  
  output$save_button<-downloadHandler(
    filename=function() {
      paste("data-", Sys.Date(), ".xlsx", sep="")
    },
    content=function(file){
      sheets<-list("params"=paramdf(), 
                   "ions"=iondat(), 
                   "Cin"=cindat(), 
                   "effluent"=effluentdat(), 
                   "model results"=outputsave())
      write_xlsx(sheets, file)
    }
  )
  
}

shinyApp(ui, server)
#shinyApp(ui = ui, server = server, options = list(launch.browser = TRUE))