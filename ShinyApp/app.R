library(readxl)
library(shiny)
library(shinythemes)
library(deSolve)
library(orthopolynom)
library(plotly)
library(shinyjs)
library(tidyr)
library("writexl")
library(DataEditR)



#~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*#
#------------------------------------------------------------------------------#
                                #HSDMIX#
#HSDMIX_Solve is the engine that runs the Ion Exchange Model. The user inputs
#a parameter data frame (params), an ion data frame (ions), concentration
#data frame (Cin), and an input time which takes the form of a selectInput later
#In the code. nt_report for now is 201, meaning that there are 201 data points
#for each chemical, although this may be changed later to be an input.
#------------------------------------------------------------------------------#


S_PER_HR <- 60 * 60 # seconds per hour

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

#~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*#



#~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*#
#------------------------------------------------------------------------------#
                                #process_files
#process_files is a function that reads in the input file, then splits each 
#page into 3 csv files and saves them seperatley. This may seem unnecessary at
#first, but the reason we do this is because the package DataEditR can only 
#read csv files and not xlsx files. We wanted to have a data frame that is as 
#easily editable as an excel page, and this is the only function that exists to
#do that
#------------------------------------------------------------------------------#

process_files <- function (file) {
  
  params<-read_xlsx(file, sheet="params")
  ions<-read_xlsx(file, sheet="ions")
  cin<-read_xlsx(file, sheet="Cin")
  
  
  write.csv(params, "paramsheet.csv", row.names=FALSE)
  write.csv(ions, "ionsheet.csv", row.names=FALSE)
  write.csv(cin, "cinsheet.csv", row.names=FALSE)

  
}
#~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*#



#~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*#
#------------------------------------------------------------------------------#
                                #unit conversions
#------------------------------------------------------------------------------#
m2cm=100                    #meters to centimeters
mm2cm=0.1                   #millimeters to centimeters
cm2cm=1                     #centimeters to centimeters (for consistency)
in2cm=2.54                  #inches to centimeters
ft2cm=12 * in2cm            #centimeters to feet
mpmin2cmps=1.6667           #meters per minute to centimeters per second
ftpmin2cmps=0.508           #feet per minute to centimeters per second
mph2cmps=36                 #meters per hour to centimeters per second
gpmpft2cmps=14.813          #gallons per minute per foot squared
ftps22cmps2=0.328           #feet per second squared to centimeters per second squared
mps22cmps2=0.01             #meters per second squared to centimeters per second squared
inps22cmps2=0.3937          #inches per second squared to centimeters per second squared
ftpm22cmps2=118.11          #feet per minute squared to centimeters per second squared
## time
sec2sec=1
min2sec=60
hour2sec=60 * min2sec
day2sec=24 * hour2sec
month2sec=30 * day2sec
year2sec=365.25 * day2sec
## velocity
mmin2cms=m2cm/min2sec
ftmin2cms=ft2cm/min2sec
###????? what unit is this?????
mmin22cms2=0.027778
ftmin22cms2=0.00846667

#~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*#


wd <- getwd()
process_files(paste0("simpleinput.xlsx"))


#~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*#
#------------------------------------------------------------------------------#
                                #BEGIN UI#
#The UI section deals with all the visuals that the user will interact with
#This section is split up into 3 tabs, Input, Output, and About
#The Input section has 2 tabs, parameters and ions, these have default values
#So that the code can be ran without inputting anything, but the user can upload
#an excel file that will overwrite this data, and the user can edit it further
#within these tabs
#The Output section displays the plots where the x and y-axis can be converted
#into other units
#The about section contains information about the tool for further details
#The side bar panel which is present throughout the app is used for minor 
#adjustments as well as the axis conversions and file inputs
#------------------------------------------------------------------------------#


ui <- fluidPage(theme=shinytheme("united"),
  
useShinyjs(),
                
  

#------------------------------------------------------------------------------#
                                 #SIDEBAR PANEL  
#------------------------------------------------------------------------------#


  sidebarLayout(
    sidebarPanel(
      fileInput("file1", "Choose .xlsx File", accept = ".xlsx"),
      textOutput("reject"),
      textOutput("OutputConcentration"),
      selectInput("OCunits", "Output Concentration Units", c("mg/L", "ug/L", "ng/L", "c/c0")),
      selectInput("timeunits","Output Time Units",c("hr", "day", "month", "bed volumes")),
      sliderInput("nrv", "Radial Collocation Points",0, 20, 7),
      sliderInput("nzv", "Axial Collocation Points", 0, 20, 13),
      actionButton("run_button", "Run Analysis", icon=icon("play")),
      downloadButton("save_button", "Save Data"),
      textOutput("ionadded"),
      textOutput("concentrationadded"),
      textOutput("analysisran")),
    

#------------------------------------------------------------------------------#
                              #MAIN PANEL
#------------------------------------------------------------------------------#
                    
   mainPanel(
    navbarPage("Ion Exchange Model",


#------------------------------------------------------------------------------#
                              #INPUT TAB
#------------------------------------------------------------------------------#

      tabPanel("Input",
       tabsetPanel(
         

         
        tabPanel("Parameters",
                                                   
         fluidRow(
            column(1,
              br(), br(), br(),
              br(), br(), 
              textOutput("RC")),
            column(2, offset=1,
              br(),
              textOutput("Q"),
              br(), br(), 
              textOutput("rb"),
              br(), br(), 
              textOutput("EBED")),
            column(3,
              numericInput("Qv", "", 1400, min=0),
              numericInput("rbv", "", 0.03375, min=0),
              numericInput("EBEDv", "", 0.35, min=0)),
            column(3,
              selectInput("qunits", "", c("meq/L")),
              selectInput("rbunits", "", c("cm", "m", "mm", "in", "ft")),
              selectInput("EBEDunits", "", c("")))),
              hr(),
     
                                                   #Parameters Row 2#
                                                   
        fluidRow(
            column(1,
              br(), br(), br(),
              br(), br(), br(),
              textOutput("CS")),
            column(2, offset=1,
              br(), 
              textOutput("Length"),
              br(), br(), br(),
              textOutput("Velocity"),
              br(), br(), br(),
              textOutput("Diameter"),
              br(), br(),
              textOutput("Flowrate")),
            column(3,
              numericInput("Lv", "",14.7646875, min=0),
              numericInput("Vv", "", 0.122857846019418, min=0),
              numericInput("Dv", "", 4, min=0),
              numericInput("Fv", "",12, min=0)),
            column(3,
              selectInput("LengthUnits", "", c("cm", "m", "mm", "in", "ft")),
              selectInput("velocityunits", "", c("cm/s", "m/s", "m/min", "m/h", "in/s","ft/s","ft/min", "gpm/ft^2")),
              selectInput("DiameterUnits","",c("cm")),
              selectInput("flowrateunits","",c("cm^3/s", "m^3/s", "ft^3/s", "mL/s", "L/min", "mL/min", "gpm", "mgd"))),
            column(2, br(), br(), br(),
              radioButtons("veloselect", "", c("Linear", "Volumetric")))),
              hr(),
                                                  
                                                   #Parameters Row 4#
        
        fluidRow(
          column(1,
            br(), br(),
            textOutput("conctime")),
            column(3, offset=6,
            selectInput("timeunits2", "", c("hr", "day")))),
        
        h4("Ion List"),
        dataEditUI("edit-1"),
        br(), br(),
        h4("Concentration Points"),
        dataEditUI("edit-2"),
        br(), br()
        
        ))),
        
        


#------------------------------------------------------------------------------#
                              #ION TAB
#------------------------------------------------------------------------------#
        
                               

      # tabPanel("Ions",
      #   h4("Ion List"),
      #   dataEditUI("edit-1"),
      #   br(), br(),
      #   h4("Concentration Points"),
      #   dataEditUI("edit-2"),
      #   br(), br()
      #   ))),
      

#------------------------------------------------------------------------------#
                                #OUTPUT TAB
#------------------------------------------------------------------------------#
      

    tabPanel("Output",
        shinycssloaders::withSpinner(
        plotlyOutput("Plot")),
        br(),
        plotlyOutput("ExtraChemicals")),


#------------------------------------------------------------------------------#
                                    #ABOUT TAB
#------------------------------------------------------------------------------#


    tabPanel("About",
        tableOutput("sum2"),
        textOutput("about"),
        br(),
        textOutput("how2use"))
                               
                    ))))



#~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*#









#~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*#
#------------------------------------------------------------------------------#
                                #BEGIN SERVER#
#------------------------------------------------------------------------------#



server <- function(input, output, session) {
  
#------------------------------------------------------------------------------#
                            #STATIC DISPLAY TEXTS#
#------------------------------------------------------------------------------#
  
  output$Q<-renderText("Resin Capacity")
  output$rb<-renderText("Bead Radius")
  output$EBED<-renderText("Bed Porosity")
  output$name<-renderText("Name")
  
  output$CS<-renderText("Column Specifications")
  output$MC<-renderText("Material Characteristics")
  output$CS3<-renderText("Solver Related")
  
  output$Length<-renderText("Length")
  output$Velocity<-renderText("Velocity")
  output$Diameter<-renderText("Diameter")
  output$Flowrate<-renderText("Flow Rate")
  
  output$kL<-renderText("Film Transfer Coefficient")
  output$Ds<-renderText("Surface Diffusion Coefficient")
  
  output$RC<-renderText("Resin Characteristics")
  
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
  
  output$about<-renderText("The Ion Exchange Model is a tool used to predict the concentration of PFAS chemicals over time
                           as a function of the water treatment apparatus. The computational model was developed by Levi Halpert,
                          and _____ at the Environmental Protection Agency. To read more about computations used in this tool,
                           one can read more here _______")
  
  output$how2use<-renderText("There are two ways to start this model. 1) Is to use an excel file to describe parameters of
                             water treatment apparatus which must follow this format ______. One can upload such file
                             by clicking 'upload xlsx' in the top left corner. 2) Is to start with the data that is
                             provided in the user interface and manipulate the data from there.
                             Once the parameters have been decided ions can be added, either in the xlsx file or on the ions tab,
                             as well as concentration points. When the user is satisfied with their settings, click 'run analysis'
                             to begin the computation. This may take a while, espeically with the more ions that have been added.

                              ")
  
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

  #When a file is uploaded, the session is reloaded. We do this because there does
  #Not seem to be any other way to overwrite the DataEditR tables. See the process_file
  #Notes for further elaboration.
  observeEvent(input$file1, {
    session$reload()
  })
  

  
#------------------------------------------------------------------------------#
                        #PARAMS DATA HANDLING#
#------------------------------------------------------------------------------#  
  
  #When the param sheet is read in, make it a reactiveVal so that the data can
  #Be saved between refreshes and edited
  paramsheet<-reactiveVal(read.csv("paramsheet.csv"))
  
  #Create a default dataframe that gets defined so an error isn't thrown and
  #Maybe there is a chance someone does not want to use an xlsx file
  paramdataframe<-reactiveVal()
  paramvals<-reactiveValues()
  
  #Creating the default values for the data frame
  observe({
    paramvals$Qv<-input$Qv
    paramvals$EBEDv<-input$EBEDv
    paramvals$Lv<-input$Lv
    paramvals$Vv<-input$Vv
    paramvals$rbv<-input$rbv
    paramvals$kLv<-input$kLv
    paramvals$Dsv<-input$Dsv
    paramvals$nrv<-input$nrv
    paramvals$nzv<-input$nzv
    paramvals$time<-1})
  
  #This Data frame is set up by default of all the default parameter values
  observe({paramdataframe(data.frame(
    name=c("Q", "EBED", "L", "v", "rb", "kL", "Ds", "nr", "nz", "time"),
    value=c(paramvals$Qv, paramvals$EBEDv, paramvals$Lv, paramvals$Vv, paramvals$rbv, NA, NA, paramvals$nrv, paramvals$nzv, NA),
    units=c(input$qunits, NA, input$LengthUnits, input$velocityunits, input$rbunits, NA, NA, NA, NA, input$timeunits)
  ))})

  #Take the data from the file that the user uploaded and overwrite the default frame
  capacity<-reactive({filter(paramsheet(), name=="Q")$value})
  observe({updateNumericInput(session, "Qv", value=capacity())})
  
  eebed<-reactive({filter(paramsheet(), name=="EBED")$value})
  observe({updateNumericInput(session, "EBEDv", value=eebed())})
  
  #Length is a protected variable
  length2<-reactive({filter(paramsheet(), name=="L")$value})
  observe({updateNumericInput(session, "Lv", value=length2())})
  
  velocity<-reactive({filter(paramsheet(), name=="v")$value})
  observe({updateNumericInput(session, "Vv", value=velocity())})
  
  beadradius<-reactive({filter(paramsheet(), name=="rb")$value})
  observe({updateNumericInput(session, "rbv", value=beadradius())})
  
  film<-reactive({filter(paramsheet(), name=="kL")$value})
  observe({updateNumericInput(session, "kLv", value=film())})
  
  diffuse<-reactive({filter(paramsheet(), name=="Ds")$value})
  observe({updateNumericInput(session, "Dsv", value=diffuse())})
  
  radial<-reactive({filter(paramsheet(), name=="nr")$value})
  observe({updateNumericInput(session, "nrv", value=radial())})
  
  axial<-reactive({filter(paramsheet(), name=="nz")$value})
  observe({updateNumericInput(session, "nzv", value=axial())})
  
  time<-reactive({1})
  

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
  lengthvector<-c("cm", "m", "mm", "in", "ft")
  velocityvector<-c("cm/s", "m/s", "m/min", "m/h", "in/s","ft/s","ft/min", "gpm/ft^2")
  
  lengthvector2<-eventReactive(input$file1, {c(paramdataframe()$units[3], lengthvector)})
  lengthvector3<-eventReactive(input$file1, {unique(lengthvector2())})
  
  velocityvector2<-eventReactive(input$file1, {c(paramdataframe()$units[4], velocityvector)})
  velocityvector3<-eventReactive(input$file1, {unique(velocityvector2())})
  
  rbvector<-eventReactive(input$file1, {c(paramdataframe()$units[5], lengthvector)})
  rbvector2<-eventReactive(input$file1, {unique(rbvector())})
  
  observe({updateSelectInput(session, "rbunits", choices=rbvector2())})
  observe({updateSelectInput(session, "LengthUnits", choices=lengthvector3())})
  observe({updateSelectInput(session, "velocityunits", choices=velocityvector3())})
  #observe({updateSelectInput(session, "flowrateunits", value=filter(paramdat(), name=="Q")$units)})
  
  observe({
    toggleState("Vv", condition=input$veloselect!="Volumetric")
    toggleState("Fv", condition=input$veloselect!="Linear")
    toggleState("Dv", condition=input$veloselect!="Linear")
  })
  
  velocityvar<-reactiveVal()
  
  observe({
    if(input$veloselect=="Linear"){
      velocityvar(input$Vv)
    }
    if(input$veloselect=="Volumetric"){
      updateNumericInput(session, "Vv", value=input$Fv/(pi*((input$Dv/2)**2)))
      
    }
  })
  
#------------------------------------------------------------------------------#
                           #IONS DATA HANDLING#
#------------------------------------------------------------------------------#  
  

  iondat<- dataEditServer("edit-1", data = "ionsheet.csv")
  dataOutputServer("output-1", data = iondat)
  

  cindat<-dataEditServer("edit-2", data="cinsheet.csv")
  dataOutputServer("output-2", data=cindat)
  

#------------------------------------------------------------------------------#
                          #INPUT UNIT CONVERSIONS#
#------------------------------------------------------------------------------#    
  
  
  observeEvent(input$rbunits,{
    if(input$rbunits=="m"){
      paramvals$rbv<-input$rbv*m2cm
    }
    else if(input$rbunits=="mm"){
      paramvals$rbv<-input$rbv*mm2cm
    }
    else if(input$rbunits=="cm"){
      paramvals$rbv<-input$rbv
    }
    else if(input$rbunits=="in"){
      paramvals$rbv<-input$rbv*in2cm
    }
    else if(input$rbunits=="ft"){
      paramvals$rbv<-input$rbv*12 * in2cm
    }
  })
  
 
  observeEvent(input$LengthUnits, {
    if(input$LengthUnits=="m"){
      paramvals$Lv<-input$Lv*m2cm
    }
    else if(input$LengthUnits=="mm"){
      paramvals$Lv<-input$Lv*mm2cm
    }
    else if(input$LengthUnits=="cm"){
      paramvals$Lv<-input$Lv
    }
    else if(input$LengthUnits=="in"){
      paramvals$Lv<-input$Lv*in2cm
    }
    else if(input$LengthUnits=="ft"){
      paramvals$Lv<-input$Lv*12 * in2cm
    }
  })
  
  observeEvent(input$velocityunits,{
    if(input$velocityunits=="ft/s"){
      paramvals$Vv<- input$Vv*12 * in2cm
    }
    else if(input$velocityunits=="m/s"){
      paramvals$Vv<- input$Vv*m2cm
    }
    else if(input$velocityunits=="cm/s"){
      paramvals$Vv<- input$Vv
    }
    else if(input$velocityunits=="in/s"){
      paramvals$Vv<- input$Vv*in2cm
    }
    else if(input$velocityunits=="m/min"){
      paramvals$Vv<- input$Vv*mpmin2cmps
    }
    else if(input$velocityunits=="ft/min"){
      paramvals$Vv<- input$Vv*ftpmin2cmps
    }
    if(input$velocityunits=="m/h"){
      paramvals$Vv<-input$Vv*mph2cmps
    }
    else if(input$velocityunits=="gpm/ft^2"){
      paramvals$Vv<-input$Vv*gpmpft2cmps
    }
  })
  
  observeEvent(input$filmunits, {
    if(input$filmunits=="ft/s"){
      paramvals$kLv<- input$kLv*12 * in2cm
    }
    else if(input$filmunits=="m/s"){
      paramvals$kLv<- input$kLv*m2cm
    }
    else if(input$filmunits=="cm/s"){
      paramvals$kLv<- input$kLv
    }
    else if(input$filmunits=="in/s"){
      paramvals$kLv<- input$kLv*in2cm
    }
    else if(input$filmunits=="m/min"){
      paramvals$kLv<- input$kLv*mpmin2cmps
    }
    else if(input$filmunits=="ft/min"){
      paramvals$kLv<- input$kLv*ftpmin2cmps
    }
  })
  
  observeEvent(input$diffusionunits,{
    if(input$diffusionunits=="ft/s^2"){
      paramvals$Dsv<-input$Dsv*ftps22cmps2
    }
    else if(input$diffusionunits=="m/s^2"){
      paramvals$Dsv<-input$Dsv*mps22cmps2
    }
    else if(input$diffusionunits=="cm/s^2"){
      paramvals$Dsv<-input$Dsv
    }
    else if(input$diffusionunits=="in/s^2"){
      paramvals$Dsv<-input$Dsv*inps22cmps2
    }
    else if(input$diffusionunits=="m/min^2"){
      paramvals$Dsv<-input$Dsv*mph2cmps
    }
    else if(input$diffusionunits=="ft/min^2"){
      paramvals$Dsv<-input$Dsv*ftpm22cmps2
    }
  })
  
  
  
  timeconverter<-reactiveVal()
  
  observeEvent(input$timeunits2, {
    if(input$timeunits2=="hr"){
      timeconverter(3600)
    }
    if(input$timeunits2=="day"){
      timeconverter(86400)
    }
  })
  
  
  
  output$params2<-renderTable(conc_convert_list)
  output$ions2<-renderTable(iondat$dat)
  output$cin2<-renderTable(cindat())
  
  
  
  # conc_convert_list<-list()
  # conc_vector<-reactiveVal()
  # correct_conc_units<-reactiveVal()
  # 
  # observe({
  #   for(unit in 1:nrow(iondat())){
  #     if(iondat()[unit,'conc_units']=="mg"){
  #       conc_convert_list[unit]<-1
  #     }
  #     else if(iondat()[unit,'conc_units']=="ug"){
  #       conc_convert_list[unit]<-1000
  #     }
  #     else if(iondat()[unit,'conc_units']=="ng"){
  #       conc_convert_list[unit]<-1000000
  #     }
  #     conc_vector(conc_convert_list)
  #   }
  # })
  # 
  # output$sum2<-renderTable(conc_vector())
  # 


  
#------------------------------------------------------------------------------#
                        #CALLING THE HSDMIX FUNCTION#
#------------------------------------------------------------------------------#    
  
  out<-reactiveVal()
  
  observeEvent(input$run_button, {
    out(HSDMIX_solve(paramdataframe(), iondat(), cindat(), timeconverter(), nt_report))
  })

  # find outlet indices
  outlet_id <- reactive({dim(out()[[2]])[4]})
  liquid_id <- reactive({dim(out()[[2]])[2]})
  
  
  
  
  
#------------------------------------------------------------------------------#
  #IEX CONCENTRATION OUTPUT DATAFRAME#
  #------------------------------------------------------------------------------#
  
  timeframe<-reactive({data.frame(hours=out()[[1]])})
  allchemicalconcs<-list()
  
  
  allchemicals<-eventReactive(input$run_button, {for (x in 1:nrow(iondat())){
    conc<-out()[[2]][, liquid_id(), x, outlet_id()]
    allchemicalconcs[[x]]<-conc
  }
    allconcdf<-data.frame(allchemicalconcs)
    colnames(allconcdf)<-iondat()$name
    allconcdf
  })
  
  
  massvector<-reactive({c(iondat()$mw/iondat()$valence)})
  allchemicalscorrected<-reactive({mapply('*', allchemicals(), massvector())})
  allchemicalscorrected2<-reactive({data.frame(allchemicalscorrected())})
  allchemicalscorrected3<-reactive({tidyr::gather(allchemicalscorrected2())})
  allchemicalscorrected4<-reactive({data.frame(name=allchemicalscorrected3()[,1],
                                               conc=allchemicalscorrected3()[,2])})
  allchems<-reactive({cbind(timeframe(), allchemicalscorrected4())})
  
  allchemicals2<-reactive({cbind(timeframe(), allchemicals())})
  
  
  
  
  
  #------------------------------------------------------------------------------#
  #END IEX CONCENTRATION OUTPUTDATAFRAME#
  #------------------------------------------------------------------------------#
  
  
  #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*#
  
  
  #------------------------------------------------------------------------------#
  #GENERATE C/C0 DATAFRAMES FROM MG/L CONCENTRATION DATAFRAMES#
  
  #The goal here is: For each column in allchemicals, divide every element
  #in each column by the first value in that column.
  #------------------------------------------------------------------------------#
  
  
  cc0vector<-reactive({c(cindat()[2,2:ncol(cindat())])})
  allchemicalscc0<-reactive({mapply('/', allchemicals(), cc0vector())})
  allchemicalscc02<-reactive(data.frame(allchemicalscc0()))
  bedvolume<-reactive({unlist(paramdataframe()$value[4])/unlist(paramdataframe()$value[5])})
  
  allchemicalscc03<-reactive({tidyr::gather(allchemicalscc02())})
  allchemicalscc04<-reactive({data.frame(name=allchemicalscc03()[,1],
                                         conc=allchemicalscc03()[,2])})
  
  
  
  #------------------------------------------------------------------------------#
  #END IEX CONCENTRATION OUTPUTDATAFRAME#
  #------------------------------------------------------------------------------#
  
  
  
  #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*#
  
  
  
  #------------------------------------------------------------------------------#
  #CONVERSION OF DATAFRAMES#
  #------------------------------------------------------------------------------#
  
  bonusdataframe<-data.frame(hours=c(), conc=c())
  
  allconcsconvert<-eventReactive(input$run_button, {for (x in 1:nrow(data_edit())){
    
    dx_frame<-data.frame(
      hours=out()[[1]], conc=out()[[2]][, liquid_id(), x, outlet_id()], name=data_edit()[x,1]
    )
    
    bonusdataframe<-rbind(bonusdataframe, dx_frame)
    
  }
    bonusdataframe
  })
  
  
  chemnames<-reactive({allconcsconvert()$name})
  allconcscc0<-reactive({rbind(allconcscc02(), chemnames())})
  
  
  #------------------------------------------------------------------------------#
  #END INITIALIZING CONVERSION FRAMES#
  #------------------------------------------------------------------------------#
  
  
  
  
  #~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*#
  
  
  #------------------------------------------------------------------------------#
  #CONVERTING OUTPUT DATAFRAMES#
  #------------------------------------------------------------------------------#
  
  
  allchemicals2<-reactive({cbind(timeframe(), allchemicals())})
  
  outputcounterions<-reactiveValues(counterion=0)
  outputions<-reactiveValues(ion=0)
  
  # 804 == 4 * nt_report
  ## counterIon_loc = 4 * nt_report
  ## addIon_loc = counterIon_loc + 1
  counteriondata<-reactive({allchems()[0:804,]})
  iondata<-reactive({allchems()[805:nrow(allchems()),]})
  
  counteriondatacc0<-reactive({ allchemicalscc04()[0:804,]})
  iondatacc0<-reactive({allchemicalscc04()[805:nrow(allchems()),]})
  
  
  outputcounterions$name<-reactive({counteriondata()$name})
  outputions$name<-reactive({iondata()$name})
  
  observe({
    req(allchemicals2())
    
    #bed
    if(input$timeunits=="hr"){
      outputcounterions$time<-counteriondata()$hours*1
      outputions$time<-iondata()$hours*1
    }
    if(input$timeunits=="day"){
      outputcounterions$time<-counteriondata()$hours/24
      outputions$time<-iondata()$hours/24
    }
    if(input$timeunits=="month"){
      outputcounterions$time<-counteriondata()$hours/720 #Assume 30 days in a month
      outputions$time<-iondata()$hours/720
    }
    if(input$timeunits=="bed volumes"){
      outputcounterions$time<-counteriondata()$hours/bedvolume()
      outputions$time<-iondata()$hours/bedvolume()
    }
  })
  
  observe({
    req(allchemicals2())
    
    if(input$OCunits=="c/c0"){
      outputcounterions$conc <-  counteriondatacc0()$conc
      outputions$conc<-iondatacc0()$conc
    }
    if(input$OCunits=="mg/L"){
      outputcounterions$conc <- counteriondata()$conc*1
      outputions$conc<-iondata()$conc*1
    }
    if(input$OCunits=="ug/L"){
      outputcounterions$conc <- counteriondata()$conc*1e6
      outputions$conc<-iondata()$conc*1e3
    }
    if(input$OCunits=="ng/L"){
      outputcounterions$conc <- counteriondata()$conc*1e6
      outputions$conc<-iondata()$conc*1e6
    }
  })
  
  
  
  processed_data <- reactive({
    
    plot_data <- counteriondata()
    plot_data$conc <- outputcounterions$conc
    plot_data$hours <- outputcounterions$time
    plot_data
    
  })
  
  processed_data2 <- reactive({
    
    plot_data2 <- iondata()
    plot_data2$conc <- outputions$conc
    plot_data2$hours <- outputions$time
    plot_data2
    
  })
  
  

  


  
  fig<-reactive({plot_ly(processed_data(), x=~hours, y=~conc,type='scatter', mode="lines", color=~name)})
  fig2<-reactive({fig()%>%layout(title="Concentration over Time",
                                 legend=list(orientation='h', y =1),
                                 xaxis=list(title=input$timeunits),
                                 yaxis=list(title=input$OCunits))})
  
  bonusfig<-reactive({plot_ly(processed_data2(), x=~hours, y=~conc,type='scatter', mode="lines", color=~name)})
  bonusfig2<-reactive({bonusfig()%>%layout(title="Concentration over Time", showlegend=TRUE,
                                           legend=list(orientation='h', y=1),
                                           xaxis=list(title=input$timeunits),
                                           yaxis=list(title=input$OCunits, showexponent='all', exponentformat='e'))})
  
  
  output$Plot<-renderPlotly(
    fig2())
  
  output$ExtraChemicals <- renderPlotly(
    bonusfig2())
  
  
  
  output$save_button<-downloadHandler(
    filename=function(){"IEX_Data.xlsx"},
    content=function(filename){
      write.xlsx(paramdataframe(), filename, sheetName="params")
      write.xlsx(iondat(), sheetName="ions", filename, append=TRUE)
      write.xlsx(cindat(), sheetName="cin", filename, append=TRUE)
      write.xlsx(allchemicals2(), sheetName="output", filename, append=TRUE)
    }
    
  )
  

  
}
#runGadget(ui, server, viewer = browserViewer(browser = getOption("browser")))

shinyApp(ui, server)