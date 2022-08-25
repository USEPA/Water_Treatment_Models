library(readxl)
library(shiny)
library(ggplot2)
library(dplyr)
library(xlsx)
library(shinythemes)
library(deSolve)
library(orthopolynom)
library(DT)

ui <- fluidPage(theme=shinytheme("united"),
                
                sidebarLayout(
                  sidebarPanel(
                    fileInput("file1", "Choose .xlsx File", accept = ".xlsx"),
                    actionButton("apply_button", "Apply Data"),
                    actionButton("run_button", "Run Analysis", icon=icon("play")),
                    tableOutput("sum"),
                    br(), br(),
                  
                  ),
                  
                  mainPanel(
                    
                    navbarPage("Ion Exchange Model",
                               
                               
                               tabPanel("Data",
                                        
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
                                                            br(), 
                                                            textOutput("rb"),
                                                            br(), br(), 
                                                            textOutput("EBED")),
                                                     column(3,
                                                           
                                                            numericInput("Qv", "", 1300),
                                                            
                                                            numericInput("rbv", "", 0.03375),
                                                          
                                                            numericInput("EBEDv", "", 0.35)),
                                                     
                                                     column(4,
                                                            selectInput("qunits", "", c("meq/L")),
                                                            selectInput("rbunits", "", c("")),
                                                            selectInput("EBEDunits", "", c("cm")))
                                                     ),

                                                    
                                                   
                                                   br(),
                                                   br(),
                                                   br(),
                                                   
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
                                                            numericInput("Lv", "",14.7646875),
                                                           
                                                            numericInput("Vv", "", 0.122857846019418),
                                                            
                                                            numericInput("Dv", "", 0),
                                                           
                                                            numericInput("Fv", "",0)),
                                                     column(4,
                                                            selectInput("LengthUnits", "", c("cm")),
                                                            selectInput("velocityunits", "", c("cm/s")),
                                                            selectInput("DiameterUnits","",c("cm^2")),
                                                            selectInput("flowrateunits","",c("cm2/s")))
                                                     
                                                     
                                                     ),
                                
                                                   
                                                   br(),
                                                   br(),
                                                   br(),
                                                   
                                                   #Parameters Row 3#
                                                   
                                                   fluidRow(
                                                     column(1,
                                                            br(), br(), br(),
                                                            br(), br(), br(),
                                                            textOutput("MC")),
                                                     column(2, offset=1,
                                                            br(), br(),
                                                            textOutput("kL"),
                                                            br(), br(), br(),
                                                            textOutput("Ds")),
                                                     column(3, 
                                                            br(),
                                                            numericInput("kLv", "",0.0021),
                                                            br(), br(),
                                                            numericInput("Dsv", "",0.0000002)),
                                                     column(4,
                                                            br(),
                                                            selectInput("filmunits","",c("cm/s")),
                                                            br(), br(),
                                                            selectInput("diffusionunits","",c("cm2/s")))
                                                     ),
                                                   
                    
                                                     
                                                   
                                                   br(),
                                                   br(),
                                                   br(),
                                                   
                                                   #Parameters Row 4#
                                                   
                                                   fluidRow(
                                                     column(1,
                                                            br(), br(),
                                                            textOutput("SR")),
                                                     column(2, offset=1,
                                                            br(),
                                                            textOutput("nr"),
                                                            br(), br(),
                                                            textOutput("nz")),
                                                     column(3,
                                                            numericInput("nrv", "",7),
                                                            br(), br(),
                                                            numericInput("nzv", "",13)),
                                                     column(4,
                                                     selectInput("radialunits", "", c("")),
                                                     br(), br(),
                                                     selectInput("axialunits","",c("")))),
                                                     
                                                   
                                                   br(), br(), br(),
                                                 
                                                   
                                                   #Parameters Row 5#
                                                   
                                                   fluidRow(
                                                    br(),
                                                     column(1,
                                                            textOutput("Time")),
                                                     column(2, offset=1,
                                                          
                                                            textOutput("TS")),
                                                     column(3, 
                                                            numericInput("tv", "",1)),
                                                    column(4,
                                                           selectInput("timeunits","",c("hr")))),
                                                   
                                                   textOutput("inputnum"),
                                            
                                                   br(), br(),
                                                  
                                                   br(), br(),
                                                   
                                                   tableOutput("paramsummary"),
                                                   tableOutput("ionsummary"),
                                                   tableOutput("cinsummary"),
                                                   tableOutput("summary"),
                                                   tableOutput("paramsummary2"),
                                                   tableOutput("ionsummary2"),
                                                   tableOutput("cinsummary2")
                                            
                                                   
                                                                   
                                          ),
                                          
                                         
                                          
                                          
                                          
                                          
                                          
                                          tabPanel("Ions",
                                                   DT::dataTableOutput("IonsTable"),
                                                   br(),
                                                   
                                                   fluidRow(
                                                     column(1, br(), br(), br(),
                                                            br(), br(),
                                                   actionButton("add", "Add Ion")),
                                                     column(2, offset=1,
                                                            textInput("name", "name"),
                                                            numericInput("mw", "mw", 1),
                                                            numericInput("initconc", "Inital Concentration", 5)),
                                                    column(3,
                                                           numericInput("KxA", "KxA", 1),
                                                           numericInput("valence", "valence", 1),
                                                           numericInput("finconc", "Final Concentration", 5)),
                                                           
                                                   column(4,
                                                          numericInput("kL", "kL", 1),
                                                          numericInput("Ds", "Ds", 1)),
                                                   tableOutput("dataview"),
                                                   tableOutput("summary2"))),
                                          
                                          tabPanel("Initial Concentration",
                                                   
                                                  DT::dataTableOutput("ICTable"),
                                                  tableOutput("summary3")
                                                   
                                                   
                                                   
                                                   
                                          ))),
                               
                               
                               
                               
                               tabPanel("Analysis",
                                        
                                        plotOutput("Plot"),
                                        plotOutput("ExtraChemicals")
                                        
                               ),
                               
                               tabPanel("Statistics",
                                        
                                        
                                        
                               ))
                    
                  )))


server <- function(input, output, session) {
  
  params2<-reactive({
    file <- input$file1
    ext <- tools::file_ext(file$datapath)
    
    req(file)
    validate(need(ext == "xlsx", "Please upload a csv file"))
    
    params<-read_xlsx(file$datapath, sheet=1)
    
    params
  })
  
  ion2<-reactive({
    file <- input$file1
    ext <- tools::file_ext(file$datapath)
    
    req(file)
    validate(need(ext == "xlsx", "Please upload a csv file"))
    
    ions<-read_xlsx(file$datapath, sheet=2)
    
    ions
  })
  
  output$IonsTable<-renderDataTable({ion2()})
  
  cin2<-reactive({
    file <- input$file1
    ext <- tools::file_ext(file$datapath)
    
    req(file)
    validate(need(ext == "xlsx", "Please upload a csv file"))
    
    cin2<-read_xlsx(file$datapath, sheet=3)
    
    cin2
  })
  
  output$ICTable<-renderDataTable({cin2()})
  
  capacity<-reactive({
    cap<-filter(params2(), name=="Q")$value
    cap})
  
  eebed<-reactive({
    val<-filter(params2(), name=="EBED")$value
    val})
  
  length2<-reactive({
    val<-filter(params2(), name=="L")$value
    val })
  
  velocity<-reactive({
    val<-filter(params2(), name=="v")$value
    val})
  
  beadradius<-reactive({
    val<-filter(params2(), name=="rb")$value
    val})
  
  
  film<-reactive({
    val<-filter(params2(), name=="kL")$value
    val})
  
  diffuse<-reactive({
    val<-filter(params2(), name=="Ds")$value
    val})
  
  
  radial<-reactive({
    val<-filter(params2(), name=="nr")$value
    val})
  
  axial<-reactive({
    val<-filter(params2(), name=="nz")$value
    val})
  
  time<-reactive({
    val<-filter(params2(), name=="time")$value
    val})
  
  observe({updateTextInput(session, "Vv", value=velocity())})
  observe({updateTextInput(session, "rbv", value=beadradius())})
  observe({updateTextInput(session, "Qv", value=capacity())})
  observe({updateTextInput(session, "EBEDv", value=eebed())})
  observe({updateTextInput(session, "Lv", value=length2())})
  observe({updateTextInput(session, "kLv", value=film())})
  observe({updateTextInput(session, "Dsv", value=diffuse())})
  observe({updateTextInput(session, "nrv", value=radial())})
  observe({updateTextInput(session, "nzv", value=axial())})
  observe({updateTextInput(session, "tv", value=time())  })
  
  
  
  
  newdataframe<-eventReactive(input$apply_button, {
    saveddataframe<-data.frame(
      name=c("Q", "EBED", "L", "v", "rb", "kL", "Ds", "nr", "nz", "time"),
      value=c(input$Qv, input$EBEDv, input$Lv, input$Vv, input$rbv, input$kLv, input$Dsv, input$nrv, input$nzv, input$tv),
      units=c(input$qunits, input$EBEDunits, input$rbunits, input$LengthUnits, input$velocityunits, input$filmunits, input$diffusionunits, input$radialunits, input$axialunits, input$timeunits)
      )
    saveddataframe
  })
  
  ionrow<-reactive({
    ionframe<-data.frame(
      name=c(input$name), mw=c(input$mw), KxA=c(input$KxA), valence=c(input$valence), kL=c(input$kL), Ds=c(input$Ds)
    )
    ionframe
  })
  
  
  iondataframe<-eventReactive(input$add, {
    newiondataframe<-rbind(ion2(), ionrow())
    newiondataframe
  })
  
  cincolumn<-reactive({
    cinframe<-data.frame(
      PFAS=c(input$initconc, input$finconc)
      )
    cinframe
  })
  
  cindataframe<-eventReactive(input$add, {
    newcindataframe<-cbind(cin2(), cincolumn())
  })
  
  #output$IonsTable<-eventReactive(input$add, {renderDataTable(iondataframe())})
  #
  
 
  output$dataview<-renderTable(ionrow())
  output$summary2<-renderTable(iondataframe())
  output$summary3<-renderTable(cindataframe())
 
  
  
  #------------------------------#
  #STATIC TEXT DISPLAYS#
  #------------------------------#
  
  
  output$Q<-renderText("Capacity of Chloride on Resin")
  output$rb<-renderText("Radius of Resin Bead")
  output$EBED<-renderText("Porosity of Bed")
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
  
  output$Time<-renderText("Time")
  output$TS<-renderText("Time Step")
  
  output$ChemicalNames<-renderText("Chemical Names")
  output$mw<-renderText("mw")
  output$KxA<-renderText("KxA")
  output$Valence<-renderText("Valence")
  
  output$Name2<-renderText("Name")
  output$InitialTime<-renderText("Inital")
  output$FinalTime<-renderText("Final")
  
  
  
  S_PER_HR <- 60 * 60 # seconds per hour
  
  
  # Inputs ----
  nt_report = 201 # number of reporting steps
  
  # Load input file ----
  
 

  
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
  
   out <- reactive({
     HSDMIX_solve(newdataframe(), iondataframe(), cindataframe(), nt_report)})


  # find outlet indices

    outlet_id <- reactive({dim(out()[[2]])[4]})
    liquid_id <- reactive({dim(out()[[2]])[2]})

 
    mytheme <-  theme(panel.background = element_rect(fill = "white", colour = NA),
                      panel.grid.major = element_line(colour = "grey70", size = 0.2),
                     panel.grid.minor = element_line(colour = "grey85", size = 0.5),
                     legend.title = element_text(colour = "black", size = 12, face = "bold", hjust = 0.5),
                      legend.text = element_text(colour = "black", size = 12),
                      legend.key.size = unit(1, "line"),
                      strip.text = element_text(colour = "black", size = 7),
                      axis.ticks = element_line(colour = "black", size = 1),
                      axis.line = element_line(colour = "black", size = 1, lineend = "square"),
                      axis.text.x = element_text(colour = "black", size = 8),
                      axis.text.y = element_text(colour = "black", size = 8),
                       axis.title.x = element_text(colour = "black", size = 15),
                      axis.title.y = element_text(colour = "black", size = 15),
                      plot.title=element_text(colour="black",size=15,face="bold", hjust=0.5))
 
 
  dat<-reactive({data.frame(hours = out()[[1]], conc = out()[[2]][, liquid_id(), 1, outlet_id()])})
  dat1<-reactive({data.frame(hours = out()[[1]], conc = out()[[2]][, liquid_id(), 2, outlet_id()])})
  dat2<-reactive({data.frame(hours = out()[[1]], conc = out()[[2]][, liquid_id(), 3, outlet_id()])})
  dat3<-reactive({data.frame(hours = out()[[1]], conc = out()[[2]][, liquid_id(), 4, outlet_id()])})
 
    bonusdataframe<-data.frame(hours=c(), conc=c())
    bonusdataframe2<-data.frame(hours=c(), conc=c())

   fulldata<-reactive({
     file <- input$file1
     ext <- tools::file_ext(file$datapath)

     req(file)
     validate(need(ext == "xlsx", "Please upload a xlsx file"))

     inputfile<-read_xlsx(file$datapath, sheet=2)
   })

   Index<-reactive({
    rownums<-nrow(ion2())
    rownums
   })


   bonusdataframe3<-eventReactive(input$apply_button, {for (x in 5:Index()){
 
     dx_frame<-data.frame(
       hours=out()[[1]], conc=out()[[2]][, liquid_id(), x, outlet_id()], chemical=ion2()[x,1]
     )

    bonusdataframe<-rbind(dx_frame, bonusdataframe2)

   }
    bonusdataframe
   })

   output$sum<-renderTable(bonusdataframe3())
  

   chlorideframe<-reactive({data.frame(
     hours=out()[[1]], conc=out()[[2]][, liquid_id(), 1, outlet_id()], Chemical=rep("Chloride", nrow(dat()))
  
   )
     })
  
   sulfateframe<-reactive({data.frame(
     hours=out()[[1]], conc=out()[[2]][, liquid_id(), 2, outlet_id()], Chemical=rep("Sulfate", nrow(dat()))
   )})
  
   bicarbonateframe<-reactive({data.frame(
     hours=out()[[1]], conc=out()[[2]][, liquid_id(), 3, outlet_id()], Chemical=rep("Bicarbonate", nrow(dat()))
   )})
  
   nitrateframe<-reactive({data.frame(
     hours=out()[[1]], conc=out()[[2]][, liquid_id(), 4, outlet_id()], Chemical=rep("Nitrate", nrow(dat()))
   )})
  
   alldata<-reactive({rbind(chlorideframe(),nitrateframe(),bicarbonateframe(),sulfateframe())})
  
   #output$dataview<-renderTable(chlorideframe())
  
  
# output$Plot<-renderPlot(
# ggplot(alldata(), mapping=aes(x=hours, y=conc, color=Chemical)) +
# geom_point() + mytheme + ggtitle("Counter Ion Concentration Over Time")
# )
#  
#  output$ExtraChemicals<-renderPlot(
#  ggplot(bonusdataframe3(), mapping=aes(x=hours, y=conc, color=name)) +
#  geom_point()  + mytheme+ ggtitle("Ion Concentration Over Time")
#   )
  
   observeEvent(input$run_button, {
     output$Plot<-renderPlot(

      ggplot(alldata(), mapping=aes(x=hours, y=conc, color=Chemical)) +
        geom_point() + mytheme + ggtitle("Counter-Ion Concentration over Time")
    )
 })

   observeEvent(input$run_button, {
    output$ExtraChemicals<-renderPlot(

      ggplot(bonusdataframe3(), mapping=aes(x=hours, y=conc, color=name)) +
        geom_point()  + mytheme+ ggtitle("Ion Concentration over Time")
    )
   })

}


shinyApp(ui, server)






