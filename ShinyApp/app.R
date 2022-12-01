library(readxl)
library(shiny)
library(ggplot2)
library(dplyr)
library(xlsx)
library(shinythemes)
library(deSolve)
library(orthopolynom)
library(DT)
library(tibble)

ui <- fluidPage(theme=shinytheme("united"),
                
                sidebarLayout(
                  sidebarPanel(
                    fileInput("file1", "Choose .xlsx File", accept = ".xlsx"),
                    textOutput("reject"),
                    #actionButton("apply_button", "Apply Inputs"),
                    textOutput("OutputConcentration"),
                    selectInput("OCunits", "Output Concentration", c("mg/L", "ug/L", "ng/L", "c/c0")),
                    numericInput("displacementtime", "Run Duration", 40),
                    selectInput("timeunits","",c("hr", "day", "month", "year")),
                    actionButton("run_button", "Run Analysis", icon=icon("play")),
                    textOutput("ionadded"),
                    textOutput("concentrationadded"),
                    textOutput("analysisran"),
                    br(), br(),
            
                    
                  ),
                  
                  mainPanel(
                    
                    navbarPage("Ion Exchange Model",
                               
                               
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
                                                            selectInput("rbunits", "", c("cm", "m", "mm", "in", "ft")),
                                                            selectInput("EBEDunits", "", c("")))
                                                   ),
                                                   
                                                   
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
                                                            numericInput("Lv", "",14.7646875),
                                                            
                                                            numericInput("Vv", "", 0.122857846019418),
                                                            
                                                            numericInput("Dv", "", 0),
                                                            
                                                            numericInput("Fv", "",0)),
                                                     column(4,
                                                            selectInput("LengthUnits", "", c("cm", "m", "mm", "in", "ft")),
                                                            selectInput("velocityunits", "", c("cm/s", "ft/s", "m/s", "in/s", "m/min", "ft/min")),
                                                            selectInput("DiameterUnits","",c("cm^2")),
                                                            selectInput("flowrateunits","",c("cm2/s")))
                                                     
                                                     
                                                   ),
                                                   
                                                   
                                                  hr(),
                                                   
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
                                                            selectInput("filmunits","",c("cm/s", "in/s", "m/min", "ft/min")),
                                                            br(), br(),
                                                            selectInput("diffusionunits","",c("cm2/s")))
                                                   ),
                                                   
                                                   
                                                   
                                                   
                                                  hr(),
                                                   
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
                                                   
                                                   
                                                  hr(),
                                                   
                                                   
                                                   #Parameters Row 5#
                                                   
                                                   
                                                   
                                                   br(), br(), br(),
                                                   
                                                   
                                                   
                                                   
                                                   
                                          ),
                                          
                                          
                                          
                                          
                                          
                                          
                                          
                                          tabPanel("Ions",
                                                   tableOutput("IonsTable"),
                                                   br(),
                                                   
                                                   fluidRow(
                                                     column(1, br(), br(), br(),
                                                            br(), br(),
                                                            actionButton("add", "Add Ion")),
                                                     column(2, offset=1,
                                                            textInput("name", "name"),
                                                            numericInput("mw", "mw", 1),
                                                            numericInput("avgconc", "Average Concentration", 5)),
                                                     column(3,
                                                            numericInput("KxA", "KxA", 1),
                                                            numericInput("valence", "valence", 1)),
                                                     
                                                     column(4,
                                                            numericInput("kL", "kL", 1),
                                                            numericInput("Ds", "Ds", 1),
                                                            textOutput("nameerror")),
                                                   )),
                                          
                                          tabPanel("Concentrations",
                                                   
                                                   DT::dataTableOutput("ICTable"),
                                                   
                                                   
                                                   
                                                   
                                                   
                                                   
                                          ))),
                               
                               
                               
                               
                               tabPanel("Output",
                                        
                                        plotOutput("Plot"),
                                        plotOutput("ExtraChemicals")
                                        
                               ),
                               
                               tabPanel("Statistics",
                               ),
                               
                               tabPanel("Summary",
                                        tableOutput("sum"),
                                        tableOutput("sum2"),
                                        tableOutput("sum3"),
                                        tableOutput("sum4"),
                                        textOutput("sum5"),
                                        #tableOutput("dataview"),
                                        #tableOutput("summary2")
                               )
                               
                    ))))


server <- function(input, output, session) {
  
  output$reject<-renderPrint({
    req(input$file1)
    
    if(input$file1$type != "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"){ stop("Please upload a .xlsx file")}
    
  })
  
  paramdataframe<-reactiveVal()
  paramvals<-reactiveValues()
  
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
  
  #This Dataframe is set up by default of all the default paramater values
    observe({paramdataframe(data.frame(
      name=c("Q", "EBED", "L", "v", "rb", "kL", "Ds", "nr", "nz", "time"),
      value=c(paramvals$Qv, paramvals$EBEDv, paramvals$Lv, paramvals$Vv, paramvals$rbv, paramvals$kLv, paramvals$Dsv, paramvals$nrv, paramvals$nzv, 1),
      units=c(input$qunits, input$EBEDunits, input$LengthUnits, input$velocityunits, input$rbunits, input$filmunits, input$diffusionunits, input$radialunits, input$axialunits, input$timeunits)
    ))})

  
  #This dataframe is created when a file is inputed, which is formatted exactly like paramdataframe  
  paramdat<-reactiveVal()
  
  observe({
    file <- input$file1
    ext <- tools::file_ext(file$datapath)
    
    req(file)
    validate(need(input$file1$type == "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet", "Please select xlsx"))
    
    params<-read_xlsx(file$datapath, sheet=1)
    
    paramdat(params)
  })
  
  
  #These values are taken from paramdat so that they can be used to update paramdataframe
  capacity<-eventReactive(input$file1, {
    val<-filter(paramdat(), name=="Q")$value
    val})

  eebed<-eventReactive(input$file1,{
    val<-filter(paramdat(), name=="EBED")$value
    val})

  length2<-eventReactive(input$file1,{
    val<-filter(paramdat(), name=="L")$value
    val })

  velocity<-eventReactive(input$file1, {
    val<-filter(paramdat(), name=="v")$value
    val})

  beadradius<-eventReactive(input$file1,{
    val<-filter(paramdat(), name=="rb")$value
    val})


  film<-eventReactive(input$file1,{
    val<-filter(paramdat(), name=="kL")$value
    val})

  diffuse<-eventReactive(input$file1,{
    val<-filter(paramdat(), name=="Ds")$value
    val})


  radial<-eventReactive(input$file1,{
    val<-filter(paramdat(), name=="nr")$value
    val})

  axial<-eventReactive(input$file1,{
    val<-filter(paramdat(), name=="nz")$value
    val})

  time<-eventReactive(input$file1,{
    val<-filter(paramdat(), name=="time")$value
    val})

  
  #once the values are defined, they update paramdataframe
  observe({updateNumericInput(session, "Vv", value=velocity())})
  observe({updateNumericInput(session, "rbv", value=beadradius())})
  observe({updateNumericInput(session, "Qv", value=capacity())})
  observe({updateNumericInput(session, "EBEDv", value=eebed())})
  observe({updateNumericInput(session, "Lv", value=length2())})
  observe({updateNumericInput(session, "kLv", value=film())})
  observe({updateNumericInput(session, "Dsv", value=diffuse())})
  observe({updateNumericInput(session, "nrv", value=radial())})
  observe({updateNumericInput(session, "nzv", value=axial())})
  observe({updateNumericInput(session, "tv", value=time())  })
  
  
  
  
  iondat<-reactiveVal()
  
  observe({
    file <- input$file1
    ext <- tools::file_ext(file$datapath)
    
    req(file)
    validate(need(input$file1$type == "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet", "Please select xlsx"))
    
    ions<-read_xlsx(file$datapath, sheet=2)
    
    iondat(ions)
  })
  
  cindat<-reactiveVal()
  
  observe({
    file <- input$file1
    ext <- tools::file_ext(file$datapath)
    
    req(file)
    validate(need(input$file1$type == "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet", "Please select xlsx"))
    
    cin2<-read_xlsx(file$datapath, sheet=3)
    cin2$time<-list(0, input$displacementtime)
    
    cindat(cin2)
  })
  
  #output$sum<-renderTable(paramdataframe())
  #output$sum2<-renderTable(iondat())
  #output$sum3<-renderTable(cindat())
  #output$sum4<-renderTable(paramdataframe())
  

  # capacity<-reactive({
  #   val<-filter(paramdat(), name=="Q")$value
  #   val})
  # 
  # eebed<-reactive({
  #   val<-filter(paramdat(), name=="EBED")$value
  #   val})
  # 
  # length2<-reactive({
  #   val<-filter(paramdat(), name=="L")$value
  #   val })
  # 
  # velocity<-reactive({
  #   val<-filter(paramdat(), name=="v")$value
  #   val})
  # 
  # beadradius<-reactive({
  #   val<-filter(paramdat(), name=="rb")$value
  #   val})
  # 
  # 
  # film<-reactive({
  #   val<-filter(paramdat(), name=="kL")$value
  #   val})
  # 
  # diffuse<-reactive({
  #   val<-filter(paramdat(), name=="Ds")$value
  #   val})
  # 
  # 
  # radial<-reactive({
  #   val<-filter(paramdat(), name=="nr")$value
  #   val})
  # 
  # axial<-reactive({
  #   val<-filter(paramdat(), name=="nz")$value
  #   val})
  # 
  # time<-reactive({
  #   val<-filter(paramdat(), name=="time")$value
  #   val})
  # 
  # observe({updateTextInput(session, "Vv", value=velocity())})
  # observe({updateTextInput(session, "rbv", value=beadradius())})
  # observe({updateTextInput(session, "Qv", value=capacity())})
  # observe({updateTextInput(session, "EBEDv", value=eebed())})
  # observe({updateTextInput(session, "Lv", value=length2())})
  # observe({updateTextInput(session, "kLv", value=film())})
  # observe({updateTextInput(session, "Dsv", value=diffuse())})
  # observe({updateTextInput(session, "nrv", value=radial())})
  # observe({updateTextInput(session, "nzv", value=axial())})
  # observe({updateTextInput(session, "tv", value=time())  })

  m2cm=100
  mm2cm=0.1
  cm2cm=1
  in2cm=2.54
  ft2cm=30.48
  mmin2cms=1.666667
  ftmin2cms=0.508
  mmin22cms2=0.027778
  ftmin22cms2=0.00846667
  sec2sec=1
  min2sec=60
  hour2sec=360
  day2sec=8640
  month2sec=259200
  year2sec=3153600

  observeEvent(input$rbunits, {
    if(input$rbunits=="m"){
      paramvals$rb<-paramvals$rb*100
    }
    if(input$rbunits=="mm"){
      paramvals$rb<-paramvals$rb*0.1
    }
    if(input$rbunits=="cm"){
      paramvals$rb<-paramvals$rb*1
    }
    if(input$rbunits=="in"){
      paramvals$rb<-paramvals$rb*2.54
    }
    if(input$rbunits=="ft"){
      paramvals$rb<-paramvals$rb*30.48
    }
  })

  observeEvent(input$LengthUnits, {
    if(input$LengthUnits=="m"){
      paramvals$Lv<-paramvals$Lv*100
    }
    if(input$LengthUnits=="mm"){
      paramvals$Lv<-paramvals$Lv*0.1
    }
    if(input$LengthUnits=="cm"){
      paramvals$Lv<-paramvals$Lv*1
    }
    if(input$LengthUnits=="in"){
      paramvals$Lv<-paramvals$Lv*2.54
    }
    if(input$LengthUnits=="ft"){
      paramvals$Lv<-paramvals$Lv*30.48
    }
  })

  observeEvent(input$velocityunits, {
    if(input$velocityunits=="ft/s"){
      paramvals$Vv<- paramvals$Vv*30.48
    }
    if(input$velocityunits=="m/s"){
      paramvals$Vv<- paramvals$Vv*100
    }
    if(input$velocityunits=="cm/s"){
      paramvals$Vv<- paramvals$Vv*1
    }
    if(input$velocityunits=="in/s"){
      paramvals$Vv<- paramvals$Vv*2.54
    }
    if(input$velocityunits=="m/min"){
      paramvals$Vv<- paramvals$Vv*1.666667
    }
    if(input$velocityunits=="ft/min"){
      paramvals$Vv<- paramvals$Vv*0.508
    }
  })

  observeEvent(input$filmunits, {
    if(input$filmunits=="ft/s"){
      paramvals$kLv<- paramvals$kLv*30.48
    }
    if(input$filmunits=="m/s"){
      paramvals$kLv<- paramvals$kLv*100
    }
    if(input$filmunits=="cm/s"){
      paramvals$kLv<- paramvals$kLv*1
    }
    if(input$filmunits=="in/s"){
      paramvals$kLv<- paramvals$kLv*2.54
    }
    if(input$filmunits=="m/min"){
      paramvals$kLv<- paramvals$kLv*1.666667
    }
    if(input$filmunits=="ft/min"){
      paramvals$kLv<- paramvals$kLv*0.508
    }
  })

  observeEvent(input$diffusionunits, {
    if(input$diffusionunits=="ft/s^2"){
      paramvals$Dsv<- paramvals$Dsv*0.328
    }
    if(input$diffusionunits=="m/s^2"){
      paramvals$Dsv<- paramvals$Dsv*0.01
    }
    if(input$diffusionunits=="cm/s^2"){
      paramvals$Dsv<- paramvals$Dsv*1
    }
    if(input$diffusionunits=="in/s^2"){
      paramvals$Dsv<- paramvals$Dsv*0.3937
    }
    if(input$diffusionunits=="m/min^2"){
      paramvals$Dsv<- paramvals$Dsv*36
    }
    if(input$diffusionunits=="ft/min^2"){
      paramvals$Dsv<- paramvals$Dsv*118.11
    }
  })

  # observeEvent(input$timeunits, {
  #   if(input$timeunits=="s"){
  #     timeconv<-input$timeunits*sec2sec
  #   }
  #   if(input$timeunits=="min"){
  #     timeconv<-input$timeunits*min2sec
  #   }
  #   if(input$timeunits=="hr"){
  #     timeconv<-input$timeunits*hour2sec
  #   }
  #   if(input$timeunits=="day"){
  #     timeconv<-input$timeunits*day2sec
  #   }
  #   if(input$timeunits=="month"){
  #     timeconv<-input$timeunits*month2sec
  #   }
  #   if(input$timeunits=="year"){
  #     timeconv<-input$timeunits*year2sec
  #   }
  # })

  timeconverter<-reactiveVal()

  observeEvent(input$timeunits, {
    if(input$timeunits=="hr"){
      timeconverter(3600)
    }
    if(input$timeunits=="day"){
      timeconverter(86400)
    }
    if(input$timeunits=="month"){
      timeconverter(2592000)
    }
  })

  # newdataframe<-eventReactive(input$apply_button, {
  #   saveddataframe<-data.frame(
  #     name=c("Q", "EBED", "L", "v", "rb", "kL", "Ds", "nr", "nz", "time"),
  #     value=c(input$Qv, input$EBEDv, lengthconv(), velocityconv(), beadradiusconv(), input$kLv, input$Dsv, input$nrv, input$nzv, 1),
  #     units=c(input$qunits, input$EBEDunits, input$LengthUnits, input$velocityunits, input$rbunits, input$filmunits, input$diffusionunits, input$radialunits, input$axialunits, input$timeunits)
  #   )
  #   saveddataframe
  # })

  observeEvent(input$add, {
    iondat(tibble::add_row(iondat(), name=input$name, mw=input$mw, KxA=input$KxA, valence=input$valence, kL=input$kL, Ds=input$Ds))
  })

  output$IonsTable<-renderTable({iondat()})



  observeEvent(input$add, {
    cindat(tibble::add_column(cindat(), !! input$name:=input$avgconc))
  })

  # cincolumn<-reactive({
  #   cinframe<-data.frame(
  #     PFAS=c(input$avgconc, input$avgconc)
  #   )
  #   cinframe
  # })
  #
  # cin2<-eventReactive(input$add, {
  #   newcindataframe<-cbind(cin2(), cincolumn())
  # })

  output$ICTable<-renderDataTable({cindat()})

  #output$IonsTable<-eventReactive(input$add, {renderDataTable(iondataframe())})
  #


  #output$dataview<-renderTable(ionrow())
  output$summary2<-renderTable(iondat())
  output$summary3<-renderTable(cindat())



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

  #output$OutputConcentration<-renderText("Graph Settings")
  output$OC<-renderText("Units")

  observeEvent(input$add, {
    output$ionadded<-renderText("Ion Added")
  })

  observeEvent(input$add, {
    output$concentrationadded<-renderText("Concentration Added")
  })

  observeEvent(input$run_button, {
    output$analysisran<-renderText("Analysis is Running")
  })




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

  out<-reactiveVal()
  
  observeEvent(input$run_button, {
    out(HSDMIX_solve(paramdataframe(), iondat(), cindat(), timeconverter(), nt_report))})


  #output$sum2<-renderTable(iondat())
  #output$sum3<-renderTable(cindat())

  # find outlet indices

  outlet_id <- reactive({dim(out()[[2]])[4]})
  liquid_id <- reactive({dim(out()[[2]])[2]})


  mytheme <-  reactive({theme(panel.background = element_rect(fill = "white", colour = NA),
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
                              plot.title=element_text(colour="black",size=15,face="bold", hjust=0.5))})


  dat<-reactive({data.frame(hours = out()[[1]], conc = out()[[2]][, liquid_id(), 1, outlet_id()])})
  dat1<-reactive({data.frame(hours = out()[[1]], conc = out()[[2]][, liquid_id(), 2, outlet_id()])})
  dat2<-reactive({data.frame(hours = out()[[1]], conc = out()[[2]][, liquid_id(), 3, outlet_id()])})
  dat3<-reactive({data.frame(hours = out()[[1]], conc = out()[[2]][, liquid_id(), 4, outlet_id()])})

  bonusdataframe<-data.frame(hours=c(), conc=c())
  bonusdataframe2<-data.frame(hours=c(), conc=c())

  bonusdataframe3<-eventReactive(input$run_button, {for (x in 5:nrow(iondat())){

    dx_frame<-data.frame(
      hours=out()[[1]], conc=out()[[2]][, liquid_id(), x, outlet_id()], chemical=iondat()[x,1]
    )

    bonusdataframe<-rbind(bonusdataframe, dx_frame)

  }
    bonusdataframe
  })


 

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

  #
  #
  #
  # # outputcc0<-reactiveValues(chloride=0, sulfate=0, bicarbonate=0, nitrate=0)
  # #
  # # sulfateconverted<- observeEvent(input$run_button, {
  # #   req(sulfateframe())
  # #
  # #   outputcc0$sulfate<-sulfateframe()$conc/sulfateframe()[[0]]
  # #
  # # })
  # #
  # # chlorideconverted<-observeEvent(input$run_button, {
  # #   req(chlorideframe())
  # #
  # #   outputcc0$chloride<-chlorideframe()$conc/chlorideframe()[[0]]
  # # })
  # #
  # # bicarbonateconverted<-observeEvent(input$run_button, {
  # #   req(bicarbonateframe())
  # #
  # #   outputcc0$bicarbonate<-bicarbonateframe()$conc/bicarbonateframe()[[0]]
  # # })
  # #
  # # nitrateconverted<-observeEvent(input$run_button, {
  # #   req(nitrateframe())
  # #
  # #   outputcc0$nitrate<-nitrateframe()$conc/nitrateframe()[[0]]
  # # })
  #
  #
  #
  #
  alldata<-reactive({rbind(chlorideframe(),nitrateframe(),bicarbonateframe(),sulfateframe())})
  #alldatacc0<-reactive({data.frame(hours=alldata()$hours, conc=1, alldata()$Chemical)})

  outputall<-reactiveValues(counterion=0)
  outputbonus<-reactiveValues(ion=0)
  #cc0chloride<-reactiveValues(conc=0)
  
  #cc0values<-reactiveValues()
  #cc0frame<-reactiveValues()
  #cc0frame<-data.frame(hours=c(), conc=c())
  

  cc0valueschloride<-reactive({unlist(cindat()[2,2])})
  cc0valuessulfate<-reactive({unlist(cindat()[2,3])})
  cc0valuesbicarbonate<-reactive({unlist(cindat()[2,4])})
  cc0valuesnitrate<-reactive({unlist(cindat()[2,5])})
  
  #cc0chloride<-reactive(chlorideframe())
  #cc0chloride<-data.frame(hours=c(), conc=c())
  

  cc0chloride2<-reactive({
    req(chlorideframe())
    
    cc0chloride<-chlorideframe()
    cc0chloride$hours<-out()[[1]]
    cc0chloride$conc<-chlorideframe()$conc/cc0valueschloride()
    cc0chloride
  })
  
  cc0sulfate2<-reactive({
    req(sulfateframe())
    
    cc0sulfate<-sulfateframe()
    cc0sulfate$hours<-out()[[1]]
    cc0sulfate$conc<-sulfateframe()$conc/cc0valuessulfate()
    cc0sulfate
  })
  
  cc0bicarbonate2<-reactive({
    req(bicarbonateframe())

    cc0bicarbonate<-bicarbonateframe()
    cc0bicarbonate$hours<-out()[[1]]
    cc0bicarbonate$conc<-bicarbonateframe()$conc/cc0valuesbicarbonate()
    cc0bicarbonate
  })
  
  cc0nitrate2<-reactive({
    req(nitrateframe())

    cc0nitrate<-nitrateframe()
    cc0nitrate$hours<-out()[[1]]
    cc0nitrate$conc<-nitrateframe()$conc/cc0valuesnitrate()
    cc0nitrate
  })

  alldatacc0<-reactive({rbind(cc0chloride2(),cc0nitrate2(),cc0bicarbonate2(),cc0sulfate2())})
  
  
  bonusdf<-reactiveValues()

  #  bonusdf2<-reactive({
  #   req(bonusdataframe3())
  # 
  #   bonusdf<-bonusdataframe3()
  #   bonusdf$conc<-bonusdataframe3()$conc[0:201]/unlist(cindat()[2,6])
  #   bonusdf
  # })
   
   # bonusdf3<-reactive({
   #   req(bonusdataframe3())
   #   
   #   for(n in 1:1){
   #     
   #     bonusdf$i<-bonusdataframe3()$conc[(i-1)*201:i*201]
   #   }
   #   bonusdf
   # })
  
  bonuscc0<-data.frame(hours=c(), conc=c())

  bonusdf2<-reactive({
    req(bonusdataframe3())
    bonuscc0<-data.frame(hours=c(), conc=c())


    for(n in 1:1){

      dxframe<-data.frame(hours=out()[[1]], conc=(bonusdataframe3()$conc[(((n-1)*200)+1):(n*201)])/unlist(cindat()[2,n+5:nrow(iondat())]))

      bonuscc0<-rbind(bonuscc0, dxframe)

     }
    bonuscc0
  })

  #output$sum<-renderTable(bonusdf3())
  output$sum2<-renderTable(bonusdf2())
  output$sum3<-renderTable(cindat()[2,4])


  observeEvent(input$run_button, {
    req(alldata())

    if(input$timeunits=="hr"){
      outputall$time<-alldata()$hours*1
    }
    if(input$timeunits=="day"){
      outputall$time<-alldata()$hours/24
    }
    if(input$timeunits=="month"){
      outputall$time<-alldata()$hours/720 #Assume 30 days in a month
    }
    if(input$timeunits=="year"){
      outputall$time<-alldata()$hours/8760
    }
  })

  observeEvent(input$run_button, {
    req(bonusdataframe3())
    if(input$timeunits=="hr"){
      outputall$time<-bonusdataframe3()$hours*1
    }
    if(input$timeunits=="day"){
      outputall$time<-bonusdataframe3()$hours/24
    }
    if(input$timeunits=="month"){
      outputall$time<-bonusdataframe3()$hours/720 #Assume 30 days in a month
    }
    if(input$timeunits=="year"){
      outputall$time<-bonusdataframe3()$hours/8760
    }
  })

  #
  observeEvent(input$run_button,{
    req(alldata())
    #req(cc0frame2())

    if(input$OCunits=="c/c0"){
      outputall$counterion <- alldatacc0()$conc
    }
    if(input$OCunits=="mg/L"){
      outputall$counterion <- alldata()$conc*1
    }
    if(input$OCunits=="ug/L"){
      outputall$counterion <- alldata()$conc*1000
    }
    if(input$OCunits=="ng/L"){
      outputall$counterion <- alldata()$conc*1000000
    }
  })

  observeEvent(input$run_button, {
    req(bonusdataframe3())

    if(input$OCunits=="c/c0"){
      outputbonus$ion<-bonusdf2()$conc
    }
    if(input$OCunits=="mg/L"){
      outputbonus$ion<-bonusdataframe3()$conc*1
    }
    if(input$OCunits=="ug/L"){
      outputbonus$ion<-bonusdataframe3()$conc*1000
    }
    if(input$OCunits=="ng/L"){
      outputbonus$ion<-bonusdataframe3()$conc*1000000
    }
  })

  # output$sum<-renderTable(cc0chloride())
  # output$sum2<-renderTable(chlorideframe())
  # output$sum3<-renderTable(cc0valueschloride())
  #output$sum4<-renderTable(cc0frame())
  #output$sum5<-renderText(cc0chloriderframe())

  processed_data <- eventReactive(input$run_button,{
    req(alldata())

    plot_data <- alldata()
    plot_data$conc <- outputall$counterion
    plot_data$hours <- outputall$time
    plot_data
  })




  processed_data2 <- eventReactive(input$run_button, {
    req(bonusdataframe3())

    plot_data2 <- bonusdataframe3()
    plot_data2$conc <- outputbonus$ion
    plot_data2$hours <- outputall$time
    plot_data2
  })


  output$Plot <- renderPlot(
    ggplot(processed_data(), mapping=aes(x=hours, y=conc, color=Chemical)) +
      geom_point() + mytheme() + xlab(input$timeunits) + ylab(input$OCunits) +xlim(0,input$displacementtime) + ggtitle("Counter-Ion Concentration over Time")
  )

  output$ExtraChemicals <- renderPlot(
    ggplot(processed_data2(), mapping=aes(x=hours, y=conc, color=name)) +
      geom_point()  + mytheme() + xlab(input$timeunits) + ylab(input$OCunits) + ggtitle("Ion Concentration over Time")
  )

}


shinyApp(ui, server)