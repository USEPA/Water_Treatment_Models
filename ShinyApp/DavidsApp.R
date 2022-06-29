library(readxl)
library(shiny)
library(ggplot2)
library(dplyr)
library(xlsx)
library(shinythemes)


ui <- fluidPage(theme=shinytheme("united"),
                
                sidebarLayout(
                  sidebarPanel(
                    fileInput("file1", "Choose CSV File", accept = ".xlsx"),
                    actionButton("apply_button", "Apply Data"),
                    actionButton("run_button", "Run Analysis", icon=icon("play")),
                  ),
                  
                  mainPanel(
                    
                    navbarPage("Ion Exchange Model",
                               
                               
                               tabPanel("Data",
                                        
                                        tabsetPanel(
                                          tabPanel("Parameters",
                                                   
                                                   fluidRow(
                                                     column(1,
                                                            br(), br(), br(),
                                                            textOutput("RC")),
                                                     column(2, offset=1,
                                                            textOutput("Q"),
                                                            br(),
                                                            textOutput("rb"),
                                                            br(),
                                                            textOutput("EBED")),
                                                     column(3,
                                                            br(),
                                                            textOutput("Qv"),
                                                            br(),
                                                            textOutput("rbv"),
                                                            br(), br(),
                                                            textOutput("EBEDv"),
                                                            textOutput("namev")),
                                                     column(4,
                                                            br(),
                                                            textOutput("Qunit"),
                                                            br(),
                                                            textOutput("rbunit"),
                                                            br(), br(),
                                                            textOutput("EBEDunit"))),
                                                   
                                                   br(),
                                                   br(),
                                                   br(),
                                                   
                                                   #Parameters Row 2#
                                                   
                                                   fluidRow(
                                                     column(1,
                                                            br(), br(), br(),
                                                            textOutput("CS")),
                                                     column(2, offset=1,
                                                            textOutput("Length"),
                                                            br(),
                                                            textOutput("Velocity"),
                                                            br(),
                                                            textOutput("Diameter"),
                                                            br(),
                                                            textOutput("Flowrate")),
                                                     column(3,
                                                            textOutput("Lv"),
                                                            br(),
                                                            textOutput("Vv"),
                                                            br(),
                                                            textOutput("Dv"),
                                                            br(),
                                                            textOutput("Fv")),
                                                     column(4,
                                                            textOutput("Lvunit"),
                                                            br(),
                                                            textOutput("Vvunit"),
                                                            br(),
                                                            textOutput("Dvunit"),
                                                            br(),
                                                            textOutput("Fvunit"))),
                                                   
                                                   br(),
                                                   br(),
                                                   br(),
                                                   
                                                   #Parameters Row 3#
                                                   
                                                   fluidRow(
                                                     column(1,
                                                            br(), br(), br(),
                                                            textOutput("MC")),
                                                     column(2, offset=1,
                                                            br(),
                                                            textOutput("kL"),
                                                            br(),
                                                            textOutput("Ds")),
                                                     column(3, 
                                                            br(),
                                                            textOutput("kLv"),
                                                            br(), br(),
                                                            textOutput("Dsv")),
                                                     column(4,
                                                            br(),
                                                            textOutput("kLunit"),
                                                            br(), br(),
                                                            textOutput("Dsunit"))),
                                                   
                                                   br(),
                                                   br(),
                                                   br(),
                                                   
                                                   #Parameters Row 4#
                                                   
                                                   fluidRow(
                                                     column(1,
                                                            br(), br(),
                                                            textOutput("SR")),
                                                     column(2, offset=1,
                                                            textOutput("nr"),
                                                            br(),
                                                            textOutput("nz")),
                                                     column(3,
                                                            textOutput("nrv"),
                                                            br(), br(),
                                                            textOutput("nzv")),
                                                     column(4,
                                                            textOutput("nrunit"),
                                                            br(), br(),
                                                            textOutput("nzunit"))),
                                                   
                                                   br(),
                                                   br(),
                                                   br(),
                                                   
                                                   #Parameters Row 5#
                                                   
                                                   fluidRow(
                                                     column(1,
                                                            textOutput("Time")),
                                                     column(2, offset=1,
                                                            textOutput("TS")),
                                                     column(3, 
                                                            textOutput("tv")),
                                                     column(4,
                                                            textOutput("tunit"))),
                                          ),
                                          
                                          
                                          
                                          
                                          
                                          tabPanel("Ions",
                                                   
                                                   fluidRow(
                                                     column(1,
                                                            textOutput("ChemicalNames")),
                                                     column(2, offset=1,
                                                            textOutput("mw")),
                                                     column(3,
                                                            textOutput("KxA")),
                                                     column(4,
                                                            textOutput("Valence")),
                                                     #column(5,
                                                     #       textOutput("kL")),
                                                     # column(6,
                                                     #        textOutput("Ds"))
                                                     
                                                   ),
                                                   
                                                   br(), br(),
                                                   
                                                   #Ion Row 2#
                                                   
                                                   fluidRow(
                                                     column(1,
                                                            textOutput("Chloride")),
                                                     column(2, offset=1,
                                                            textOutput("Chloridemw")),
                                                     column(3,
                                                            textOutput("ChlorideKxA")),
                                                     column(4,
                                                            textOutput("ChlorideValence"))
                                                   ),
                                                   
                                                   #Ion Row 3#
                                                   
                                                   fluidRow(
                                                     column(1,
                                                            textOutput("Sulfate")),
                                                     column(2, offset=1,
                                                            textOutput("Sulfatemw")),
                                                     column(3,
                                                            textOutput("SulfateKxA")),
                                                     column(4,
                                                            textOutput("SulfateValence"))
                                                   ),
                                                   
                                                   #Ion Row 4#
                                                   
                                                   fluidRow(
                                                     column(1,
                                                            textOutput("Bicarbonate")),
                                                     column(2, offset=1,
                                                            textOutput("Bicarbonatemw")),
                                                     column(3,
                                                            textOutput("BicarbonateKxA")),
                                                     column(4,
                                                            textOutput("BicarbonateValence"))
                                                   ),
                                                   
                                                   #Ion Row 5#
                                                   
                                                   fluidRow(
                                                     column(1,
                                                            textOutput("Nitrate")),
                                                     column(2, offset=1,
                                                            textOutput("Nitratemw")),
                                                     column(3,
                                                            textOutput("NitrateKxA")),
                                                     column(4,
                                                            textOutput("NitrateValence"))
                                                     
                                                     
                                                   )),
                                          
                                          tabPanel("Initial Concentration",
                                                   
                                                   fluidRow(
                                                     column(1,
                                                            textOutput("Name2")),
                                                     column(2, offset=1,
                                                            textOutput("InitialTime")),
                                                     column(3,
                                                            textOutput("FinalTime"))
                                                   ),
                                                   
                                                   textOutput("Chloride2"),
                                                   textOutput("Sulfate2"),
                                                   textOutput("Bicarbonate2"),
                                                   textOutput("Nitrate2"),
                                                   textOutput("PFOA"),
                                                   
                                                   #Initial Concentration Row 2#
                                                   
                                                   # fluidRow(
                                                   #     column(1,
                                                   #        textOutput("Chloride2")),
                                                   #     column(2, offset=1,
                                                   #        textOutput("Chlorideti")),
                                                   #     column(3,
                                                   #        textOutput("Chloridetf"))
                                                   #   ),
                                                   
                                                   #Initial Concentration Row 3#
                                                   
                                                   # fluidRow(
                                                   #      column(1,
                                                   #        textOutput("Sulfate2")),
                                                   #      column(2,
                                                   #        textOutput("Sulfateti")),
                                                   #      column(3,
                                                   #        textOutput("Sulfatetf"))
                                                   #    ),
                                                   
                                                   #Initial Concentration Row 4#
                                                   
                                                   # fluidRow(
                                                   #     column(1,
                                                   #       textOutput("Bicarbonate2")),
                                                   #     column(2,
                                                   #       textOutput("Bicarbonateti")),
                                                   #     column(3,
                                                   #       textOutput("Bicarbonatetf"))
                                                   #   ),
                                                   
                                                   #Inital Concentration Row 5#
                                                   
                                                   #     fluidRow(
                                                   #         column(1,
                                                   #           textOutput("Nitrate2")),
                                                   #         column(2,
                                                   #           textOutput("Nitrateti")),
                                                   #         column(3,
                                                   #           textOutput("Nitratetf"))
                                                   # ),
                                                   
                                                   # fluidRow(
                                                   #   column(1,
                                                   #          textOutput("PFOA2")),
                                                   #   column(2,
                                                   #          textOutput("PFOAi")),
                                                   #   column(3,
                                                   #          textOutput("PFOAf"))
                                                   # )
                                                   
                                                   
                                                   
                                                   
                                          ))),
                               
                               
                               
                               
                               tabPanel("Analysis",
                                        
                                        plotOutput("Plot"),
                                        plotOutput("ExtraChemicals")
                                        
                               ),
                               
                               tabPanel("Statistics",
                                        
                                        
                                        
                               ))
                    
                  )))


server <- function(input, output) {
  
  params2<-reactive({
    file <- input$file1
    ext <- tools::file_ext(file$datapath)
    
    req(file)
    validate(need(ext == "xlsx", "Please upload a csv file"))
    
    params<-read_xlsx(file$datapath, sheet=1)
  })
  
  # observe({
  #   assign(
  #     fulldata = "your_global_variable", 
  #     value = input$file1$datapath, 
  #     envir = .GlobalEnv
  #   )
  # })
  
  
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
  
  
  
  
  output$Qv <- renderText({
    file <- input$file1
    ext <- tools::file_ext(file$datapath)
    
    req(file)
    validate(need(ext == "xlsx", "Please upload a csv file"))
    
    params<-read_xlsx(file$datapath, sheet=1)
    Q<-filter(params, name== "Q")$value
    
  })
  
  output$Qunit <- renderText({
    file <- input$file1
    ext <- tools::file_ext(file$datapath)
    
    req(file)
    validate(need(ext == "xlsx", "Please upload a csv file"))
    
    params<-read_xlsx(file$datapath, sheet=1)
    Qunitv<-filter(params, name=="Q")$units
    
  })
  
  
  output$rbv <- renderText({
    file <- input$file1
    ext <- tools::file_ext(file$datapath)
    
    req(file)
    validate(need(ext == "xlsx", "Please upload a csv file"))
    
    rb<-filter(params, name=="rb")$value
  })
  
  output$rbunit <- renderText({
    file <- input$file1
    ext <- tools::file_ext(file$datapath)
    
    req(file)
    validate(need(ext == "xlsx", "Please upload a csv file"))
    
    params<-read_xlsx(file$datapath, sheet=1)
    rbunit<-filter(params, name=="rb")$units
    
  })
  
  
  
  output$EBEDv <- renderText({
    file <- input$file1
    ext <- tools::file_ext(file$datapath)
    
    req(file)
    validate(need(ext == "xlsx", "Please upload a csv file"))
    
    EBED<-filter(params, name=="EBED")$value
  })
  
  output$EBEDunit <- renderText({
    file <- input$file1
    ext <- tools::file_ext(file$datapath)
    
    req(file)
    validate(need(ext == "xlsx", "Please upload a csv file"))
    
    params<-read_xlsx(file$datapath, sheet=1)
    EBEDunit<-filter(params, name=="EBED")$units
    
  })
  
  output$Lv <- renderText({
    file <- input$file1
    ext <- tools::file_ext(file$datapath)
    
    req(file)
    validate(need(ext == "xlsx", "Please upload a csv file"))
    
    Length<-filter(params, name=="L")$value
  })
  
  output$Lvunit <- renderText({
    file <- input$file1
    ext <- tools::file_ext(file$datapath)
    
    req(file)
    validate(need(ext == "xlsx", "Please upload a csv file"))
    
    params<-read_xlsx(file$datapath, sheet=1)
    Lvunits<-filter(params, name=="L")$units
    
  })
  
  output$Vv <- renderText({
    file <- input$file1
    ext <- tools::file_ext(file$datapath)
    
    req(file)
    validate(need(ext == "xlsx", "Please upload a csv file"))
    
    Velocity<-filter(params, name=="v")$value
  })
  
  output$Vvunit <- renderText({
    file <- input$file1
    ext <- tools::file_ext(file$datapath)
    
    req(file)
    validate(need(ext == "xlsx", "Please upload a csv file"))
    
    params<-read_xlsx(file$datapath, sheet=1)
    vunits<-filter(params, name=="v")$units
    
  })
  
  output$Dv <- renderText({
    file <- input$file1
    ext <- tools::file_ext(file$datapath)
    
    req(file)
    validate(need(ext == "xlsx", "Please upload a csv file"))
    
    Diameter<-filter(params, name=="Dv")$value
  })
  
  output$Dvunit <- renderText({
    file <- input$file1
    ext <- tools::file_ext(file$datapath)
    
    req(file)
    validate(need(ext == "xlsx", "Please upload a csv file"))
    
    params<-read_xlsx(file$datapath, sheet=1)
    Dvunits<-filter(params, name=="diameter")$units
    
  })
  
  output$Fv <- renderText({
    file <- input$file1
    ext <- tools::file_ext(file$datapath)
    
    req(file)
    validate(need(ext == "xlsx", "Please upload a csv file"))
    
    Flowrate<-filter(params, name=="Fv")$value
  })
  
  output$Fvunit <- renderText({
    file <- input$file1
    ext <- tools::file_ext(file$datapath)
    
    req(file)
    validate(need(ext == "xlsx", "Please upload a csv file"))
    
    params<-read_xlsx(file$datapath, sheet=1)
    Fvunits<-filter(params, name=="Flowrate")$units
    
  })
  
  output$kLv <- renderText({
    file <- input$file1
    ext <- tools::file_ext(file$datapath)
    
    req(file)
    validate(need(ext == "xlsx", "Please upload a csv file"))
    
    kLvv<-filter(params, name=="kL")$value
  })
  
  output$kLunit <- renderText({
    file <- input$file1
    ext <- tools::file_ext(file$datapath)
    
    req(file)
    validate(need(ext == "xlsx", "Please upload a csv file"))
    
    params<-read_xlsx(file$datapath, sheet=1)
    kLunits<-filter(params, name=="kL")$units
  })
  
  output$Dsv <- renderText({
    file <- input$file1
    ext <- tools::file_ext(file$datapath)
    
    req(file)
    validate(need(ext == "xlsx", "Please upload a csv file"))
    
    Dsvv<-filter(params, name=="Ds")$value
  })
  
  output$Dsunit <- renderText({
    file <- input$file1
    ext <- tools::file_ext(file$datapath)
    
    req(file)
    validate(need(ext == "xlsx", "Please upload a csv file"))
    
    params<-read_xlsx(file$datapath, sheet=1)
    Dsunits<-filter(params, name=="Ds")$units
  })
  
  output$nrv <- renderText({
    file <- input$file1
    ext <- tools::file_ext(file$datapath)
    
    req(file)
    validate(need(ext == "xlsx", "Please upload a csv file"))
    
    nrvv<-filter(params, name=="nr")$value
  })
  
  output$nrunit <- renderText({
    file <- input$file1
    ext <- tools::file_ext(file$datapath)
    
    req(file)
    validate(need(ext == "xlsx", "Please upload a csv file"))
    
    params<-read_xlsx(file$datapath, sheet=1)
    nrunits<-filter(params, name=="nr")$units
  })
  
  output$nzv <- renderText({
    file <- input$file1
    ext <- tools::file_ext(file$datapath)
    
    req(file)
    validate(need(ext == "xlsx", "Please upload a csv file"))
    
    nzvv<-filter(params, name=="nz")$value
  })
  
  output$nzunit <- renderText({
    file <- input$file1
    ext <- tools::file_ext(file$datapath)
    
    req(file)
    validate(need(ext == "xlsx", "Please upload a csv file"))
    
    params<-read_xlsx(file$datapath, sheet=1)
    nrunits<-filter(params, name=="nz")$units
  })
  
  output$tv <- renderText({
    file <- input$file1
    ext <- tools::file_ext(file$datapath)
    
    req(file)
    validate(need(ext == "xlsx", "Please upload a csv file"))
    
    tvv<-filter(params, name=="time")$value
  })
  
  output$tunit <- renderText({
    file <- input$file1
    ext <- tools::file_ext(file$datapath)
    
    req(file)
    validate(need(ext == "xlsx", "Please upload a csv file"))
    
    params<-read_xlsx(file$datapath, sheet=1)
    tvunits<-filter(params, name=="time")$units
  })
  
  output$Chloride <- renderText({
    file <- input$file1
    ext <- tools::file_ext(file$datapath)
    
    req(file)
    validate(need(ext == "xlsx", "Please upload a csv file"))
    
    ions<-read_xlsx(file$datapath, sheet=2)
    chloridename<-filter(ions, name=="CHLORIDE")$name
  })
  
  output$Chloridemw <- renderText({
    file <- input$file1
    ext <- tools::file_ext(file$datapath)
    
    req(file)
    validate(need(ext == "xlsx", "Please upload a csv file"))
    
    ions<-read_xlsx(file$datapath, sheet=2)
    chloridename<-filter(ions, name=="CHLORIDE")$mw
  })
  
  output$ChlorideKxA <- renderText({
    file <- input$file1
    ext <- tools::file_ext(file$datapath)
    
    req(file)
    validate(need(ext == "xlsx", "Please upload a csv file"))
    
    ions<-read_xlsx(file$datapath, sheet=2)
    chloridename<-filter(ions, name=="CHLORIDE")$KxA
  })
  
  output$ChlorideValence <- renderText({
    file <- input$file1
    ext <- tools::file_ext(file$datapath)
    
    req(file)
    validate(need(ext == "xlsx", "Please upload a csv file"))
    
    ions<-read_xlsx(file$datapath, sheet=2)
    chloridename<-filter(ions, name=="CHLORIDE")$valence
  })
  
  output$Sulfate <- renderText({
    file <- input$file1
    ext <- tools::file_ext(file$datapath)
    
    req(file)
    validate(need(ext == "xlsx", "Please upload a csv file"))
    
    ions<-read_xlsx(file$datapath, sheet=2)
    sulfatename<-filter(ions, name=="SULFATE")$name
  })
  
  output$Sulfatemw <- renderText({
    file <- input$file1
    ext <- tools::file_ext(file$datapath)
    
    req(file)
    validate(need(ext == "xlsx", "Please upload a csv file"))
    
    ions<-read_xlsx(file$datapath, sheet=2)
    sulfatename<-filter(ions, name=="SULFATE")$mw
  })
  
  output$SulfateKxA <- renderText({
    file <- input$file1
    ext <- tools::file_ext(file$datapath)
    
    req(file)
    validate(need(ext == "xlsx", "Please upload a csv file"))
    
    ions<-read_xlsx(file$datapath, sheet=2)
    sulfatename<-filter(ions, name=="SULFATE")$KxA
  })
  
  output$SulfateValence <- renderText({
    file <- input$file1
    ext <- tools::file_ext(file$datapath)
    
    req(file)
    validate(need(ext == "xlsx", "Please upload a csv file"))
    
    ions<-read_xlsx(file$datapath, sheet=2)
    sulfatename<-filter(ions, name=="SULFATE")$valence
  })
  
  
  output$Bicarbonate <- renderText({
    file <- input$file1
    ext <- tools::file_ext(file$datapath)
    
    req(file)
    validate(need(ext == "xlsx", "Please upload a csv file"))
    
    ions<-read_xlsx(file$datapath, sheet=2)
    bicarbonatename<-filter(ions, name=="BICARBONATE")$name
  })
  
  output$Bicarbonatemw <- renderText({
    file <- input$file1
    ext <- tools::file_ext(file$datapath)
    
    req(file)
    validate(need(ext == "xlsx", "Please upload a csv file"))
    
    ions<-read_xlsx(file$datapath, sheet=2)
    bicarbonatename<-filter(ions, name=="BICARBONATE")$mw
  })
  
  output$BicarbonateKxA <- renderText({
    file <- input$file1
    ext <- tools::file_ext(file$datapath)
    
    req(file)
    validate(need(ext == "xlsx", "Please upload a csv file"))
    
    ions<-read_xlsx(file$datapath, sheet=2)
    bicarbonatename<-filter(ions, name=="BICARBONATE")$KxA
  })
  
  output$BicarbonateValence <- renderText({
    file <- input$file1
    ext <- tools::file_ext(file$datapath)
    
    req(file)
    validate(need(ext == "xlsx", "Please upload a csv file"))
    
    ions<-read_xlsx(file$datapath, sheet=2)
    bicarbonatename<-filter(ions, name=="BICARBONATE")$valence
  })
  
  
  output$Nitrate <- renderText({
    file <- input$file1
    ext <- tools::file_ext(file$datapath)
    
    req(file)
    validate(need(ext == "xlsx", "Please upload a csv file"))
    
    ions<-read_xlsx(file$datapath, sheet=2)
    nitratename<-filter(ions, name=="NITRATE")$name
  })
  
  output$Nitratemw <- renderText({
    file <- input$file1
    ext <- tools::file_ext(file$datapath)
    
    req(file)
    validate(need(ext == "xlsx", "Please upload a csv file"))
    
    ions<-read_xlsx(file$datapath, sheet=2)
    nitratename<-filter(ions, name=="NITRATE")$mw
  })
  
  output$NitrateKxA <- renderText({
    file <- input$file1
    ext <- tools::file_ext(file$datapath)
    
    req(file)
    validate(need(ext == "xlsx", "Please upload a csv file"))
    
    ions<-read_xlsx(file$datapath, sheet=2)
    nitratename<-filter(ions, name=="NITRATE")$KxA
  })
  
  output$NitrateValence <- renderText({
    file <- input$file1
    ext <- tools::file_ext(file$datapath)
    
    req(file)
    validate(need(ext == "xlsx", "Please upload a csv file"))
    
    ions<-read_xlsx(file$datapath, sheet=2)
    nitratename<-filter(ions, name=="NITRATE")$valence
  })
  
  data2<-eventReactive(input$file1, {
    file<-input$file1
    ext<-tools::file_ext(file$datapath)
    y<-read.xlsx(file$datapath, 3, header=FALSE)
    
    g<-noquote(y[,3])
    return(g)
  })
  
  data3<-eventReactive(input$file1, {
    file<-input$file1
    ext<-tools::file_ext(file$datapath)
    y<-read.xlsx(file$datapath, 3, header=FALSE)
    
    g<-noquote(y[,4])
    return(g)
  })
  
  data4<-eventReactive(input$file1, {
    file<-input$file1
    ext<-tools::file_ext(file$datapath)
    y<-read.xlsx(file$datapath, 3, header=FALSE)
    
    g<-noquote(y[,5])
    return(g)
  })
  
  data5<-eventReactive(input$file1, {
    file<-input$file1
    ext<-tools::file_ext(file$datapath)
    y<-read.xlsx(file$datapath, 3, header=FALSE)
    
    g<-noquote(y[,6])
    return(g)
  })
  
  data6<-eventReactive(input$file1, {
    file<-input$file1
    ext<-tools::file_ext(file$datapath)
    y<-read.xlsx(file$datapath, 3, header=FALSE)
    
    g<-noquote(y[,7])
    return(g)
  })
  
  output$Chloride2<-renderPrint({data2()})
  output$Sulfate2<-renderPrint({data3()})
  output$Bicarbonate2<-renderPrint({data4()})
  output$Nitrate2<-renderPrint({data5()})
  output$PFOA<-renderPrint({data6()})
  
  
  
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
  
  
  dat<-data.frame(hours = out[[1]], conc = out[[2]][, liquid_id, 1, outlet_id])
  dat1<-data.frame(hours = out[[1]], conc = out[[2]][, liquid_id, 2, outlet_id])
  dat2<-data.frame(hours = out[[1]], conc = out[[2]][, liquid_id, 3, outlet_id])
  dat3<-data.frame(hours = out[[1]], conc = out[[2]][, liquid_id, 4, outlet_id])
  
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
    rownums<-nrow(fulldata())
  })
  
  
  bonusdataframe3<-reactive({for (x in 5:Index()){
    
    dx_frame<-data.frame(
      hours=out[[1]], conc=out[[2]][, liquid_id, x, outlet_id], chemical=fulldata()[x,1]
    )
    
    bonusdataframe<-rbind(dx_frame, bonusdataframe2)
    
  }
    bonusdataframe
  })
  
  
  chlorideframe<-data.frame(
    hours=out[[1]], conc=out[[2]][, liquid_id, 1, outlet_id], Chemical=rep("Chloride", nrow(dat))
  )
  
  sulfateframe<-data.frame(
    hours=out[[1]], conc=out[[2]][, liquid_id, 2, outlet_id], Chemical=rep("Sulfate", nrow(dat))
  )
  
  bicarbonateframe<-data.frame(
    hours=out[[1]], conc=out[[2]][, liquid_id, 3, outlet_id], Chemical=rep("Bicarbonate", nrow(dat))
  )
  
  nitrateframe<-data.frame(
    hours=out[[1]], conc=out[[2]][, liquid_id, 4, outlet_id], Chemical=rep("Nitrate", nrow(dat))
  )
  
  alldata<-rbind(chlorideframe,nitrateframe,bicarbonateframe,sulfateframe)
  
  
  
  observeEvent(input$run_button, {
    output$Plot<-renderPlot(
      
      ggplot(alldata, mapping=aes(x=hours, y=conc, color=Chemical)) +
        geom_point() + mytheme + ggtitle("Ion Concentration over Time")
    )
  })
  
  observeEvent(input$run_button, {
    output$ExtraChemicals<-renderPlot(
      
      ggplot(bonusdataframe3(), mapping=aes(x=hours, y=conc, color=name)) +
        geom_point()  + mytheme+ ggtitle("Counter-Ion Concentration over Time")
    )
  })
  
}
    
  
  shinyApp(ui, server)
  
  


