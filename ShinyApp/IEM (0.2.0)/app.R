library(readxl)
library(shiny)
library(ggplot2)
library(dplyr)
library(xlsx)
library(shinythemes)
library(deSolve)
library(orthopolynom)


ui <- fluidPage(theme=shinytheme("united"),
                
                sidebarLayout(
                  sidebarPanel(
                    fileInput("file1", "Choose CSV File", accept = ".xlsx"),
                    actionButton("apply_button", "Apply Data"),
                    actionButton("run_button", "Run Analysis", icon=icon("play")),
                    tableOutput("summary"),
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
                                                            
                                                            textInput("Qv", ""),
                                                            
                                                            textInput("rbv", ""),
                                                            
                                                            textInput("EBEDv", "")),
                                                     
                                                     column(4,
                                                            selectInput("qunits", "", c("meq/L")),
                                                            selectInput("rbunits", "", c("cm")),
                                                            selectInput("EBEDunits", "", c("N/A")))
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
                                                            textInput("Lv", ""),
                                                            
                                                            textInput("Vv", ""),
                                                            
                                                            textInput("Dv", ""),
                                                            
                                                            textInput("Fv", "")),
                                                     column(4,
                                                            selectInput("LengthUnits", "", c("cm")),
                                                            selectInput("velocityunits", "", c("cm/s")),
                                                            selectInput("DiameterUnits","",c("cm^2")),
                                                            selectInput("flowrateunits","",c("cm^2/s")))
                                                     
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
                                                            textInput("kLv", ""),
                                                            br(), br(),
                                                            textInput("Dsv", "")),
                                                     column(4,
                                                            br(),
                                                            selectInput("filmunits","",c("cm/s")),
                                                            br(), br(),
                                                            selectInput("diffusionunits","",c("cm^2/s")))
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
                                                            textInput("nrv", ""),
                                                            br(), br(),
                                                            textInput("nzv", "")),
                                                     column(4,
                                                            selectInput("radialunits", "", c("N/A")),
                                                            br(), br(),
                                                            selectInput("axialunits","",c("N/A")))),
                                                   
                                                   
                                                   br(), br(), br(),
                                                   
                                                   
                                                   #Parameters Row 5#
                                                   
                                                   fluidRow(
                                                     br(),
                                                     column(1,
                                                            textOutput("Time")),
                                                     column(2, offset=1,
                                                            
                                                            textOutput("TS")),
                                                     column(3, 
                                                            textInput("tv", "")),
                                                     column(4,
                                                            selectInput("timeunits","",c("hr")))),
                                                   
                                                   br(), br(),
                                                   actionButton("save_data", "Save Data", offset=2),
                                                   br(), br()
                                                   
                                                   
                                                   
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
  
  cin2<-reactive({
    file <- input$file1
    ext <- tools::file_ext(file$datapath)
    
    req(file)
    validate(need(ext == "xlsx", "Please upload a csv file"))
    
    cin2<-read_xlsx(file$datapath, sheet=3)
    
    cin2
  })
  
  capacity<-eventReactive(input$apply_button, {
    cap<-filter(params2(), name=="Q")$value
    cap
  })
  
  eebed<-eventReactive(input$apply_button, {
    val<-filter(params2(), name=="EBED")$value
    val
  })
  
  length<-eventReactive(input$apply_button, {
    val<-filter(params2(), name=="L")$value
    val
  })
  
  velocity<-eventReactive(input$apply_button, {
    val<-filter(params2(), name=="v")$value
    val
  })
  
  beadradius<-eventReactive(input$apply_button, {
    val<-filter(params2(), name=="rb")$value
    val
  })
  
  
  film<-eventReactive(input$apply_button, {
    val<-filter(params2(), name=="kL")$value
    val
  })
  
  diffuse<-eventReactive(input$apply_button, {
    val<-filter(params2(), name=="Ds")$value
    val
  })
  
  
  radial<-eventReactive(input$apply_button, {
    val<-filter(params2(), name=="nr")$value
    val
  })
  
  axial<-eventReactive(input$apply_button, {
    val<-filter(params2(), name=="nz")$value
    val
  })
  
  time<-eventReactive(input$apply_button, {
    val<-filter(params2(), name=="time")$value
    val
  })
  
  observe({
    updateTextInput(session, "Vv", value=velocity())
  })
  
  observe({
    updateTextInput(session, "rbv", value=beadradius())
  })
  
  
  observe({
    updateTextInput(session, "Qv", value=capacity())
  })
  
  observe({
    updateTextInput(session, "EBEDv", value=eebed())
  })
  
  observe({
    updateTextInput(session, "Lv", value=length())
  })
  
  observe({
    updateTextInput(session, "kLv", value=film())
  })
  
  observe({
    updateTextInput(session, "Dsv", value=diffuse())
  })
  
  observe({
    updateTextInput(session, "nrv", value=radial())
  })
  
  observe({
    updateTextInput(session, "nzv", value=axial())
  })
  
  observe({
    updateTextInput(session, "tv", value=time())
  })
  
  
  newdataframe<-eventReactive(input$save_data, {
    saveddataframe<-data.frame(
      name=c("Q", "EBED", "L", "v", "rb", "kL", "Ds", "nr", "nz", "time"),
      value=c(input$Qv, input$rbv, input$Lv, input$Vv, input$rbv, input$kLv, input$Dsv, input$nrv, input$nzv, input$tv),
      units=c(input$qunits, input$rbunits, input$EBEDunits, input$LengthUnits, input$velocityunits, input$filmunits, input$diffusionunits, input$radialunits, input$axialunits, input$timeunits)
    )
    saveddataframe
  })
  
  
  
  output$summary<-renderTable(newdataframe())
  
  
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




