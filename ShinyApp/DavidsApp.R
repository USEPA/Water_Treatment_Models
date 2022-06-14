library(readxl)
library(shiny)
library(ggplot2)
library(dplyr)

if (interactive()) {
  
  ui <- fluidPage(
    sidebarLayout(
      sidebarPanel(
        fileInput("file1", "Choose CSV File", accept = ".xlsx"),
        actionButton("run_button", "Run Analysis", icon=icon("play")),
      ),
      
      
      mainPanel(
        tabsetPanel(
          tabPanel(
            
            title="Data",
            
            tabsetPanel(
              tabPanel(
                title="Parameters",
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
            
            tabPanel(
              title='Ions',
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
                
                fluidRow(
                  column(1,
                         textOutput("Nitrate")),
                  column(2, offset=1,
                         textOutput("Nitratemw")),
                  column(3,
                         textOutput("NitrateKxA")),
                  column(4,
                         textOutput("NitrateValence"))
                ),
            
          
                
              
             ),
            
            tabPanel(
              title='Inital Concentration',
              fluidRow(
                column(1,
                       textOutput("Name2")),
                column(2, offset=1,
                       textOutput("InitialTime")),
                column(3,
                       textOutput("FinalTime"))
              ),
              
              fluidRow(
                column(1,
                       textOutput("Chloride2")),
                column(2,
                       textOutput("Chlorideti")),
                column(3,
                       textOutput("Chloridetf"))
              ),
              
              fluidRow(
                column(1,
                       textOutput("Sulfate2")),
                column(2,
                       textOutput("Sulfateti")),
                column(3,
                       textOutput("Sulfatetf"))
              ),
              
              fluidRow(
                column(1,
                       textOutput("Bicarbonate2")),
                column(2,
                       textOutput("Bicarbonateti")),
                column(3,
                       textOutput("Bicarbonatetf"))
              ),
              
              fluidRow(
                column(1,
                       textOutput("Nitrate2")),
                column(2,
                       textOutput("Nitrateti")),
                column(3,
                       textOutput("Nitratetf"))
              )
              
              
            )
            )),
          
          
          
          tabPanel(
            title="Analysis",
            plotOutput("Plot"),
            plotOutput("ExtraChemicals")
          ),
          
          
          tabPanel(
            title="Statistics",
            
            fluidRow(
              column(1,
                     tableOutput("Group1")),
              column(2, offset=2,
                     tableOutput("Char1"),
                     tableOutput("Char2"),
                     tableOutput("Char3"))),
            
          br(),
          
            fluidRow(
              column(1,
                     textOutput("Group2")),
              column(2, offset=2,
                     textOutput("Char4"),
                     textOutput("Char5"),
                     textOutput("Char6")))
          
        )
      )
    )
  )
  )
  
  server <- function(input, output) {
    
    params2<-reactive({
      file <- input$file1
      ext <- tools::file_ext(file$datapath)
      
      req(file)
      validate(need(ext == "xlsx", "Please upload a csv file"))
      
      params<-read_xlsx(file$datapath, sheet=1)
    })
    
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
    output$InitialTime<-renderText("Inital Time")
    output$FinalTime<-renderText("Final Time")
    
    
    
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
    
    
    
    
    
    mytheme <-  theme(panel.background = element_rect(fill = "white", colour = NA),
                      panel.grid.major = element_line(colour = "grey70", size = 0.2),
                      panel.grid.minor = element_line(colour = "grey85", size = 0.5),
                      legend.position = "top",
                      legend.title = element_text(colour = "black", size = 8, face = "bold", hjust = 0.5),
                      legend.text = element_text(colour = "black", size = 8),
                      legend.key.size = unit(0.5, "line"),
                      strip.text = element_text(colour = "black", size = 7),
                      axis.ticks = element_line(colour = "black", size = 1),
                      axis.line = element_line(colour = "black", size = 1, lineend = "square"),
                      axis.text.x = element_text(colour = "black", size = 8),
                      axis.text.y = element_text(colour = "black", size = 8),
                      axis.title.x = element_text(colour = "black", size = 8),
                      axis.title.y = element_text(colour = "black", size = 8))
    
    
    
    dat<-data.frame(hours = out[[1]], conc = out[[2]][, liquid_id, 1, outlet_id])
    dat1<-data.frame(hours = out[[1]], conc = out[[2]][, liquid_id, 2, outlet_id])
    dat2<-data.frame(hours = out[[1]], conc = out[[2]][, liquid_id, 3, outlet_id])
    dat3<-data.frame(hours = out[[1]], conc = out[[2]][, liquid_id, 4, outlet_id])
    
    bonusdataframe<-data.frame(hours=c(), conc=c())
    bonusdataframe2<-data.frame(hours=c(), conc=c())
    
    
    for (x in 5:nrow(fulldata)){
      dx_frame<-data.frame(
        hours=out[[1]], conc=out[[2]][, liquid_id, x, outlet_id], chemical=fulldata[x,1]
      )
      bonusdataframe<-rbind(dx_frame, bonusdataframe2)
    }
    
    
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
    
    
    
    output$Plot<-renderPlot(ggplot(alldata, mapping=aes(x=hours, y=conc, color=Chemical)) +
                              geom_point()
    )
    
    output$ExtraChemicals<-renderPlot(ggplot(bonusdataframe, mapping=aes(x=hours, y=conc, color=name)) +
                                        geom_point())
    
  }
  
  shinyApp(ui, server)
}

