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
                     textOutput("tunit")))
              
            
            
            ),
          
          
          
          tabPanel(
            title="Analysis",
            plotOutput("Plot"),
            plotOutput("ExtraChemicals")
          ),
          
          
          tabPanel(
            title="Statistics",
            
            fluidRow(
              column(1,
                     textOutput("Group1")),
              column(2, offset=2,
                     textOutput("Char1"),
                     textOutput("Char2"),
                     textOutput("Char3"))),
            
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
    
    
    
    
    dat<-data.frame(hours = out[[1]], conc = out[[2]][, liquid_id, 1, outlet_id])
    dat1<-data.frame(hours = out[[1]], conc = out[[2]][, liquid_id, 2, outlet_id])
    dat2<-data.frame(hours = out[[1]], conc = out[[2]][, liquid_id, 3, outlet_id])
    dat3<-data.frame(hours = out[[1]], conc = out[[2]][, liquid_id, 4, outlet_id])
    
    #bonusdat<-data.frame(hours = out[[1]], conc = out[[2]][, liquid_id, 5, outlet_id])
    
    #inputdat<-data.frame(hours = out[[1]], conc = out[[2]][, liquid_id, selection, outlet_id])
    
    
    output$Plot<-renderPlot(ggplot()+
                              geom_point(dat, mapping=aes(x = hours, y = conc))+
                              geom_point(dat1, mapping=aes(x = hours, y = conc), color="blue")+
                              geom_point(dat2, mapping=aes(x = hours, y = conc), color="red")+
                              geom_point(dat3, mapping=aes(x = hours, y = conc), color="green")
    )
    
    output$ExtraChemicals<-renderPlot(ggplot()+
                                        geom_point(bonusdat, mapping=aes(x = hours, y = conc)))
    
  }
  
  shinyApp(ui, server)
}

