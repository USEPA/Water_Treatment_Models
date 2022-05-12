library(readxl)

if (interactive()) {
  
  ui <- fluidPage(
    sidebarLayout(
      sidebarPanel(
        fileInput("file1", "Choose CSV File", accept = ".xlsx"),
        selectInput("ChemicalList", "Chemicals",
                    c("Chloride" = 1, "Sulfate" = 2, "Bicarbonate" = 3, "Nitrate" = 4)),
        actionButton("run_button", "Run Analysis", icon=icon("play")),
      ),
      mainPanel(
        tabsetPanel(
          tabPanel(
            title="Data",  
            tableOutput("contents")
          ),
          tabPanel(
            title="Analysis",
            plotOutput("Plot")
          )
        )
      )
    )
  )
  
  server <- function(input, output) {
    output$contents <- renderTable({
      file <- input$file1
      ext <- tools::file_ext(file$datapath)
      
      req(file)
      validate(need(ext == "xlsx", "Please upload a csv file"))
      
      read_excel(file$datapath, sheet=1)
    })
    
    
    
    dat<-data.frame(hours = out[[1]], conc = out[[2]][, liquid_id, 1, outlet_id])
    dat1<-data.frame(hours = out[[1]], conc = out[[2]][, liquid_id, 2, outlet_id])
    dat2<-data.frame(hours = out[[1]], conc = out[[2]][, liquid_id, 3, outlet_id])
    dat3<-data.frame(hours = out[[1]], conc = out[[2]][, liquid_id, 4, outlet_id])
    
    
    output$Plot<-renderPlot(ggplot()+
                              geom_point(dat, mapping=aes(x = hours, y = conc))+
                              geom_point(dat1, mapping=aes(x = hours, y = conc), color="blue")+
                              geom_point(dat2, mapping=aes(x = hours, y = conc), color="red")+
                              geom_point(dat3, mapping=aes(x = hours, y = conc), color="green")
                              )
    
          
  }
  
  shinyApp(ui, server)
