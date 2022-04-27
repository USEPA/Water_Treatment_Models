library(shiny)
library(readxl)

runApp(
  list(
    ui = fluidPage(
      titlePanel("Use readxl"),
      sidebarLayout(
        sidebarPanel(
          fileInput('file1', 'Choose xlsx file',
                    accept = c(".xlsx")
          )
        ),
        mainPanel(
          tableOutput('contents'))
      )
    ),
    
    
    server = function(input, output){
    
        output$contents <- renderTable({
        
        req(input$file1)
        
        inFile <- input$file1
        
        read_excel(inFile$datapath, 1)
      })
      
        output$contents <- renderTable({
        
        req(input$file1)
        
        inFile <- input$file1
        
        read_excel(inFile$datapath, 2)
      })
    }
  )
)
