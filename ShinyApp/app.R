library(shiny)
library(ggplot2)  # for the diamonds dataset
library(DT)



ui <- fluidPage(
  

  
  
  title = "Examples of DataTables",
  sidebarLayout(
    sidebarPanel(
      
      fileInput('file1', 'Choose xlsx file',
                accept = c(".xlsx")),
      
      params<-read_excel('file1', sheet = 'params')
      
      conditionalPanel(
        'input.dataset === "Params"'
        #checkboxGroupInput("show_vars", "Columns in diamonds to show:",
        #                   names(diamonds), selected = names(diamonds))
      ),
      conditionalPanel(
        'input.dataset === "Ions"',
        helpText("Click the column header to sort a column.")
      ),
      conditionalPanel(
        'input.dataset === "Cin"',
        helpText("Display 5 records by default.")
      )
    ),
    mainPanel(
      tabsetPanel(
        id = 'dataset',
        tabPanel("Params", DT::dataTableOutput("mytable1")),
        tabPanel("Ions", DT::dataTableOutput("mytable2")),
        tabPanel("Cin", DT::dataTableOutput("mytable3"))
      )
    )
  )
)

server <- function(input, output) {
  
  # choose columns to display
  diamonds2 = diamonds[sample(nrow(diamonds), 1000), ]
  output$mytable1 <- DT::renderDataTable({
    DT::datatable(diamonds2[, input$show_vars, drop = FALSE])
  })
  
  # sorted columns are colored now because CSS are attached to them
  output$mytable2 <- DT::renderDataTable({
    DT::datatable(mtcars, options = list(orderClasses = TRUE))
  })
  
  # customize the length drop-down menu; display 5 rows per page by default
  output$mytable3 <- DT::renderDataTable({
    DT::datatable(iris, options = list(lengthMenu = c(5, 30, 50), pageLength = 5))
  })
  
}

shinyApp(ui, server)