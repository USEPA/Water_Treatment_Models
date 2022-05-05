library(shiny)

ui<-navbarPage(
  title='Ion Exchange Model',
  main_page<-tabPanel(title="Analysis",
                      titlePanel("Analysis"),
                      sidebarLayout(
                        sidebarPanel(
                          title="Inputs",
                          fileInput("csv_input", "Select CSV File to Import", accept=".csv"),
                          checkboxGroupInput("ChemicalList", "Chemicals",
                                             c("Chloride", "Sulfate", "Bicarbonate", "Nitrate")),
                          actionButton("run_button", "Run Analysis", icon=icon("play"))
                        ),
                        mainPanel(
                          tabsetPanel(
                            tabPanel(
                              title="Plot",
                              plotOutput("plot_1")
                            ),
                            tabPanel(
                              title="Statistics",
                              
                            )
                          )
                        )
                      )
                    ),
  about_page<-tabPanel(title="About")
)

server<-function(input,output){
  
  data_input<-reactive({
    req(input$csv_input)
    fread(input$csv_input$datapath)
  })
  
  observeEvent(data_input(), {
    choices<-c(not_sel,names(data_input()))
    updateSelectInput(inputId="ChemicalList", choices=choices)
  })
  
  ChemicalList<-eventReactive(input$run_button, input$ChemicalList)
  
  plot_1<-eventReactive(input$run_button, {
    draw_plot_1(data_input(), ChemicalList())
  })
  
  output$plot_1<-renderPlot(plot_1())
  
}

shinyApp(ui=ui, server=server)
