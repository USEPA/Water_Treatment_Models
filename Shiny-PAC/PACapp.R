library(reticulate)
library(plotly)
library(shiny)
library(readxl)
library(shinyjs)
library(DataEditR)
library(tidyr)
library(dplyr)
library(stringr)
library(colorBlindness)
library(writexl)
library(shinyalert)

#~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*#
#------------------------------------------------------------------------------#
#unit conversions
#------------------------------------------------------------------------------#
## length
m2cm<-100                               #meters to centimeters
mm2cm<-0.1                              #millimeters to centimeters
cm2cm<-1                                #centimeters to centimeters (for consistency)
in2cm<-2.54                             #inches to centimeters
ft2cm<-12 * in2cm                       #centimeters to feet
## time
sec2sec<-1
min2sec<-60
S_PER_HR <- 60 * 60                     # 
hour2sec<-60 * min2sec
day2sec<-24 * hour2sec
month2sec<-30 * day2sec                 #assumes 30 day month
year2sec<-365.25 * day2sec
month2day<-1/30
day2day<-1
hour2day<-24
second2day<- 60 * 60 * hour2day
year2day<-month2day/12
min2min <- 1
min2hour <- 1/60
## velocity
mpmin2cmps<-m2cm/min2sec                #meters per minute to centimeters per second
ftpmin2cmps<-ft2cm/min2sec              #feet per minute to centimeters per second
mph2cmps<-m2cm/hour2sec                 #meters per hour to centimeters per second
mmin2cms<-m2cm/min2sec
ftmin2cms<-ft2cm/min2sec
gal2ft3<-0.133680555556
gpmpft2cmps<-gal2ft3 * ft2cm / min2sec  #gallons per minute per foot squared
ft2ps2cm2ps<-(ft2cm)^2                  #feet squared per second to centimeters squared per second
m2ps2cm2ps<-(m2cm)^2                    #meters per second squared to centimeters per second squared
in2ps2cm2ps<-(in2cm)^2                  #inches per second squared to centimeters per second squared
ft2pm2cm2ps<-(ft2cm)^2 / (min2sec)      #feet per minute squared to centimeters per second squared
m2min2cm2s<-(m2cm^2) / (min2sec) 
## volume
gal2ml<-3785.411784
mgd2mlps<-1e6 * gal2ml/day2sec          #mgd to ml/sec
l2ml <- 1000.

#~~~~~~~~~~~~~~~~~~~ end unit conversions

#------------------------------------------------------------------------------#
#conversion dictionaries
#------------------------------------------------------------------------------#
##set up dictionaries   ### IF new values are added to drop-downs, must also be added here
length_conv <- c("m"=m2cm, "cm"=cm2cm, "mm"=mm2cm, "in"=in2cm, "ft"=ft2cm)

velocity_conv <- c("cm/s"=cm2cm, "m/s"=m2cm, "m/min"=mpmin2cmps, "m/h"=mph2cmps,
                   "m/hr"=mph2cmps, "in/s"=in2cm, "ft/s"=ft2cm, "ft/min"=ftpmin2cmps,
                   "gpm/ft^2"=gpmpft2cmps)

volumetric_conv <- c("cm^3/s"=cm2cm*min2sec, "m^3/s"=min2sec*m2cm^3, "ft^3/s"=min2sec*ft2cm^3,
                     "mL/s"=min2sec*cm2cm, "L/min"=l2ml, "mL/min"=1,
                     "gpm"=gal2ml, "mgd"=1e6 * gal2ml)


# time_conv <- c("Hours"=hour2day, "Days"=day2day, "Months"=month2day, "Years"=year2day,
#                "hr"=hour2day, "day"=day2day, "month"=month2day, "year"=year2day,
#                "hours"=hour2day, "days"=day2day, "hrs"=hour2day)

time_conv <- c("Minutes"=min2min, "Hours"=min2hour)

kL_conv <- c("ft/s"=ft2cm, "m/s"=m2cm, "cm/s"=cm2cm, "in/s"=in2cm, 
             "m/min"=mpmin2cmps, "ft/min"=ftpmin2cmps, "m/h"=mph2cmps,
             "m/hr"=mph2cmps)
ds_conv <- c("ft^2/s"=ft2ps2cm2ps, "m^2/s"=m2ps2cm2ps, "cm^2/s"=cm2cm,
             "in^2/s"=in2ps2cm2ps)

mass_conv <- c("meq"=1000, "meq/L"=1000, "mg"=1000, "ug"=1, "ng"=1e-3, "mg/L"=1000, "ug/L"=1, "ng/L"=1e-3) ### changed

density_conv<-c("g/ml"=1)

weight_conv<-c("kg"=1000, "g"=1, "lb"=1000/2.204)

formatvector <- c("square")
ldvector <- c("m")
# ldvector <- c("cm", "m", "mm", "in", "ft")
tempvector <- c("C")
heightvector <- c("m")
volvector <- c("L")
# volvector <- c("ml", "cm3", "l", "gal", "ft3", "m3")
flowvector <- c("m3/s")
# flowvector <- c("ml/min", "ml/s", "gpm", "gal/min", "lpm", "l/min", "mgd", "ft3/s", "cfs", "m3/s")
radvector <- c("cm")
HRTvector <- c("min")
CRTvector <- c("min")
dosagevector <- c("mg/L")
denvector <- c("gm/ml")

dosage_range <- range(c(5, 150))

notificationDuration <- 10 # Number of seconds to display the notification

reticulate::source_python("PAC.py")

read_in_files <- function(input, file) {
  tryCatch({
    Properties <- read_excel(file, sheet = 'Contactor')
    write.csv(Properties, 'temp_file/Contactor.csv', row.names=FALSE)
  },
  error=function(err){
    print(err)
    showNotification("Warning: Contactor sheet doesn't exist. Reverting to default values.", duration = notificationDuration, closeButton = TRUE, type = "warning")
  })
  tryCatch({
    Kdata <- read_excel(file, sheet = 'PAC Characteristics')
    write.csv(Kdata, 'temp_file/PAC.csv', row.names=FALSE)
  },
  error=function(err){
    print(err)
    showNotification("Warning: PAC sheet doesn't exist. Reverting to default values.", duration = notificationDuration, closeButton = TRUE, type = "warning")
  })
  tryCatch({
    columnSpecs <- read_excel(file, sheet= 'Compounds')
    write.csv(columnSpecs, 'temp_file/Compounds.csv', row.names=FALSE)
  },
  error=function(err){
    print(err)
    showNotification("Warning: Compounds sheet doesn't exist. Reverting to default values.", duration = notificationDuration, closeButton = TRUE, type = "warning")
  })
  tryCatch({
    name<-read_excel(file, sheet='name')
    write.csv(name, 'temp_file/filename.csv', row.names=FALSE)
  },
  warning=function(war){
  },
  error=function(err){
    tryCatch({
      print(err)
      namedata<-data.frame(name=c(input$file1$name))
      write.csv(namedata, "temp_file/filename.csv", row.names=FALSE)
    }, 
    error=function(e){
      print(e)
      print(file)
      file_name <- data.frame(name=c(file))
      write.csv(file_name, 'temp_file/filename.csv', row.names=FALSE)
    })
  })
}

read_in_files(input, paste0("PAC_config.xlsx"))

ui <- fluidPage(
    # Adds the EPA banner to the model
    HTML("<html lang = 'en'>"),

    tags$body(class = "html wide-template"),
    tags$head(tags$link(rel = "stylesheet",
                        type = "text/css",
                        href = "style.css")),

    # Header
    HTML("<header class='masthead clearfix' role='banner'>
        <img alt='' class='site-logo' src='https://www.epa.gov/sites/all/themes/epa/logo.png'>
        <div class='site-name-and-slogan'>
        <h1 class='site-name'><a href='https://www.epa.gov' rel='home' title='Go to the home page'><span>US EPA</span></a></h1>
        <div class='site-slogan'>
        United States Environmental Protection Agency
        </div>
        </div>
        <div class='region-header'>
        <div class='block-epa-core-gsa-epa-search' id='block-epa-core-gsa-epa-search'>"),

    HTML("</div>
        </div>
        </header>
        <nav class='nav main-nav clearfix' role='navigation'>
        <div class='nav__inner'>
        <h2 class='element-invisible'>Main menu</h2>
        <ul class='menu' role='menu'>
        <li class='expanded active-trail menu-item' role='presentation'>
        <a class='active-trail menu-link' href='https://www.epa.gov/environmental-topics' role='menuitem' title='View links to the most popular pages for each of EPA&#8217s top environmental topics.'>Environmental Topics</a></li>
        <li class='menu-item' role='presentation'>
        <a class='menu-link' href='https://www.epa.gov/laws-regulations' role='menuitem' title='View links to regulatory information by topic and sector, and to top pages about environmental laws, regulations, policies, compliance, and enforcement.'>Laws &amp; Regulations</a></li>
        <li class='expanded menu-item' role='presentation'>
        <a class='menu-link' href='https://www.epa.gov/aboutepa' role='menuitem' title='Learn more about our mission, organization, and locations.'>About EPA</a></li>
        </ul>
        </div>
        </nav>
        <div class='mobile-nav' id='mobile-nav'>
        <div class='mobile-bar clearfix'>
        <label class='menu-button' for='mobile-nav-toggle'>Menu</label>
        </div><input checked id='mobile-nav-toggle' type='checkbox'>
        <div class='mobile-links element-hidden' id='mobile-links' style='height:2404px;'>
        <ul class='mobile-menu'>
        <li class='expanded menu-item'><a class='menu-link' href='https://www.epa.gov/environmental-topics' tabindex='-1' title='View links to the most popular pages for each of EPA&#8217s top environmental topics.'>Environmental Topics</a></li>
        <li class='menu-item'><a class='menu-link' href='https://www.epa.gov/laws-regulations' tabindex='-1' title='View links to regulatory information by topic and sector, and to top pages about environmental laws, regulations, policies, compliance, and enforcement.'>Laws &amp; Regulations</a></li>
        <li class='expanded menu-item'><a class='menu-link' href='https://www.epa.gov/aboutepa' tabindex='-1' title='Learn more about our mission, organization, and locations.'>About EPA</a></li>
        </ul>
        </div>
        </div>
        <section class='main-content clearfix' id='main-content' lang='en' role='main' tabindex='-1'>
        <div class='region-preface clearfix'>
        <div class='block-views-revision-hublinks-block' id='block-views-revision-hublinks-block'>
        <div class='view view-revision-hublinks view-id-revision_hublinks'>
        <span class='related-info'><strong>Related Topics:</strong></span>
        <ul class='menu pipeline'>
        <li class='menu-item'><a href='https://www.epa.gov/environmental-topics'>Environmental Topics</a></li>
        </ul>
        </div>
        </div>
        <div class='block block-pane block-pane-epa-web-area-connect' id='block-pane-epa-web-area-connect'>
        <ul class='menu utility-menu'>
        <li class='menu-item'><a class='menu-link' href='https://www.epa.gov/water-research/forms/contact-us-about-water-research'>Contact Us</a></li>
        </ul>
        </div>
        </div>
        <div class='main-column clearfix'><!--googleon:all-->
        <h1  class='page-title'>Powdered Activated Carbon Model</h1>
        <div class='panel-pane pane-node-content'>
        <div class='pane-content'>
        <div class='node node-page clearfix view-mode-full'>"),
    ####Added from EPA template######################################################

    tags$head(
    tags$style(HTML('.navbar-default .navbar-nav > li > a:hover, .navbar-default .navbar-nav > li > a:focus {
        color: #000; /*Sets the text hover color on navbar*/
    }

    .navbar-default .navbar-nav > .active > a, .navbar-default .navbar-nav > .active >
        a:hover, .navbar-default .navbar-nav > .active > a:focus {
        color: white; /*BACKGROUND color for active*/
            background-color: #0e6cb6;
        }

    .navbar-default {
        background-color: #0e6cb6;
        border-color: #030033;
    }

    .navbar-nav > li > a, .navbar-brand {
        padding-top:15px !important;
        padding-bottom:0 !important;
        height: 25px;
    }
    .navbar {min-height:25px !important;}


    .navbar-default .navbar-nav > li > a {
        color: white; /*Change active text color here*/
    }'))),

    tags$style(HTML("
        .tabbable > .nav > li > a                  {background-color: #D3D3D3;  color:black}
    # ")),

    useShinyjs(),

    navbarPage("", id = "inTabset", # Allows for automatic switching between tab panels
        tabPanel("Input",
            sidebarLayout(
                sidebarPanel(
                    fileInput("file1", "Choose .xlsx File", accept = ".xlsx"),
                    tableOutput("selectedfile"),

                    br(),

                    sliderInput("nrv", "Radial Collocation Points",3, 18, 7),

                    br(),
                    
                    actionButton("run_button", "Run Analysis", icon=icon("play")),
                    
                    br(), br(),
                    
                    actionButton("Stop", "Stop App", icon=icon("square"), style="color: #000000; background-color: #ff0000; border-color: #e60000")
                ),                     
                mainPanel(
                    tabsetPanel(
                        tabPanel("Column Parameters",
                    
                            br(), br(),

                            fluidRow(
                                # This radio button toggles volume
                                column(3, HTML(paste0("<h4>","<strong>", "Contactor", "</strong>", "</h4>"))),
                                column(3, shinyWidgets::autonumericInput(
                                    inputId = "temp",
                                    label="Temperature",
                                    value = "10",
                                    currencySymbolPlacement = "p",
                                    decimalPlaces = 3,
                                    digitGroupSeparator = ",",
                                    decimalCharacter = "."
                                )),
                                column(3, selectInput("tempunits", "Temperature Units", c("C")))
                            ),
                        
                            fluidRow(
                                column(3, radioButtons("volselect", "", c("Dimensions", "Volume"))),                                        
                                column(3,
                                    selectInput("format", "Format", list('Square'))
                                ),
                                column(3,)
                            ),

                            fluidRow(
                                column(3,),                                        
                                column(3, shinyWidgets::autonumericInput(
                                    inputId = "ld",
                                    label="Length/diameter",
                                    value = 10,
                                    decimalPlaces = 4,
                                    digitGroupSeparator = ",",
                                    decimalCharacter = "."
                                )),
                                column(3, selectInput("ldunits", "Length/diameter Units", c("m")))
                            ),
                        
                            fluidRow(
                                column(3,),
                                column(3, shinyWidgets::autonumericInput(
                                    inputId = "height",
                                    label="Height",
                                    value = 5,
                                    currencySymbolPlacement = "p",
                                    decimalPlaces = 3,
                                    digitGroupSeparator = ",",
                                    decimalCharacter = "."
                                )),
                                column(3, selectInput("heightunits", "Height Units", c("m")))
                            ),
                        
                            fluidRow(
                                column(3,),
                                column(3, shinyWidgets::autonumericInput(
                                    inputId = "vol",
                                    label="Volume",
                                    value = 500000,
                                    currencySymbolPlacement = "p",
                                    decimalPlaces = 3,
                                    digitGroupSeparator = ",",
                                    decimalCharacter = "."
                                )),
                                column(3, selectInput("volunits", "Volume Units", c("L")))
                            ),

                            fluidRow(
                                column(3,),
                                column(3, shinyWidgets::autonumericInput(
                                    inputId = "flow",
                                    label="Flow",
                                    value = 5,
                                    currencySymbolPlacement = "p",
                                    decimalPlaces = 3,
                                    digitGroupSeparator = ",",
                                    decimalCharacter = "."
                                )),
                                column(3, selectInput("flowunits", "Flow Units", c("m3/s")))
                            ),

                            fluidRow(
                                column(3,),
                                column(3,
                                    sliderInput("hrt", "HRT", 0, 300, 0)
                                ),
                                column(3, selectInput("HRTunits","HRT Units",c("min"))),
                            ),

                            fluidRow(
                                column(3,),
                                column(3,
                                    sliderInput("crt", "CRT", 0, 300, 0),
                                ),
                                column(3, selectInput("CRTunits","CRT Units",c("min"))),
                            ),

                            fluidRow(
                                column(3,),
                                column(3,
                                    sliderInput("dosage", "PAC Dosage", 0, 20, 0),
                                ),
                                column(3, selectInput("dosageunits","Dosage Units",c("mg/L")))
                            ),

                            hr(),

                            fluidRow(
                                column(3, HTML(paste0("<h4>","<strong>", "PAC Characteristics", "</strong>", "</h4>"))),
                                column(3, shinyWidgets::autonumericInput(
                                    inputId = "den",
                                    label="Density",
                                    value = 0.8034,
                                    currencySymbolPlacement = "p",
                                    decimalPlaces = 4,
                                    digitGroupSeparator = ",",
                                    decimalCharacter = "."
                                )),
                                column(3, selectInput("denunits", "Density Units", c("gm/ml")))
                            ),

                            fluidRow(
                                column(3,),
                                column(3, shinyWidgets::autonumericInput(
                                    inputId = "por",
                                    label="Porosity",
                                    value = 0.641,
                                    decimalPlaces = 3,
                                    digitGroupSeparator = ",",
                                    decimalCharacter = "."
                                ))
                            ),
                            
                            fluidRow(
                                column(3,),
                                column(3, shinyWidgets::autonumericInput(
                                    inputId = "rad",
                                    label="Radius",
                                    value = 0.00103,
                                    decimalPlaces = 5,
                                    digitGroupSeparator = ",",
                                    decimalCharacter = "."
                                )),
                                column(3, selectInput("radunits","Radius Units",c("cm")))
                            ),
                        ),
                        tabPanel("Compounds",

                            br(), br(),

                            h4("Compound List"),
                            dataEditUI("edit-1") 
                        )
                    )
                )
            )                
        ),

        tabPanel("Concentration Output",
            sidebarLayout(
                sidebarPanel(
                    selectInput("OCunits", "Output Concentration Units", c("mg/L", "ug/L")),
                    selectInput("timeunits","Output Time Units",c("Minutes", "Hours")),

                    br(), br(), br(),
          
                    downloadButton("save_button1", "Save Data")
                ),
                mainPanel(
                    br(), br(),
                    
                    shinycssloaders::withSpinner(plotlyOutput("plot1"))
                )
            ),
        ),

        tabPanel(HTML("Dosage Calculator</a></li><li><a href='https://github.com/USEPA/Water_Treatment_Models/tree/master/PAC' target='_blank'>Help"),
            sidebarLayout(
                sidebarPanel(
                    sliderInput("dosagerange", label = "Dosage Range", min = dosage_range[1], max = dosage_range[2], value = dosage_range, step = 5),
                    sliderInput("dosageinterval", "Dosage Intervals", 3, 7, 5),
                    selectInput("OCunits2", "Output Concentration Units", c("mg/L", "ug/L")),

                    br(),

                    shinyWidgets::autonumericInput(inputId = "target", label="Target Concentration", value = 4.0, currencySymbolPlacement = "p", decimalPlaces = 1, digitGroupSeparator = ",", decimalCharacter = "."),
                    selectInput("targetunits", "Target Units", c("ng")),

                    actionButton("calculate_by_target", "Calculate", icon=icon("play")),
                    
                    br(), br(), br(),

                    downloadButton("save_button2", "Save Data")
                ),                     

                mainPanel(
                    tabsetPanel(
                        tabPanel("Concentration by Dosage",
                            br(), br(),

                            fluidRow(
                                column(3, selectInput("compound", "Compound", c("MIB"))),
                                column(3,),
                                column(3,)
                            ),
                            
                            # Suppress Shiny errors
                            tags$style(type="text/css",
                            ".shiny-output-error { visibility: hidden; }",
                            ".shiny-output-error:before { visibility: hidden; }"
                            ),

                            hr(),
                            
                            shinycssloaders::withSpinner(plotlyOutput("plot2")),
                        ),
                        tabPanel("HRT for Target",
                            br(), br(),

                            shinycssloaders::withSpinner(plotlyOutput("plot3"))
                        )
                    )
                )
            )
        )
    )
)

# Define server logic ----
server <- function(input, output, session) {
    observeEvent(input$file1, {
        file <- input$file1
        read_in_files(input, paste0(file$datapath))
    })
    
    observeEvent(input$file1, {
        session$reload()
    })
    
    observeEvent(input$Stop, {
        stopApp()
    })
    
    fileuploadedname <- read.csv("temp_file/filename.csv")
    output$selectedfile <- renderTable(fileuploadedname)

    file_direc <- paste(getwd(), '/temp_file/', sep = '')
    contactor<-reactiveVal(read.csv(paste(file_direc, "Contactor.csv", sep = '')))
    pac<-reactiveVal(read.csv(paste(file_direc, "PAC.csv", sep = '')))

    test_df <- data.frame(C = c('format', 'length/diameter', 'height', 'volume'))
    flags <- reactive({test_df$C %in% contactor()$name})

    observe({
        # if (flags()[1] & flags()[2] & flags()[3]) {
        if (flags()[1]) {
            updateSelectInput(session, "format", choices = unique(c(filter(contactor(), name == 'format')$value, formatvector)))
            updateNumericInput(session, "ld", value = filter(contactor(), name == 'length/diameter')$value)
            updateSelectInput(session, "ldunits", choices = unique(c(filter(contactor(), name == 'length/diameter')$units, ldvector)))
            updateNumericInput(session, "height", value = filter(contactor(), name == 'height')$value)
            updateSelectInput(session, "heightunits", choices = unique(c(filter(contactor(), name == 'height')$units, heightvector)))
            updateRadioButtons(session, "volselect", selected = "Dimensions")
        } else if (flags()[4]) {                         
            updateNumericInput(session, "vol", value = filter(contactor(), name == 'volume')$value)
            updateSelectInput(session, "vol", choices = unique(c(filter(contactor(), name == 'volume')$units, volvector)))
            updateRadioButtons(session, "volselect", selected = "Volume")
        }
    })

    observe({
        toggleState("format", condition = input$volselect != "Volume")
        toggleState("ld", condition = input$volselect != "Volume")
        toggleState("ldunits", condition = input$volselect != "Volume")
        toggleState("height", condition = input$volselect != "Volume")
        toggleState("heightunits", condition = input$volselect != "Volume")
        toggleState("vol", condition = input$volselect != "Dimensions")
        toggleState("volunits", condition = input$volselect != "Dimensions")
    })

    # format <- reactive({filter(contactor(), name == "format")$value})
    # lengthdiameter <- reactive({filter(contactor(), name == "length/diameter")$value})
    temperature <- reactive({filter(contactor(), name == "temperature")$value})
    # height <- reactive({filter(contactor(), name == 'height')$value})
    # volume <- reactive({filter(contactor(), name == 'volume')$value})
    flow <- reactive({filter(contactor(), name == 'flow')$value})
    HRT <- reactive({filter(contactor(), name == 'HRT')$value})
    CRT <- reactive({filter(contactor(), name == 'CRT')$value})
    dosage <- reactive({filter(contactor(), name == 'PAC Dosage')$value})

    dens <- reactive({filter(pac(), name == 'density')$value})
    porosity <- reactive({filter(pac(), name == 'porosity')$value})
    radius <- reactive({filter(pac(), name == 'radius')$value})

    nrv <- reactive(7)
    
    # formatvec <- reactive({
    #     formatv <- c(filter(contactor(), name == 'format')$value, formatvector)

    #     return(unique(formatv))
    # })

    # ldvec <- reactive({
    #     ldv <- c(filter(contactor(), name == 'length/diameter')$units, ldvector)

    #     return(unique(ldv))
    # })
    
    tempvec <- reactive({
        tempv <- c(filter(contactor(), name == 'temperature')$units, tempvector)

        return(unique(tempv))
    })
    
    # heightvec <- reactive({
    #     heightv <- c(filter(contactor(), name == 'height')$units, heightvector)

    #     return(unique(heightv))
    # })

    # volvec <- reactive({
    #     volv <- c(filter(contactor(), name == 'volume')$units, volvector)

    #     return(unique(volv))
    # })
    
    flowvec <- reactive({
        flowv <- c(filter(contactor(), name == 'flow')$units, flowvector)

        return(unique(flowv))
    })
    
    denvec <- reactive({
        denv <- c(filter(pac(), name == 'density')$units, denvector)

        return(unique(denv))
    })
    
    radvec <- reactive({
        radv <- c(filter(pac(), name == 'radius')$units, radvector)

        return(unique(radv))
    })

    HRTvec <- reactive({
        HRTv <- c(filter(contactor(), name == 'HRT')$units, HRTvector)

        return(unique(HRTv))
    })

    CRTvec <- reactive({
        CRTv <- c(filter(contactor(), name == 'CRT')$units, CRTvector)

        return(unique(CRTv))
    })

    dosagevec <- reactive({
        dosagev <- c(filter(contactor(), name == 'PAC Dosage')$units, dosagevector)

        return(unique(dosagev))
    })

    compoundvec <- reactive({
        compounds <- compounddat()[,-1]
        compoundv <- colnames(compounds)
    })
  
    #------------------------------------------------------------------------------#
            #Updating default values with the values that were uploaded#
    #------------------------------------------------------------------------------#    
    observe({
        # updateNumericInput(session, "ld", value = format(lengthdiameter(), digits = 4, scientific = FALSE))
        updateNumericInput(session, "temp", value = temperature())
        # updateNumericInput(session, "height", value = format(height(), digits = 4, scientific = FALSE))
        # updateNumericInput(session, "vol", value = format(volume(), digits = 4, scientific = FALSE))
        updateNumericInput(session, "flow", value = format(flow(), digits = 4, scientific = FALSE))
        updateNumericInput(session, "den", value = dens())
        updateNumericInput(session, "por", value = porosity())
        updateNumericInput(session, "rad", value = format(radius(), digits = 4, scientific = FALSE))
        updateNumericInput(session, "nrv", value = nrv())

        updateSliderInput(session, "hrt", value = HRT())
        updateSliderInput(session, "crt", value = CRT())
        updateSliderInput(session, "dosage", value = dosage())
        
        # updateSelectInput(session, "format", choices = formatvec())
        # updateSelectInput(session, "ldunits", choices = ldvec())
        updateSelectInput(session, "tempunits", choices = tempvec())
        # updateSelectInput(session, "heightunits", choices = heightvec())
        # updateSelectInput(session, "volunits", choices = volvec())
        updateSelectInput(session, "flowunits", choices = flowvec())
        updateSelectInput(session, "denunits", choices = denvec())
        updateSelectInput(session, "radunits", choices = radvec())
        updateSelectInput(session, "HRTunits", choices = HRTvec())
        updateSelectInput(session, "CRTunits", choices = CRTvec())
        updateSelectInput(session, "dosageunits", choices = dosagevec())
        updateSelectInput(session, "compound", choices = compoundvec())
    })

    compounddat <- dataEditServer("edit-1",  data = paste(file_direc, 'Compounds.csv', sep = ''))
    dataOutputServer("output-1", data = compounddat)

    pac_obj <- reactiveVal(data.frame())
    sub_data <- reactiveVal(data.frame())
    sub_data2 <- reactiveVal(data.frame())

    observeEvent(input$run_button, {
        # Pass inputs to Python helper
        contactor_df <- contactor()[,-1]
        rownames(contactor_df) <- contactor()[,1]
        pac_df <- pac()[,-1]
        rownames(pac_df) <- pac()[,1]
        compounds_df <- compounddat()[,-1]
        rownames(compounds_df) <- compounddat()[,1]
        PAC_instance <- PAC_CFPSDM(contactor_df, pac_df, compounds_df, input$nrv)

        # Process concentration output
        df <- as.data.frame(PAC_instance$run_PAC_PSDM())
        df$time <- as.numeric(rownames(df))
        df <- pivot_longer(df, cols = !time, names_to = "name", values_to = "conc")
        df$conc <- as.numeric(unlist(df$conc))
        df <- df[order(factor(df$name, levels = colnames(compounddat()))), ]
        pac_obj(df)

        showNotification("Starting model run.", duration = notificationDuration, closeButton = TRUE, type = "message")
        updateTabsetPanel(session, "inTabset", selected = "Concentration Output")
    })

    observeEvent(input$calculate_by_target, {
        # Pass inputs to Python helper
        contactor_df <- contactor()[,-1]
        rownames(contactor_df) <- contactor()[,1]
        pac_df <- pac()[,-1]
        rownames(pac_df) <- pac()[,1]
        compounds_df <- compounddat()[,-1]
        rownames(compounds_df) <- compounddat()[,1]
        PAC_instance <- PAC_CFPSDM(contactor_df, pac_df, compounds_df, input$nrv)

        # Calculate dosage intervals and set target
        vector <- numeric()
        index <- 1
        interval <- (input$dosagerange[2] - input$dosagerange[1])/(input$dosageinterval - 1)
        for (i in seq(input$dosagerange[1], input$dosagerange[2], by = interval)) {
            vector[index] <- i
            index <- index + 1 
        }
        data_dict <- PAC_instance$run_multi_dosage(vector)
        target_HRT <- c(30, 60, 90, 120)

        # Process concentration by dosage output
        processed_dict <- PAC_instance$multi_dosage_analyzer(data_dict, target_HRT)
        df <- as.data.frame(processed_dict)
        df <- cbind(dosage = as.numeric(rownames(df)), df)
        df <- data.frame(lapply(df, unlist))
        names(df)[names(df) == "dosage"] <- "dosage (mg/L)"
        df <- df %>% mutate(across(!`dosage (mg/L)`, ~ .x / mass_conv[input$OCunits2]))
        df <- df %>% rename_with(
            .fn = ~ str_replace(., "^(.*)\\.(\\d+)\\.0$", "HRT: \\2 min \\1"),
            .cols = !all_of("dosage (mg/L)")
        )
        sub_data(df)

        # Process HRT for target output
        df2 <- as.data.frame(PAC_instance$`_R_HRT_calculator_for_dosage`(input$target, conc_units=input$targetunits))
        df2$`dosage (mg/L)` <- as.numeric(rownames(df2))
        df2 <- pivot_longer(df2, cols = !`dosage (mg/L)`, names_to = "name", values_to = "HRT to below Target (Minutes)")
        df2$`HRT to below Target (Minutes)` <- as.numeric(unlist(df2$`HRT to below Target (Minutes)`))
        df2 <- df2[order(factor(df2$name, levels = colnames(compounddat()))), ]
        sub_data2(df2)
    })

    # Apply conversions and rename columns
    pac_obj_processed <- reactive({
        df <- data.frame(
            time = pac_obj()$time * time_conv[input$timeunits],
            conc = (pac_obj()$conc / 1000) / mass_conv[input$OCunits],
            name = pac_obj()$name
        )

        colnames(df) <- c(paste0("time (", input$timeunits, ")"), paste0("conc (", input$OCunits, ")"), "name")
        
        df
    })
    sub_data_processed <- reactive({
        if (ncol(sub_data()) > 0) {
            df <- sub_data()
            df %>% mutate(across(!`dosage (mg/L)`, ~ .x / mass_conv[input$OCunits2]))
            df$`conc units` <- input$OCunits2
            return(df)
        } else {
            return(data.frame())
        }
    })

 
    # Prepare concentration by dosage plot
    HRT_obj <- reactive(data.frame(
        dosage = unname(sub_data_processed()["dosage (mg/L)"]),
        HRT30 = unname(sub_data_processed()[paste0("HRT: 30 min ", input$compound)]),
        HRT60 = unname(sub_data_processed()[paste0("HRT: 60 min ", input$compound)]),
        HRT90 = unname(sub_data_processed()[paste0("HRT: 90 min ", input$compound)]),
        HRT120 = unname(sub_data_processed()[paste0("HRT: 120 min ", input$compound)])
    ))

    # Format plots
    p1 <- reactive({plot_ly(pac_obj_processed(), x = ~get(colnames(pac_obj_processed())[1]), y = ~get(colnames(pac_obj_processed())[2]), color = ~name, type = 'scatter', mode = 'lines') %>% layout(title = "Concentration over Time", showlegend = TRUE,
                                                 legend = list(orientation = 'h', x=0.5, y=1), hovermode = 'x unified',
                                                 xaxis = list(title=paste0("Time (", input$timeunits, ")"), gridcolor = 'ffff'),
                                                 yaxis = list(title=paste0("Concentration (", input$OCunits, ")"), rangemode = "tozero"))
    })
    p2 <- reactive({plot_ly(HRT_obj(), x = ~dosage, y = ~HRT30, type = 'scatter', mode = 'lines+markers', name = paste0("HRT: 30 min ", input$compound)) %>% layout(title = input$compound, showlegend = TRUE,
                                       legend = list(orientation = 'h', y=1), hovermode = 'x unified',
                                       xaxis = list(title="Dosage (mg/L)", gridcolor = 'ffff'),
                                       yaxis = list(title=paste0("Concentration (", input$OCunits2, ")"), rangemode = "tozero")) %>%
                                       add_trace(data = HRT_obj(), x = ~dosage, y = ~HRT60, type = 'scatter', mode = 'lines+markers', name = paste0("HRT: 60 min ", input$compound)) %>%
                                       add_trace(data = HRT_obj(), x = ~dosage, y = ~HRT90, type = 'scatter', mode = 'lines+markers', name = paste0("HRT: 90 min ", input$compound)) %>%
                                       add_trace(data = HRT_obj(), x = ~dosage, y = ~HRT120, type = 'scatter', mode = 'lines+markers', name = paste0("HRT: 120 min ", input$compound))
    })
    p3 <- reactive({plot_ly(sub_data2(), x = ~`dosage (mg/L)`, y = ~`HRT to below Target (Minutes)`, color = ~name, type = 'scatter', mode = 'lines+markers') %>% layout(title = "10.0 ng/L Target - Geosmin", showlegend = TRUE,
                                         legend = list(orientation = 'h', y=1), hovermode = 'x unified',
                                         xaxis = list(title="Dosage (mg/L)", gridcolor = 'ffff'),
                                         yaxis = list(title="HRT to below Target (Minutes)", rangemode = "tozero"))
    })

    # Output plots
    output$plot1 <- renderPlotly(p1())
    output$plot2 <- renderPlotly(p2())
    output$plot3 <- renderPlotly(p3())

    # Save data
    callDownloadHandler <- function() {
        downloadHandler(
            filename = function() {
            paste("data-", Sys.Date(), ".xlsx", sep="")
            },
            content=function(file) {
                sheets <- list("Contactor" = contactor(),
                            "PAC Characteristics" = pac(),
                            "Compounds" = compounddat(),
                            "Concentration Output" = pac_obj_processed(),
                            "Concentration by Dosage" = sub_data_processed(),
                            "HRT for Target" = sub_data2()
                )

                write_xlsx(sheets, file)
            }
        )
    }
    output$save_button1 <- callDownloadHandler()
    output$save_button2 <- callDownloadHandler()
}

shinyApp(ui = ui, server = server)