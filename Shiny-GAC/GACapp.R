library(reticulate)
library(plotly)
library(shiny)
library(readxl)
library(shinyjs)
library(DataEditR)
library(tidyr)
## Commented out for running locally
# renv::install("bioconductor-source/BiocVersion")## needed for colorBlindness on remote
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


time_conv <- c("Hours"=hour2day, "Days"=day2day, "Months"=month2day, "Years"=year2day,
               "hr"=hour2day, "day"=day2day, "month"=month2day, "year"=year2day,
               "hours"=hour2day, "days"=day2day, "hrs"=hour2day)

kL_conv <- c("ft/s"=ft2cm, "m/s"=m2cm, "cm/s"=cm2cm, "in/s"=in2cm, 
             "m/min"=mpmin2cmps, "ft/min"=ftpmin2cmps, "m/h"=mph2cmps,
             "m/hr"=mph2cmps)
ds_conv <- c("ft^2/s"=ft2ps2cm2ps, "m^2/s"=m2ps2cm2ps, "cm^2/s"=cm2cm,
             "in^2/s"=in2ps2cm2ps)

mass_conv <- c("meq"=1000, "meq/L"=1000, "mg"=1000, "ug"=1, "ng"=1e-3, "mg/L"=1000, "ug/L"=1, "ng/L"=1e-3) ### changed

density_conv<-c("g/ml"=1)

weight_conv<-c("kg"=1000, "g"=1, "lb"=1000/2.204)

prvector<-c("cm", "m", "mm", "in", "ft")
pdvector<-c("g/ml")
advector<-c("g/ml")
lengthvector<-c("cm", "m", "mm", "in", "ft")
velocityvector<-c("cm/s", "m/s", "m/min", "m/h", "in/s","ft/s","ft/min", "gpm/ft^2")
timevector <- c("hrs","days")
flowratevector<-c("cm^3/s", "m^3/s", "ft^3/s", "mL/s", "L/min", "mL/min", "gpm", "mgd")
diametervector<-c("cm", "m", "mm", "in", "ft")

weightvector<-c("kg", "g", "lb")
concentrationvector<-c("ug", "ng", "mg")
wfoulingvector <- c("Organic Free", "Rhine", "Portage", "Karlsruhe", "Wausau", "Houghton") # Used to store accepted water types
cfoulingvector <- c("halogenated alkenes", "halogenated alkanes", "halogenated alkanes QSPR", "trihalo-methanes", "aromatics", "nitro compounds", "chlorinated hydrocarbon", "phenols", "PNAs", "pesticides", "PFAS") # Used to store accepted chemical types

notificationDuration <- 10 # Number of seconds to display the notification

reticulate::source_python("GAC_Shiny_helper.py")

#------------------------------------------------------------------------------#
                              #read_in_files
#reads in a file that the user selects, reads pages 'Properties', 'Kdata',
#'columnSpecs', 'dat'. The influent and effluent data get separated and then
#pivoted to a data frame shape that is much easier to use with plotly
#------------------------------------------------------------------------------#
read_in_files <- function(input, file) {
  # Attempts to read-in sheets from Excel file, if sheet doesn't exist it reverts to default values
  tryCatch({
    Properties <- read_excel(file, sheet = 'Properties')
    write.csv(Properties, 'temp_file/Properties.csv', row.names=FALSE)
  },
  error=function(err){
    print(err)
    showNotification("Warning: Properties sheet doesn't exist. Reverting to default values.", duration = notificationDuration, closeButton = TRUE, type = "warning")
  })
  tryCatch({
    Kdata <- read_excel(file, sheet = 'Kdata')
    write.csv(Kdata, 'temp_file/Kdata.csv', row.names=FALSE)
    write.csv(Kdata, 'temp_file/Kdata2.csv', row.names=FALSE)
  },
  error=function(err){
    print(err)
    showNotification("Warning: Kdata sheet doesn't exist. Reverting to default values.", duration = notificationDuration, closeButton = TRUE, type = "warning")
  })
  tryCatch({
    columnSpecs <- read_excel(file, sheet= 'columnSpecs')
    write.csv(columnSpecs, 'temp_file/columnSpecs.csv', row.names=FALSE)
  },
  error=function(err){
    print(err)
    showNotification("Warning: Properties sheet doesn't exist. Reverting to default values.", duration = notificationDuration, closeButton = TRUE, type = "warning")
  })
  tryCatch({
    dat <- read_excel(file, sheet = 'data')
    pivoted_influent <- pivot_wider(filter(dat, type == 'influent')[, 2:ncol(dat)], names_from = 'compound', values_from = 'concentration')
    pivoted_effluent <- pivot_wider(filter(dat, type == 'effluent')[, 2:ncol(dat)], names_from = 'compound', values_from = 'concentration')
    write.csv(pivoted_influent, 'temp_file/dat_influent.csv', row.names=FALSE)
    write.csv(pivoted_effluent, 'temp_file/dat_effluent.csv', row.names=FALSE)
  },
  error=function(err){
    print(err)
    showNotification("Warning: data sheet doesn't exist. Reverting to default values.", duration = notificationDuration, closeButton = TRUE, type = "warning")
  })
  
  # Attempts to read-in name from Excel file, if it doesn't exist it sets it to the name of the Excel file
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

  # Attempts to read-in fouling data from Excel file, if it doesn't exist it sets it to default values
  tryCatch({
    foulingdata<-read_excel(file, sheet='Fouling Data')
    write.csv(foulingdata, 'temp_file/Foulingdata.csv', row.names=FALSE)
  },
  error=function(err){
    print(err)
    foulingdata <- data.frame(WaterFouling=c('Organic Free'), ChemicalFouling=c('halogenated alkenes'))
    write.csv(foulingdata, "temp_file/Foulingdata.csv", row.names=FALSE)
  })
}

#------------------------------------------------------------------------------#
                                      #column_data
#Column data takes input from within the app and creates a data frame that gives
#the GAC modeling function consistent units
#------------------------------------------------------------------------------#
column_data <- function(input) {
  if (input$veloselect == 'Linear') {
    Fv = pi/4 * ((input$Dv * length_conv[input$DiameterUnits])**2) * input$Vv*velocity_conv[input$VelocityUnits]
  } else {
    Fv = input$Fv * volumetric_conv[input$FlowrateUnits]
  }
  if (input$conc_units == 'ug') {
    mass_mult = 1.0
  } else if (input$conc_units == 'ng') {
    mass_mult = 0.001
  } else {
    mass_mult = 1000. ## mg/L
  }
  if (input$tunits2 == 'days') {
    t_mult = 1440. ## minutes/day
  } else {
    t_mult = 60. ## minutes/hour
  }
  
  columndataframe <- data.frame(
    name = c('carbonID',
             'rad', ## media radius (cm)
             'epor', ## porosity (unitless)
             'psdfr', 
             'rhop', ## carbon density
             'rhof', ## packed bed density, apparent density
             'L',  ## column length (cm)
             'wt', ## media mass (g)
             'flrt', ## flow rate (ml/min)
             'diam', ## column diameter (cm)
             'tortu', ## tortuosity
             'influentID',
             'effluentID',
             'units',
             'time',
             'mass_mul',
             'flow_type',
             'flow_mult',
             't_mult'     
    ),
    value = c('Carbon',
              input$prv*length_conv[input$prunits],
              input$EPORv,
              input$psdfrv,
              input$pdv*density_conv[input$pdunits],
              input$adv*density_conv[input$adunits],
              input$Lv*length_conv[input$LengthUnits],
              input$wv*weight_conv[input$wunits],
              input$Fv*volumetric_conv[input$FlowrateUnits],
              input$Dv*length_conv[input$DiameterUnits],
              input$tortuv,
              'influent',
              'Carbon',
              input$conc_units,
              input$timeunits,
              mass_mult,
              'ml',
              0.001,
              t_mult
    )
  )
  
  return(columndataframe)
}

#------------------------------------------------------------------------------#
                             #effluent_data_processor
#Takes the effluent data and creates a unique name for the chemicals so that
#When the effluent data is plotted with influent and computed data they will
#all be distinguishable. THe gather function then takes the current shape that
#it is in and changes it to a shape that is friendlier with plotly
#------------------------------------------------------------------------------#
effluent_data_processor <- function(effluent) {
  if (nrow(effluent) > 1) {                                      #If effluent data is not empty
    mydata <- effluent
    colnames(mydata) <- paste(colnames(mydata), "effluent", sep = "_")#Distinguish the names from the simulated data
    timevec <- data.frame(hours = c(mydata[, 1]))
    cframe <- data.frame(data.frame(mydata[, 2:ncol(mydata)]))
    colnames(cframe) <- paste(colnames(mydata))[2:length(colnames(mydata))]
    concframe <- gather(cframe)                 #Gather into shape that is easy to convert and plot
    effframe <- cbind(timevec, concframe)
    colnames(effframe) <- c("hours", "name", "conc")
  } else {
    effframe <- data.frame(hours = NA, name = NA, conc = NA)
  }

  return(effframe)
}

#------------------------------------------------------------------------------#
                              #influent_chemical_renamer
#Changes the names of the influent chemical data so that when the data is 
#plotted the computed data, influent data, and effluent data can be 
#distinguished
#------------------------------------------------------------------------------#
influent_chemical_renamer <- function(influent) {
  names <- colnames(influent)
  cindata <- data.frame(influent[, 2:ncol(influent)])
  time <- influent[, 1]
  colnames(cindata) <- paste(names[2:length(names)], "influent", sep = "_")
  alldat <- cbind(time, cindata)

  return(alldat)
}

#------------------------------------------------------------------------------#
                                #fitted_chemical_renamer
#This function renames the chemicals in the fitted chemical data so that they
#Are distinct from the computational, effluent, and influent data on the graph
#This function is currently not being used because the fitted data is not being
#plotted.
#------------------------------------------------------------------------------#
fitted_chemical_renamer <- function(fitted_data) {
  for (chemical in 1:nrow(fitted_data)) {
    fitted_data[chemical, 'name'] <- paste(fitted_data[chemical, 'name'], 'fitted', sep = "_")
  }
  
  return(fitted_data)
}

#------------------------------------------------------------------------------#
                                        #influent_organizer
#influent_organizer takes the influent data frame and changes the shape
#of the data frame into something that is easier to manipulate to and plot
#------------------------------------------------------------------------------#
influent_organizer <- function(influent) {
  cindat_organized <- tidyr::gather(influent[2:ncol(influent)])
  cin_time <- influent[, 1]
  cin_prepped <- cbind(cin_time, cindat_organized)
  colnames(cin_prepped) <- c("hours", "name", "conc")
  
  return(cin_prepped)
}

#------------------------------------------------------------------------------#
                                    #process_output
#If the output data frame is not empty, or in other words if the analysis has 
#been ran, then the output data changes shape into something that is easier
#to manipulate and plot
#------------------------------------------------------------------------------#
process_output <- function(dat, input) {
  if (nrow(dat) > 1) {
    dat2 <- dat %>% pivot_longer(!time, names_to = "name", values_to = "conc")
    colnames(dat2) <- c("hours", "name", "conc")
    totaldat <- dat2
  } else {
    totaldat <- dat
  }

  return(totaldat)
}

#------------------------------------------------------------------------------#
                                  #output_conv
#This function takes the data from the analysis and puts the data into ngl units
#------------------------------------------------------------------------------#
output_conv <- function(dat, input) {
  if (nrow(dat) > 1) {
    dat$conc <- dat$conc * mass_conv[input$conc_units]
  } else {
    dat<-dat
  }
  
  return(dat)
}

#------------------------------------------------------------------------------#
                                #get_bv_in_sec
#This function calculates the bed volume time. This is the time it takes for 
#concentrated water to pass through 1000 beds.
#------------------------------------------------------------------------------#
get_bv_in_sec <- function(input) {
  #get number of seconds per bv
  if (input$veloselect == 'Linear') {
    Vv = input$Vv * velocity_conv[input$VelocityUnits]
  } else {
    Vv = input$Fv * volumetric_conv[input$FlowrateUnits] / (pi / 4 * ((input$Dv * length_conv[input$DiameterUnits]) ** 2))
  }
  
  ## divide converted length by velocity to get BV in seconds
  return(input$Lv * length_conv[input$LengthUnits] / Vv)
}

#------------------------------------------------------------------------------#
                            #create_plotly
#This function creates the plot that is outputted on the output tab
#Frame 1 is the computed data
#Frame 2 is the effluent data
#Frame 3 is the influent data
#------------------------------------------------------------------------------#
create_plotly<-function(frame1, frame2, frame3) {
  #Create a subset of data that 
  computationaldata <- frame1
  effluentdata <- frame2
  influentdata <- frame3
  
  #Using the curated data, plot
  counterionfig <- plot_ly(computationaldata, x = ~ hours, y = ~ conc, type = 'scatter', mode = 'lines', color = ~ name, colors = SteppedSequential5Steps) %>%
                   add_trace(data = effluentdata, x = ~ hours, y = ~ conc, mode = 'markers') %>%
                   add_trace(data = influentdata, x = ~hours, y = ~ conc, mode = 'lines+markers')
  
  return(counterionfig)
}

#------------------------------------------------------------------------------#
                                #cc0_conv_ngl
#This function divides the computed data by the inital influent concentration of 
#that chemical to put the chemical into units of c/c0
#------------------------------------------------------------------------------#  
cc0_conv_ngl <- function(concdata, dataoutput) {
  if (nrow(dataoutput) > 1) {
    df <- concdata
    output <- dataoutput[, colnames(dataoutput)[colnames(dataoutput) != 'time']]
    output_time <- dataoutput['time']
    c0_values <- df[1, 2:ncol(df)]
    output_in_cc0_dat <- mapply('/', output, c0_values)
    output_in_cc0 <- cbind(output_time, output_in_cc0_dat)
    colnames(output_in_cc0) <- colnames(concdata)
  } else {
    output_in_cc0 <- data.frame(time = c(NA), ame = c(NA), conc = c(NA))
  }

  return(output_in_cc0)
}  

#------------------------------------------------------------------------------#
                                  #c_points_cc0
#This function takes the effluent data of a chemical and divides the effluent 
#data by the inital influent concentration of that chemical to put the chemical
#into units of c/c0
#------------------------------------------------------------------------------#  
c_points_cc0 <- function(concdata, effluent) {
  if (nrow(effluent) > 1) {
    df <- concdata
    c0_values <- df[1, 2:ncol(df)]
    effluentdata <- effluent
    just_effluent_conc <- effluentdata[, colnames(effluentdata)[colnames(effluentdata) != 'time']]
    effluent_time <- effluentdata['time']
    effluent_cc0 <- mapply('/', just_effluent_conc, c0_values)
    effluent_cc0_df <- cbind(effluent_time, effluent_cc0)
    colnames(effluent_cc0_df) <- colnames(concdata)
  } else {
    effluent_cc0_df <- data.frame(hours = c(NA), name = c(NA), conc = c(NA))
  }

  return(effluent_cc0_df)
}

#------------------------------------------------------------------------------#
                              #infdat_prep
#This function orders influent data in the format used in the data sheet of
#the Excel file
#------------------------------------------------------------------------------#
infdat_prep <- function(inf_pivoted) {
    inf_pivoted <- cbind("influent", inf_pivoted)
    colnames(inf_pivoted) <- c('type', 'time', 'compound', 'concentration')
    inf_pivoted_ordered <- inf_pivoted[, c('type', 'time', 'concentration', 'compound')]

    return(inf_pivoted_ordered)
}

#------------------------------------------------------------------------------#
                              #effdat_prep
#This function orders effluent data in the format used in the data sheet of
#the Excel file
#------------------------------------------------------------------------------#
effdat_prep <- function(eff_pivoted) {
  eff_pivoted <- cbind("effluent", eff_pivoted)
  colnames(eff_pivoted) <- c('type', 'time', 'compound', 'concentration')
  eff_pivoted_ordered <- eff_pivoted[, c('type', 'time', 'concentration', 'compound')]

  return(eff_pivoted_ordered)
}

#------------------------------------------------------------------------------#
                                  #columnspecs_prep
#This function combines column name, values, and units into one data frame and
#orders it in the format used in the columnSpecs sheet of the Excel file
#------------------------------------------------------------------------------#  
columnspecs_prep <- function(prunits, LengthUnits, wunits, FlowrateUnits, DiameterUnits, prv, EPORv, psdfrv, pdv, adv, Lv, wv, Fv, Dv, tortuv, conc_units, tunits2) {
  data.frame(
    name = c('CarbondID', 'radius', 'porosity', 'psdfr', 'particleDensity', 'apparentDensity', 'length', 'weight', 'flowrate', 'diameter', 'tortuosity', 'influentID', 'effluentID', 'units', 'time'),
    values = c('F400', prv, EPORv, psdfrv, pdv, adv, Lv, wv, Fv, Dv, tortuv, 'influent', 'effluent', conc_units, tunits2),
    units = c(NA, prunits, NA, NA, 'g/ml', 'g/ml', LengthUnits, wunits, FlowrateUnits, DiameterUnits, NA, NA, NA, NA, NA)
  )
}

read_in_files(input, paste0("GAC_config.xlsx"))

#==============================================================================#
#------------------------------------------------------------------------------#
                                  #UI SECTION#
#------------------------------------------------------------------------------#
#==============================================================================#
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
     <h1  class='page-title'>Granular Activated Carbon Model</h1>
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

#------------------------------------------------------------------------------#
                                    #Input Tab#
#------------------------------------------------------------------------------#             
    tabPanel("Input",

#------------------------------------------------------------------------------#
                              #Side Bar on Input Tab#
#------------------------------------------------------------------------------#                      
      sidebarLayout(
        sidebarPanel(
          fileInput("file1", "Choose .xlsx File", accept = ".xlsx"),
          tableOutput("selectedfile"),

          br(),

          h4("Fouling"),  
            
          selectInput("WFouling", "Water Type", list(
            'Organic Free',
            'Rhine',
            'Portage',
            'Karlsruhe',
            'Wausau',
            'Houghton')),
          
          selectInput("CFouling", "Chemical Type", list(
            'halogenated alkenes',
            'halogenated alkanes',
            'halogenated alkanes QSPR',
            'trihalo-methanes',
            'aromatics',
            'nitro compounds',
            'chlorinated hydrocarbon',
            'phenols',
            'PNAs',
            'pesticides',
            'PFAS')),
            
          br(),
          
          sliderInput("nrv", "Radial Collocation Points",3, 18, 7),
          sliderInput("nzv", "Axial Collocation Points", 3, 18, 13),
          
          br(),
          
          actionButton("run_button", "Run Analysis", icon=icon("play")),
          
          br(), br(),
          
          actionButton("Stop", "Stop App", icon=icon("square"), style="color: #000000; background-color: #ff0000; border-color: #e60000")
        ),

#------------------------------------------------------------------------------#
                            #Main Panel on Input Tab#
#------------------------------------------------------------------------------#                      
        mainPanel(
          tabsetPanel(
            tabPanel("Column Parameters",
              
              br(), br(),
              
              fluidRow(
                column(3, HTML(paste0("<h4>","<strong>", "Media Characteristics", "</strong>", "</h4>"))),
                column(3, shinyWidgets::autonumericInput(
                  inputId = "prv",
                  label="Particle Radius",
                  value = 0.0513,
                  decimalPlaces = 4,
                  digitGroupSeparator = ",",
                  decimalCharacter = "."
                )),
                column(3, selectInput("prunits", "Particle Radius Units", c("cm", "m", "mm", "in", "ft")))),
              
              fluidRow(
                column(3,),
                column(3, shinyWidgets::autonumericInput(
                  inputId = "EPORv",
                  label="Bed Porosity",
                  value = "0.641",
                  currencySymbolPlacement = "p",
                  decimalPlaces = 3,
                  digitGroupSeparator = ",",
                  decimalCharacter = "."
                ))),
              
              fluidRow(
                column(3,),
                column(3, shinyWidgets::autonumericInput(
                  inputId = "pdv",
                  label="Particle Density",
                  value = 0.803,
                  currencySymbolPlacement = "p",
                  decimalPlaces = 3,
                  digitGroupSeparator = ",",
                  decimalCharacter = "."
                )),
                column(3, selectInput("pdunits", "Particle Density Units", c("g/ml")))),
              
              fluidRow(
                column(3,),
                column(3, shinyWidgets::autonumericInput(
                  inputId = "adv",
                  label="Apparent Density",
                  value = 0.5,
                  currencySymbolPlacement = "p",
                  decimalPlaces = 3,
                  digitGroupSeparator = ",",
                  decimalCharacter = "."
                )),
                column(3, selectInput("adunits", "App. Density Units", c("g/ml")))),
              
              fluidRow(
                column(3,),
                column(3, shinyWidgets::autonumericInput(
                  inputId = "psdfrv",
                  label="PSDFR",
                  value = 5,
                  currencySymbolPlacement = "p",
                  decimalPlaces = 3,
                  digitGroupSeparator = ",",
                  decimalCharacter = "."
                ))),
              
              hr(),
              
              fluidRow(
                column(3, HTML(paste0("<h4>","<strong>", "Column Specifications", "</strong>", "</h4>"))),                                        
                column(3, shinyWidgets::autonumericInput(
                  inputId = "Lv",
                  label="Length",
                  value = 8,
                  decimalPlaces = 3,
                  digitGroupSeparator = ",",
                  decimalCharacter = "."
                )),
                column(3, selectInput("LengthUnits", "Length Units", c("cm", "ft", "m", "mm", "in")))),

              fluidRow(
                #This radio button toggles between Linear and volumetric flowrate
                column(3, radioButtons("veloselect", "", c("Volumetric", "Linear"))),                                        
                column(3, shinyWidgets::autonumericInput(
                  inputId = "Vv",
                  label="Velocity",
                  value = 0.123,
                  decimalPlaces = 3,
                  digitGroupSeparator = ",",
                  decimalCharacter = "."
                )),
                column(3, selectInput("VelocityUnits", "Velocity Units", c("cm/s", "m/s", "m/min", "m/h", "in/s","ft/s","ft/min", "gpm/ft^2")))),
                
              fluidRow(
                column(3,),
                column(3, shinyWidgets::autonumericInput(
                  inputId = "Dv",
                  label="Diameter",
                  value = 10,
                  decimalPlaces = 3,
                  digitGroupSeparator = ",",
                  decimalCharacter = "."
                )),
                column(3, selectInput("DiameterUnits","Diameter Units",c("cm", "ft","mm", "m", "in")))),

              fluidRow(
                column(3,),
                column(3, shinyWidgets::autonumericInput(
                  inputId = "Fv",
                  label="Flow Rate",
                  value = 500,
                  decimalPlaces = 2,
                  digitGroupSeparator = ",",
                  decimalCharacter = "."
                )),
                column(3, selectInput("FlowrateUnits","Flow Rate Units",c("cm^3/s", "m^3/s", "ft^3/s", "mL/s", "L/min", "mL/min", "gpm", "mgd")))),

              fluidRow(
                column(3,),
                column(3, shinyWidgets::autonumericInput(
                  inputId = "wv",
                  label="Weight",
                  value = 8500,
                  currencySymbolPlacement = "p",
                  decimalPlaces = 3,
                  digitGroupSeparator = ",",
                  decimalCharacter = "."
                )),
                column(3, selectInput("wunits", "Weight Units", c("g", "kg", "lbs", "oz")))),
                                
              fluidRow(
                column(3,),
                column(3, shinyWidgets::autonumericInput(
                  inputId = "tortuv",
                  label="Tortuosity",
                  value = 1,
                  currencySymbolPlacement = "p",
                  decimalPlaces = 3,
                  digitGroupSeparator = ",",
                  decimalCharacter = "."
                ))),
              
              fluidRow(
                column(3,),
                column(3, selectInput("conc_units", "Concentration Units", c("ug", "ng", "mg")))),

              fluidRow(
                column(3,),
                column(3, selectInput("tunits2", "Time Units", c("days", "hours")))),       
            ),

            #Properties Tab
            tabPanel("Compounds",
              h4("Compound List"),
              dataEditUI("edit-1"),
              br(),
              h4("K Data"),
              dataEditUI("edit-2"),
              br(),
              h4("Influent Concentration Data"),
              dataEditUI("edit-3"),
              br(),
              h4("Effluent Concentration Data"),
              dataEditUI("edit-4")      
            ),
          )
        )
      )                
    ),

#------------------------------------------------------------------------------#
                                #Output Tab#
#------------------------------------------------------------------------------#             
    tabPanel("Output",

#------------------------------------------------------------------------------#
                          #Side Bar on Output Tab#
#------------------------------------------------------------------------------#                        
      sidebarLayout(
        sidebarPanel(
          selectInput("OCunits", "Output Concentration Units", c("mg/L", "ug/L", "ng/L", "c/c0")),
          selectInput("timeunits","Output Time Units",c("Days", "Bed Volumes (x1000)", "Hours", "Months", "Years")),
          
          br(),
          
          checkboxInput("computeddata", "Computed Data", TRUE),
          checkboxInput("effluentdata", "Effluent Data", FALSE),
          checkboxInput("influentdata", "Influent Data", FALSE),
          
          br(),
          
          HTML(paste0("<h5>","<strong>", "Effluent Fitting", "</strong>", "</h5>")),
          
          radioButtons("xn", "Options for 1/n increment", choices=c(0.01, 0.025, 0.05), inline=TRUE),
          sliderInput("pm", "Range of K values to test (± %)",0, 50, 30, step=5),
          
          actionButton('fitting', 'Fit Data'),
          
          br(), br(), br(),
          
          downloadButton("save_button", "Save Data"),
          
          actionButton("Stop", "Stop App", icon=icon("square"), style="color: #000000; background-color: #ff0000; border-color: #e60000")                    
        ),

#------------------------------------------------------------------------------#
                      #Main Panel on Output Tab#
#------------------------------------------------------------------------------#                             
        mainPanel(
          shinycssloaders::withSpinner(plotlyOutput("Plot")),#Counterions

          br(),

          textOutput("CounterIonPlot"),    
          plotlyOutput('Plot2'),      
        )
      )
    ),

#------------------------------------------------------------------------------#
                                #Fitted Data Tab#
#------------------------------------------------------------------------------#               
    tabPanel(HTML("Fitted Data</a></li><li><a href='https://github.com/USEPA/Water_Treatment_Models/blob/master/Shiny-GAC/README.md' target='_blank'>Help"), shinycssloaders::withSpinner(uiOutput('FitK')),
      actionButton('Use', 'Use Data'),
      h6('Note: This will replace the K Data in the compounds tab on the Input tab and the modeled output on the Output tab. The model must be run again to view the updated output.')           
    )     
  )
)

#==============================================================================#
#------------------------------------------------------------------------------#
                            #SERVER SECTION#
#------------------------------------------------------------------------------#
#==============================================================================#
server <- function(input, output, session) {
  #Read in file 
  observeEvent(input$file1, {
    file <- input$file1
    read_in_files(input, paste0(file$datapath))
  })
   
  #When a file is uploaded, the session is reloaded. We do this because there does
  #Not seem to be any other way to overwrite the DataEditR tables. See the process_file
  #Notes for further elaboration.
  observeEvent(input$file1, {
    session$reload()
  })
  
  #For some reason when the reticulate package and tidyverse have a strange bug
  #When a person uses the "stop" button within RStudio to stop the app R crashes.
  #This can be avoided by added in-app stop buttons here
  observeEvent(input$Stop, {
    stopApp()
  })
  
  fileuploadedname <- read.csv("temp_file/filename.csv")
  output$selectedfile <- renderTable(fileuploadedname)
  
  #These sheets are used to store data for the analysis
  file_direc <- paste(getwd(), '/temp_file/', sep = '')
  properties<-reactiveVal(read.csv(paste(file_direc, "Properties.csv", sep = '')))
  columnSpecs<-reactiveVal(read.csv(paste(file_direc, "columnSpecs.csv", sep = '')))
  Kdata<-reactiveVal(read.csv(paste(file_direc, "Kdata.csv", sep = '')))
  
  #------------------------------------------------------------------------------#
                      #Volumetric vs Linear Velocities#
  #Here we look at the file that has been uploaded and see if the user has
  #uploaded data with linear or volumetric flow rates. Whatever they do not 
  #use will be grayed out
  #------------------------------------------------------------------------------#
  test_df <- data.frame(C = c('v', 'flowrate'))
  flags <- reactive({test_df$C %in% columnSpecs()$name})

  observe({
    if (flags()[1]) {
      updateNumericInput(session, "Vv", value = filter(columnSpecs(), name == 'v')$value)
      updateSelectInput(session, "VelocityUnits", choices = unique(c(filter(columnSpecs(), name == 'v')$units, velocityvector)))
      updateRadioButtons(session, "veloselect", selected = "Linear")
    } else if (flags()[2]) {                         
      updateNumericInput(session, "Fv", value = filter(columnSpecs(), name == 'flowrate')$value)
      updateSelectInput(session, "FlowrateUnits", choices = unique(c(filter(columnSpecs(), name == 'flowrate')$units, flowratevector)))
      updateRadioButtons(session, "veloselect", selected = "Volumetric")
    }
  })

  observe({
    toggleState("Vv", condition = input$veloselect != "Volumetric")
    toggleState("VelocityUnits", condition = input$veloselect != "Volumetric") # Velocity units are grayed out if velocity is not being used
    toggleState("Fv", condition = input$veloselect != "Linear")
    toggleState("FlowrateUnits", condition = input$veloselect != "Linear") # Flowrate units are grayed out if flowrate is not being used
  })
  
  #------------------------------------------------------------------------------#
  #Gathering Uploaded Values#
  #------------------------------------------------------------------------------#  
  #Numeric Values
  CarbonID <- reactive({filter(columnSpecs(), name == "carbonID")$value})
  radius <- reactive({filter(columnSpecs(), name == "radius")$value})
  porosity <- reactive({filter(columnSpecs(), name == 'porosity')$value})
  psdfr <- reactive({filter(columnSpecs(), name == 'psdfr')$value})
  particleDensity <- reactive({filter(columnSpecs(), name == 'particleDensity')$value})
  apparentDensity <- reactive({filter(columnSpecs(), name == 'apparentDensity')$value})
  length <- reactive({filter(columnSpecs(), name == 'length')$value})
  weight <- reactive({filter(columnSpecs(), name == 'weight')$value})
  flowrate <- reactive({filter(columnSpecs(), name == 'flowrate')$value})
  diameter <- reactive({filter(columnSpecs(), name == 'diameter')$value})
  tortuosity <- reactive({filter(columnSpecs(), name == 'tortuosity')$value})
  influentID <- reactive({filter(columnSpecs(), name == 'influentID')$value})
  effluentID <- reactive({filter(columnSpecs(), name == 'effluentID')$value})
  time <- reactive({filter(columnSpecs(), name == 'time')$value})
  nrv <- reactive(7)
  nzv <- reactive(12)
  
  #String Values
  prvec <- reactive({
    prv <- c(filter(columnSpecs(), name == 'radius')$units, prvector)

    return(unique(prv))
  })

  pdvec <- reactive({
    pdv <- c(filter(columnSpecs(), name == 'particleDensity')$units, pdvector)

    return(unique(pdv))
  })

  advec <- reactive({
    adv <- c(filter(columnSpecs(), name == 'apparentDensity')$units, advector)

    return(unique(adv))
  })
  
  lengthvec <- reactive({
    lenvec <- c(filter(columnSpecs(), name == 'length')$units, lengthvector)

    return(unique(lenvec))
  })
  
  weightvec <- reactive({
    wvec <- c(filter(columnSpecs(), name == 'weight')$units, weightvector)

    return(unique(wvec))
  })

  concvec <- reactive({
    concv <- c(filter(columnSpecs(), name == 'units')$value, concentrationvector)

    return(unique(concv))
  })
  
  timevec <- reactive({
    timevec <- c(filter(columnSpecs(), name == 'time')$value, timevector)

    return(unique(timevec))
  })
  
  flowvec <- reactive({
    flowv <- c(filter(columnSpecs(), name == 'flowrate')$units, flowratevector)

    return(unique(flowv))
  })
  
  diamvec <- reactive({
    diamv <- c(filter(columnSpecs(), name == 'diameter')$units, diametervector)

    return(unique(diamv))
  })

  # Updates water type select input choiecs
  wfoulingvec <- reactive({
    wfoulingv <- c(select(read.csv("temp_file/Foulingdata.csv"), WaterFouling), wfoulingvector)

    return(unique(wfoulingv))
  })

  # Updates chemical type select input choices
  cfoulingvec <- reactive({
    cfoulingv <- c(select(read.csv("temp_file/Foulingdata.csv"), ChemicalFouling), cfoulingvector)

    return(unique(cfoulingv))
  })
  
  #------------------------------------------------------------------------------#
          #Updating default values with the values that were uploaded#
  #------------------------------------------------------------------------------#    
  observe({
    updateNumericInput(session, "prv", value = format(radius(), digits = 4, scientific = FALSE))
    updateNumericInput(session, "EPORv", value = porosity())
    updateNumericInput(session, "pdv", value = format(particleDensity(), digits = 4, scientific = FALSE))
    updateNumericInput(session, "adv", value = format(apparentDensity(), digits = 4, scientific = FALSE))
    updateNumericInput(session, "psdfrv", value = psdfr())
    updateNumericInput(session, "Lv", value = length())
    updateNumericInput(session, "wv", value = format(weight(), digits = 4, scientific = FALSE))
    updateNumericInput(session, "Dv", value = diameter())
    updateNumericInput(session, "tortuv", value = tortuosity())
    updateNumericInput(session, "nrv", value = nrv())
    updateNumericInput(session, "nzv", value = nzv())
    updateSelectInput(session, "prunits", choices = prvec())
    updateSelectInput(session, "pdunits", choices = pdvec())
    updateSelectInput(session, "adunits", choices = advec())
    updateSelectInput(session, "LengthUnits", choices = lengthvec())
    updateSelectInput(session, "conc_units", choices = concvec())
    updateSelectInput(session, "tunits2", choices = timevec())
    updateSelectInput(session, "wunits", choices = weightvec())
    updateSelectInput(session, "FlowrateUnits", choices = flowvec())
    updateSelectInput(session, "DiameterUnits", choices = diamvec())
    updateSelectInput(session, "WFouling", choices = wfoulingvec(), selected = select(read.csv("temp_file/Foulingdata.csv"), WaterFouling)) # Updates water type select input
    updateSelectInput(session, "CFouling", choices = cfoulingvec(), selected = select(read.csv("temp_file/Foulingdata.csv"), ChemicalFouling)) # Updates chemical type select input
  })
  
  #------------------------------------------------------------------------------#
                            #Dynamic Data Frames
  #Here we have data frames that act like excel tables. These tables can be
  #manipulated to change values and add/remove rows and columns. These data frames
  #are populated with default values from the config file and get overwritten
  #when a file is uploaded
  #------------------------------------------------------------------------------#  
  #data frame of chemicals and their properties
  iondat <- dataEditServer("edit-1",  data = paste(file_direc, 'Properties.csv', sep = ''))
  dataOutputServer("output-1", data = iondat)
  
  #data frame of k data for each chemical
  kdat <- dataEditServer("edit-2", data = paste(file_direc, 'Kdata.csv', sep = ''))
  dataOutputServer("output-2", data = kdat) 
  
  #influent data for each chemical
  infdat <- dataEditServer("edit-3", data = paste(file_direc, 'dat_influent.csv', sep = ''))
  dataOutputServer("output-3", data = infdat)
  
  #effluent data for each chemical, this is optional
  effdat <- dataEditServer("edit-4", data = paste(file_direc, 'dat_effluent.csv', sep = ''))
  dataOutputServer("output-4", data = effdat)
  
  ##Column_data_converted = column_info
  column_data_converted <- reactive({column_data(input)})

  ##chem_data = properties
  chem_data <- reactive({properties()})
  
  #------------------------------------------------------------------------------#
                        #Running the Analysis/PSDM function
  #Here is the function that runs the analysis and does the heaviest lifting
  #in the app. It takes the arguments: columndata, chem_data, kdata, infdat, 
  #effdat, nr, nz, water_type, and chem_type
  #------------------------------------------------------------------------------#  
  out <- reactiveVal(data.frame(Chemicals = c(0, 0), time = c(0, 0)))
  
  # Detects errors in input values
  error_handling <- eventReactive(input$run_button, {
    errorflag <- 0
    coldensity <- (input$wv * weight_conv[input$wunits]) / (pi * (input$Dv*length_conv[input$DiameterUnits] / 2) ^ 2 * (input$Lv * length_conv[input$LengthUnits]))
    appdensity <- input$adv * density_conv[input$adunits]

    if (coldensity > appdensity) {
      errorflag <- 1
      shinyalert("Error", "Apparent Density value is too low.", type = "error")
    } 
    if (!(input$WFouling %in% wfoulingvector)) {
      errorflag <- 1
      shinyalert("Error", "Water type is not accepted. Please select one from the list.", type = "error")
    } 
    if (!(input$CFouling %in% cfoulingvector)) {
      errorflag <- 1
      shinyalert("Error", "Chemical type is not accepted. Please select one from the list.", type = "error")
    }

    return(errorflag)
  })
  
  observeEvent(input$run_button, {
    if (error_handling() != 1) {
      tryCatch({
        showNotification("Starting model run.", duration = notificationDuration, closeButton = TRUE, type = "message") # Notifies the user that the model is being run
        out(run_PSDM(column_data_converted(), chem_data(), kdat(), infdat(), effdat(), nrv(), nzv(), input$WFouling, input$CFouling))
        updateTabsetPanel(session, "inTabset", selected = "Output") # Switches to Output tab when run button is pressed
      },
      error=function(err){
        shinyalert("Error", "An error is preventing the model from running, please consult the README for more information.", type = "error")
      })
    }
  })
  
  computed_data_prep <- reactive({process_output(out())})
  computed_data <- reactive({output_conv(computed_data_prep(), input)})
  
  #------------------------------------------------------------------------------#
              #Running the fit to the Analysis/PSDM function
  #This function returns computed data fitted to the effluent data and kdata that
  #is fitted to the effluent data. Currently, on the kdata part is being used. 
  #------------------------------------------------------------------------------#  
  out_fit <- reactiveVal(data.frame(hours = c(NA), name = c(NA), conc = c(NA))) #Stores PSDM function
  output_fit <- reactiveVal(data.frame(hours = c(NA), name = c(NA), conc = c(NA))) #Stores first value of PSDM function
  kdata_fit <- reactiveVal(data.frame(Chemical = c(0, 0, 0, 0, 0))) #Stores second value of PSDM function
  kdata_fit_save <- reactiveVal(data.frame(Chemical = c(0, 0, 0, 0, 0))) #Used to Export Save File
  kdataframe <- data.frame(name = c('K', '1/n', 'q', 'brk', 'AveC')) #Used label values
  colnames(kdataframe) <- c('...1') #For some reason if this isn't here it crashes
  output$FitK <- renderTable(kdataframe)
  
  observeEvent(input$fitting, {
    if (nrow(out()) > 2) {
      showNotification("This might take several minutes.", duration = notificationDuration, closeButton = TRUE, type = "message")
      showNotification("Fitting data.", duration = notificationDuration, closeButton = TRUE, type = "message")
      out_fit(run_PSDM_fitter(column_data_converted(), chem_data(), kdat(), infdat(), effdat(), nrv(), nzv(), input$WFouling, input$CFouling, input$pm, input$xn))
      output_fit(out_fit()[[1]])
      kdata_fit(out_fit()[[2]])
      kdataframe <- cbind(kdataframe, kdata_fit())
      kdata_fit_save(kdataframe)
      output$FitK <- renderTable({kdataframe})
      updateTabsetPanel(session, "inTabset", selected = HTML("Fitted Data</a></li><li><a href='https://github.com/USEPA/Water_Treatment_Models/blob/master/Shiny-GAC/README.md' target='_blank'>Help")) # Switches to Fitted Data tab when fit button is pressed
    } else {
      output$FitK <- renderText({'No fitting data available'})
    }
  })
  
  #Changing the shape of fitted kdata to be able to run in the analysis again
  kdat_fitted <- reactive({
    df <- kdata_fit()
    df <- cbind(kdataframe, df)
    rownames(df) <- 1:nrow(df)

    return(df)
  })
  
  #Rerunning the analysis with fitted kdata
  observeEvent(input$Use, {
    write.csv(kdat(), paste(file_direc, 'Kdata2.csv', sep = ''), row.names = FALSE)
    write.csv(kdat_fitted(), paste(file_direc, 'Kdata.csv', sep = ''), row.names = FALSE)
    kdat<- dataEditServer("edit-2", data = paste(file_direc, 'Kdata.csv', sep = ''))
    dataOutputServer("output-2", data = kdat)
    if ((read.csv("temp_file/Kdata.csv")[1, 2]) == 0) {
      shinyalert("Error", "K Data is empty.", type = "error")
    } else {
      showNotification("Updating K Data.", duration = notificationDuration, closeButton = TRUE, type = "message")
      out(out_fit()[[1]])
      write.csv(data.frame(WaterFouling = c(input$WFouling), ChemicalFouling = c(input$CFouling)), paste(file_direc, 'Foulingdata.csv', sep = ''), row.names = FALSE) # Saves water fouling data to an Excel file to be read-in after the session is reloaded
      session$reload() # When K Data is used, the session is reloaded in order to update the DataEditR tables
    }
  })
  
  #Putting Data into correct shapes
  effdat_plot <- reactive({
    output_conv(effluent_data_processor(effdat()), input)
  })
  
  influent_plot <- reactive({
    dat <- influent_chemical_renamer(infdat())
    output_conv(influent_organizer(dat), input)
  })
  
  #------------------------------------------------------------------------------#
                                #cc0 conversions#
  #------------------------------------------------------------------------------# 
  computed_data_cc0 <- reactiveVal(data.frame(time = c(NA), name = c(NA), conc = c(NA)))
  cc0data_ngl <- eventReactive(input$run_button, {cc0_conv_ngl(infdat(), out())})
  observe({computed_data_cc0(process_output(cc0data_ngl()))})
  
  effluent_data_cc0 <- reactiveVal(data.frame(time = c(NA), name = c(NA), conc = c(NA)))
  effluentcc0data <- reactive({c_points_cc0(infdat(), effdat())})
  observe({effluent_data_cc0(process_output(effluentcc0data()))})
  
  influent_data_cc0 <- reactiveVal(data.frame(time = c(NA), name = c(NA), conc = c(NA)))
  influentcc0data <- reactive({c_points_cc0(infdat(), infdat())})
  observe({influent_data_cc0(process_output(influentcc0data()))})
  
  #------------------------------------------------------------------------------#
                          #Converting the output data
  #------------------------------------------------------------------------------#   
  outputeffluent<-reactiveValues()
  outputinfluent<-reactiveValues()
  outputchemicals<-reactiveValues()
  outputfit<-reactiveValues()
  
  observe({
    ## convert time units for graphing
    # calculating kBV
    if (input$timeunits == "Bed Volumes (x1000)") {
      bv_conv <- get_bv_in_sec(input)
      outputchemicals$hours <- computed_data()$hours / (bv_conv / hour2sec) / 1e3
      outputeffluent$hours <- effdat_plot()$hours / (bv_conv / hour2sec) / 1e3
      outputinfluent$hours <- influent_plot()$hours  / (bv_conv / hour2sec) / 1e3  
    } else {
      outputchemicals$hours <- computed_data()$hours * (time_conv[input$timeunits])
      outputeffluent$hours <- effdat_plot()$hours * (time_conv[input$timeunits])
      outputinfluent$hours <- influent_plot()$hours * (time_conv[input$timeunits])
    }
  })
  
  observe({
    ### convert y-axis/mass units for graphing
    if (input$OCunits == "c/c0") {
      ## just replicates the returned data
      outputchemicals$conc <- computed_data_cc0()$conc
      outputeffluent$conc <- effluent_data_cc0()$conc
      outputinfluent$conc <- influent_data_cc0()$conc
    } else {
      outputchemicals$conc <- computed_data()$conc / mass_conv[input$OCunits]
      outputeffluent$conc <- effdat_plot()$conc / mass_conv[input$OCunits]
      outputinfluent$conc <- influent_plot()$conc / mass_conv[input$OCunits]
    }
  })
  
  #------------------------------------------------------------------------------#
                #Putting together the data frames to plot
  #Now we put together all of output data that has been converted into plots
  #This is the easiest way to handle plotting the data with plotly in the next
  #section.
  #------------------------------------------------------------------------------#
  computational_processed <- reactive({
    if (input$computeddata == TRUE) {
      plotdata <- computed_data()
      plotdata$conc <- outputchemicals$conc
      plotdata$hours <- outputchemicals$hours
    } else {
      plotdata <- data.frame(hours = c(NA), name = c(NA), conc = c(NA))
    }

    return(plotdata)
  })
  
  effluent_processed <- reactive({
    if (input$effluentdata == TRUE) {
      plot_data3 <- effdat_plot()
      plot_data3$conc <- outputeffluent$conc
      plot_data3$hours <- outputeffluent$hours
    } else {
      plot_data3 <- data.frame(hours = c(NA), name = c(NA), conc = c(NA))
    }

    return(plot_data3)
  })
  
  influent_processed <- reactive({
    if (input$influentdata == TRUE) {
      plot_data4 <- influent_plot()
      plot_data4$conc <- outputinfluent$conc
      plot_data4$hours <- outputinfluent$hours
    } else {
      plot_data4 <- data.frame(hours=c(NA), name = c(NA), conc = c(NA))
    }

    return(plot_data4)
  })
  
  #------------------------------------------------------------------------------#
                        #Saving Output Data To .xlsx
  #This section is where everything comes together, the computational, effluent,
  #and influent data are plotted
  #------------------------------------------------------------------------------#
  fig <- reactive({create_plotly(computational_processed(), effluent_processed(), influent_processed())})
  counterionfigure <- reactive({fig() %>% layout(title = "Concentration over Time", showlegend = TRUE,
                                                 legend = list(orientation = 'h', y=1), hovermode = 'x unified',
                                                 xaxis = list(title=input$timeunits, gridcolor = 'ffff'),
                                                 yaxis = list(title=paste0("Concentration (", input$OCunits, ")"), showexponent = 'all', exponentformat = 'e', gridcolor = 'ffff'))
  })
  
  output$Plot <- renderPlotly(counterionfigure())
 
  #------------------------------------------------------------------------------#
                          #Saving Output Data To .xlsx
  #This sections takes, the Properties, Kdata, columnSpecs, data (influent and
  #effluent data), Model Results, Fit Data, and Fouling Data data frames and
  #saves them to an excel file with the exact values that were inputted into
  #The model.
  #------------------------------------------------------------------------------#
  infdatsave <- reactive({
    infdat_prep(infdat() %>% pivot_longer(!time, names_to = 'compound', values_to = 'concentration'))
  })

  effdatsave <- reactive({
    effdat_prep(effdat() %>% pivot_longer(!time, names_to = 'compound', values_to = 'concentration'))
  })
  
  columnspecssave <- reactive({
    df <- columnspecs_prep(input$prunits, input$LengthUnits, input$wunits, input$FlowrateUnits, input$DiameterUnits, input$prv, input$EPORv, input$psdfrv, input$pdv, input$adv, input$Lv, input$wv, input$Fv, input$Dv, input$tortuv, input$conc_units, input$tunits2)
    colnames(df) <- c('name', 'value', 'units')

    return(df)
  })
  
  #The influent and effluent data were split up to treat them independently earlier
  #Now they are being brought back together 
  outputconcentrations <- reactive({rbind(infdatsave(), effdatsave())})
  
  #Fouling data is saved in the excel file just in case
  foulingdata <- reactive({data.frame(WaterFouling = c(input$WFouling), ChemicalFouling = c(input$CFouling))})
  
  #Fouling data is saved in the excel file just in case  
  output$save_button <- downloadHandler(
    filename = function() {
      paste("data-", Sys.Date(), ".xlsx", sep="")
    },
    content=function(file) {
      sheets <- list("Properties" = chem_data(),
                     "Kdata" = kdat(),
                     "columnSpecs" = columnspecssave(),
                     "data" = outputconcentrations(),
                     "Model Results" = computational_processed(),
                     'Fit Data' = kdata_fit(),
                     'Fouling Data' = foulingdata()
      )

      write_xlsx(sheets, file)
    }
  )
}

shinyApp(ui, server)