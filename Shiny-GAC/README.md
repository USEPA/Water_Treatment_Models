# Granular Activated Carbon Tool

The Granular Activated Carbon Modeling Tool is used to model a __________________ in a drinking water treatment plant. This model relies on ______________________ and predicts the breakthrough behavior for unit operation design. To read a further in-depth analysis of the theory behind this model please reference 

1. [Requirements](#requirements)
2. [Excel Formatting](#excel-based-input-file)
3. [Quick Start](#quick-start)
4. [Appendix](#appendix)
5. [Notes to the User](#notes-to-the-user)
6. [Development Team](#development-team)

## Requirements 
1. R/R Studio (At least version 2022.7)
2. Python (At least version 3.9.7)
3. Excel (recommended)
4. Files: GACapp.R, config.xlsx
5. Optional: Example_TCE.xlsx

## Excel-based Input File
The input for the Excel-based input file must be formatted like the one shown in the figure below if the user wants to import data. The Shiny GAC app looks for sheetnames of "Properties", "Kdata", "columnSpecs", "data", "data_variable" (optional) and "data_optimize" (optional). If one or more of those sheets are not found then the app cannot be run using that input file. The app is loaded with default data if the user does not want to use an Excel-based file, and additional changes can be made within the GUI. There is a fourth and fifth optional sheet "data_variable", which represents varying influent and effluent concentration. There is also "data_optimize" which ___________. These pages do not need to be present when the file is ran and is not required to exist in the excel sheet.

<figure>
    <img src="DocumentPics/excelsheet.PNG"
         alt="Excel Input">
</figure>

<figure>
     <img src="DocumentPics/excelsheet2.PNG"
         alt="Excel Input">
     <figcaption>The Excel file consists of three sheets: parameters of the system, the list of ions that the user is interested in along with their properties, and the list of concentrations for the ions at a given time. Each tab is broken down in detail in the features section of this document.
    * Dp and Dp_units are provided in example input files. These are required for PSDM modeling but will be ignored for HSDM modeling.
    </figcaption>
</figure>

&nbsp;

## Set Up

In order for the tool to work the user must point their R Studio to a Python Interpreter

1. Open RStudio

<figure>
    <img src="DocumentPics/selecttool.PNG"
         alt="Excel Input">
</figure>

2. Select the "Tools" tab at the top of the page and then select "Global Variables"

<figure>
    <img src="DocumentPics/selectglobal.PNG"
         alt="Excel Input">
</figure>

3. Go to the Python tab in Global Variables, then select 'Browse' and select the file where Python is installed locally.

<figure>
    <img src="DocumentPics/pythonselect.PNG"
         alt="Excel Input">
</figure>

4. Click "Ok" then "Apply"


## Quick Start

1. In RStudio, click the "Run App" button in the top right corner of the window that contains the code

![Start](DocumentPics/start.PNG)

2. Application opens to Input & Column Parameter view. User can select "Browse" below "Choose .xlsx File" to import preconfigured input file in the upper left, as shown. Or, a user can begin editing input values following steps 4 and 5 below.

![Import File](DocumentPics/Slide1.PNG)


3. (Optional) Change the parameters to match the specifications of your Ion Exchange apparatus

![Inputs](DocumentPics/Slide2.PNG)

4. (Optional) In the ions tab, add chemicals and concentration points to match your interest. These can be added or edited by right clicking the data table.

![Adjust](DocumentPics/IonsTab.PNG)
![IonEdit](DocumentPics/Slide4.PNG)

5.	Click the Run Analysis button that’s at the bottom of the same side panel as the file import. 3 total chemicals with 2 concentration points takes about 10-30 seconds.

![Run](DocumentPics/Slide5.PNG)

6.	Switch to the Output tab (There should be a loading “spinner” to let you know it’s running)

![waiting](DocumentPics/Slide6.PNG)

7.	Your graph will appear. You can export the data as an xlsx file along with the conditions you input.

![plot](DocumentPics/Slide7.PNG)

8. The user can use the check boxes on the left to toggle the influent data and the effluent data (if available). Note that the user can toggle on and off individual traces on the graph by clicking on the desired data on the legend (data will be grayed out if it isn't displayed).

![plottraces](DocumentPics/Slide8.PNG)




## Appendix

### Column Parameters

The Column Parameters tab (Input>Column Parameters) is used to describe the resin characteristics and column specifications. Some of these parameters, like tortuosity, can be nontrivial to measure so appropriate references are provided where a user can find additional information.

The parameters tab is used to describe the physical constraints of the resin characteristics and column specifications. Some of these measurements, like Resin Capacity, can be nontrivial to measure so we have tried to supply a source where the user can find the information if they do not have it already.

| Input                        | Description | Units |
|---                           |---          |---    |
|Carbon ID                     |             |N/A    |
|radius                        |Bead radius is the measurement of the distance of the bead resin from the center to the surface. | cm, m, mm, in, ft |
|Bead Porosity | The bead porosity is the measure of the bead volume occupied by a solvent, usually water. The factor is between 0 and 1, where 0 represents a bead absent of a solvent and 1 is a bead where all the available space is filled with a solvent. A well packed bead will typically have and EPOR of 0.2. | N/A|
|psdfr                          |            |N/A    |
|particle density               | Mass per unit volume of bead particle. | g/ml |
|apparent density               | Mass per unit volume of bead particle measured within a medium. | g/ml|
|length                         |The depth of the media in packed column. Some vessels may only be filled partially, so this number may be shorter than the height of the contractor.| m, cm, mm, in, ft|
|weight                         | Weight of carbon. | lb, kg, g |
|flowrate                       | The average flow rate through the column. HSDM only considers and average or steady-state condition, not variable flow. | gpm, cm^3/s, m^3/s, ft^3/s, mL/s, L/min, mL/min, mgd |
|diameter                     | The diameter of a cylindrical column.  | cm, mm, m, ft, in |
|tortuosity                   | Parameter of flow between curve and length. | N/A |
|units                        | Influent and Effluent Concentrations | ug, ng, mg |
|time                         | The units for time in the corresponding "data" sheet in the Excel-based files or "Concentration Points" table under Input>Ion's tab. |




### Ions Tab

The Ions tab contains information about the ions to be simulated. The Ions tab should contain any anionic species found in the system, as HSDM always models a simultaneous competitive exchange process (counterions like sulfate, nitrate, chloride, and bicarbonate must also be modeled in simulations). Ion characteristics (molecular weight, selectivity, valence, mass transfer coefficients) are stored in the "ions" tab in Excel-based input files and listed under "Ion List" within HSDM under the Input>Ions tab. information on concentrations within the system over time are stored in the "Cin" tab in Excel-based input files or "Concentration Points" within HSDM under the Input>Ions tab.

**Note:** Order of rows in the ions table should match order of columns in the Cin table.


### Ion List
|  Input        	                |Column Name   | Description                                                                      |
|---            				    |---        |---                                                                               |
|Molecular Weight               |MW         |Molecular weight of ionic species. Gram per mol.                                          |           |
|Molar Volume                   |MolarVol       |Volume occupied by one mole of ion substance. |                     
|Boiling Point                           | BP          |Temperature in Celsius that the given chemical boils. |
|Density      |        |Mass per unit volume of ion.     |
|Solubility|    | Amount of substance that will disolve in water. |
|Vapor Pressure  |         |Pressure of vapor form of substance. |      


### K Data

|Name                   |Variable      | Description            |
|---                    |---           |---                     |
|K                      | K            |                        |
|                       | 1/n          |                        |
|                       | q            |                        |
|                       | brk          |                        |
|                       | AveC         |                        |

### Concentration Points / Influent and Effluent Data

The Concentration Points table (stored in Excel-based files in the "data" sheet) represents time series of concentrations in the system over time (conc_units for an ion defines the concentration units for each column). NOTE: The duration of a simulation is specified by the largest time/last row in the Time column. Times in this table should be specified in ascending order, and a minimum of two (2) rows are needed to run a simulation (time 0, and run duration). If concentrations are the same in all rows, the simulation has a constant concentration, but variable concentrations can be modeled (where linear interpolation is used between points). Time units for the times in this table is specified on the "Input>Column Parameters" tab as time (day or hr). This data can be plotted as well by simply clicking the "Influent Data" button on the 'Output' page.

The effluent data represents the results from a previous experiment or results from a previous model. This data does not need to be present for the model to run. If the user would like to see this data on the plot, the data must be presentwith the influent in the 'data' sheet. When the analysis has been ran, the user can toggle on the "effluent data" button on the left hand side to see the data ploted. The units that the effluent data is in is represented by the time units as well i nthe column parameters. Whatever unit corresponds to the chemical in the 'concentration units' selection will be converted to the 'Output Units' on the output page. Moreover, if the user wants to plot the effluent data then the chemical the user wishes to plot must exist in both the 'ions' page, the 'cin page', and the 'effluent page'.

### Output
The Output tab provides a graphical output of results after a simulation is completed. The units for time and concentration can be adjusted via the dropdown menus in the left-hand column.

![output](DocumentPics/output.png)

The graphs will dynamically update to reflect the current selection. If c/c0 is selected, the output is scaled relative to c0 (the initial concentrations), which is the first row of the "Concentration Points" table. "Bed volumes x1000" presents time in thousands of bed volumes treated, where bed volume is the empty bed volume of the media in the system (relating to bed length if velocity is specified, or bed length and diameter if flow rate is specified).

The data can be exported by clicking the save button. This saves the data points to an excel file where the chemicals inputted into the analysis are the header and the columns are the corresponding concentration points with the first column being the corresponding time.

To export the graph, the user can hover over the graph with their cursor which will display the graph settings in the top right. Clicking the camera icon will open up a file explorer where the user can save the graphic.

![savegraph](DocumentPics/saveplot.png)

## Notes to the User

#### Selection of Collocation Points (nr and nz)
The parameters nr and nz control the size of the grid used to numerically solve the underlying differential equations during the simulation. Increasing nr and nz may increase the accuracy of simulations but doing so also makes them take longer to run. No analytical expression has been found for determining optimal grid dimensions for this class of problems, so selecting nr and nz may take some experimentation. Generally, the sharper the ion exchange zone is relative to the column length, the higher nz will need to be and the sharper the diffusion gradient in the resin beads becomes, the higher nr will need to be. In practice, it is rare for nr to be the controlling parameter for grid size, with nr=7 being accurate enough for most cases without unduly increasing computational cost. The parameter nz is more likely to need attention. Setting nz too low will often produce erroneous oscillations in the breakthrough curves. The illustration below shows simulations with (a) and without (b) these erroneous oscillations.

![nrnz](DocumentPics/Picture1.png)

If increasing nz does not smooth out the erroneous oscillations, there may be other problems with the simulation. In this case, the user is advised to double check the inputs for errors. If there are no errors in the inputs, it is possible the ion exchange zones in the requested simulations are simply too sharp for this numerical approximation to handle. Faced with this problem, the user may wish to consider reducing the empty bed contact time of the simulation or seek out an alternate method of solution such as an equilibrium-based column model.

#### Specifying Column Size and Flow Rate

The underlying model equations in this code use column length (bed depth), L, to define filter size and superficial (linear) flow velocity, v, to define flow rate of simulated systems. If both parameters are readily available to the user, they can be input directly selecting the “linear” radio button on the left side of the column specification section. In practice, flow in adsorption systems is often specified as a hydraulic loading rate (or “surface loading rate”) given in units of volumetric flow rate divided by area (for instance, with units of gpm/ft2 or gallons per minute per square foot). This specification is ultimately equivalent to specifying a superficial flow velocity and can also be entered directly for v provided appropriate units are selected from the corresponding drop-down menu. 

Occasionally design specification may include bed dimensions and empty bed contact time (EBCT) but omit flow information. In this case, the user can obtain a superficial flow velocity from the following formula:

![eq3](DocumentPics/eq3.png)

The column size and flow rate may also be defined in terms of L, bed diameter (d), and volumetric flow rate (flrt) by selecting the “volumetric” radio button.
A note on selection of flow convention: The entry field for the two conventions are independent. The values shown in disabled fields (gray backgrounds) are not updated to correspond to values entered using the other convention. Thus, switching between the radio buttons usually results in switching between two different systems.

**Note:** If flow rate is provided, diameter is required


## References
ACS EST Water 2023, 3, 2, 576–587
Publication Date:January 19, 2023
https://doi.org/10.1021/acsestwater.2c00572

Environ. Sci. Technol. 2021, 55, 8, 5001–5011
Publication Date:March 22, 2021
https://doi.org/10.1021/acs.est.1c00769

## Development Team
David Colantonio

Levi Haupert

Jonathan Burkhardt
