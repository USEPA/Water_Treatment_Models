# Granular Activated Carbon Tool

The Shiny-GAC Modeling Tool is used to model a granular activated carbon (GAC) unit operation in a drinking water treatment plant. Shiny-GAC includes a pore and surface diffusion model (PSDM) to predict the breakthrough behavior for unit operation design. To read a further in-depth analysis of the theory behind this model please reference the [AdDesignS User Manual](https://github.com/USEPA/Environmental-Technologies-Design-Option-Tool/blob/master/VS2019_AdDesignS/ADS_VBnet_MV/help/ads.doc).

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
4. Files: GACapp.R, GAC_config.xlsx
5. Optional: Example_TCE.xlsx

## Excel-based Input File
The input for the Excel-based input file must be formatted like the one shown in the figure below if the user wants to import data. The Shiny-GAC App looks for sheetnames of "Properties", "Kdata", "columnSpecs", and "data". If one or more of those sheets are not found then the app cannot be run using that input file. The Shiny-GAC App also looks for a sheetname of "Fouling Data" and if it is not found it sets the water type to organic free and the chemical type to halogenated alkenes. The App is loaded with default data (found in 'GAC_config.xlsx') if the user does not want to use an Excel-based file, which can be modified as needed within the GUI.

Find out where influent concentration gets converted in IEX

<figure>
    <img src="DocumentPics/Properties.png"
         alt="Excel Input">
</figure>

<figure>
     <img src="DocumentPics/Kdata.png"
         alt="Excel Input">
</figure>

<figure>
    <img src="DocumentPics/columnSpecs.png"
         alt="Excel Input">
</figure>

<figure>
    <img src="DocumentPics/data.png"
         alt="Excel Input">
    <figcaption>The Excel file consists of four sheets: the first is the properties of each chemical the user is interested in, the second is the 'k data' of each of those ions. The number of ions must be the same for both properties and k data. Third is the column specifications of the GAC apparatus. The fourth is a list of the influent and effluent points of each of the chemicals listed in properties and k data tab. Each tab is broken down in detail in the features section of this document.
    </figcaption>
</figure>

&nbsp;

## Set Up

In order for the tool to work the user must point their R Studio to a Python Interpreter. If you are using the web version you can jump to Quick Start step 2.

1. Open RStudio

<figure>
    <img src="DocumentPics/selecttool.PNG"
         alt="Excel Input">
</figure>

2. Select the "Tools" tab at the top of the page and then select "Global Options"

<figure>
    <img src="DocumentPics/selectglobal.PNG"
         alt="Excel Input">
</figure>

3. Go to the Python tab in Global Variables, then select 'Browse' and select the file where Python is installed locally.

<figure>
    <img src="DocumentPics/pythonselect.PNG"
         alt="Excel Input">
</figure>

4. Click "Apply" then "OK"

NOTE: The following packages must be installed in the Python version being used; numpy, scipy, pandas, and matplotlib. 


## Quick Start

1. In RStudio, click the "Run App" button in the top right corner of the window that contains the code

![Start](DocumentPics/start.PNG)

2. Application opens to Input & Column Parameter view. User can select "Browse" below "Choose .xlsx File" to import preconfigured input file in the upper left, as shown. Or, a user can begin editing input values following steps 4 and 5 below.

![Import File](DocumentPics/Slide1.PNG)


3. (Optional) Change the parameters to match the specifications of your Granular Activated Carbon apparatus. This can be done by typing in a number or using the scroll wheel to increment the number up or down.

![Inputs](DocumentPics/Slide2.PNG)

4. (Optional) In the compounds tab, add chemicals and concentration points to match your interest. These can be added or edited by right clicking the data table.

![Adjust](DocumentPics/IonsTab.PNG)
![IonEdit](DocumentPics/Slide4.PNG)

5.	Click the Run Analysis button that’s at the bottom of the same side panel as the file import. A notification window will appear in the bottom right corner to either indicate a successful run or give an error message. 3 total chemicals with 2 concentration points takes about 10-30 seconds.

   *NOTE*:  Rows "kf", "dp", and "ds" can be added to the "Compound List" and will manually override the mass transfer coefficients. If you set their value to 0 or omit them, they will be recalculated automatically.  

![Run](DocumentPics/Slide5.PNG)

6.	The program will automatically switch to the Output tab when it is finished running or you can manually switch to the Output tab by clicking on it (There should be a loading “spinner” to let you know it’s running)

![waiting](DocumentPics/Slide6.PNG)

7.	Your graph will appear. You can export the data as an xlsx file along with the conditions you input.

![plot](DocumentPics/Slide7.PNG)
Output for Example_Multi.xlsx

8. The user can use the check boxes on the left to toggle the influent data and the effluent data (if available). Note that the user can toggle on and off individual traces on the graph by clicking on the desired data on the legend (data will be grayed out if it isn't displayed).

![plottraces](DocumentPics/Slide8.PNG)
Output for Example_Multi.xlsx

9. (Optional: Fitting the Data) The user can fit the K data to the any effluent data that is provided. To do this, simply click the 'Fit Data' button above the 'Save Data' button on the side bar panel. The fitted k data will then display on the 'Fitted Data' panel that the user can select on the header. Finally, a user may click 'Use Data' and this will replace the data in the Compounds tab with the fitted k data. The user can then rerun the program to view the updated results. There is an example fitter input file, called Example_fitter.xlsx, that can be found in the Examples folder.

![Fit](DocumentPics/Fit.PNG)
![Use](DocumentPics/Use.PNG)
![FitResults](DocumentPics/FitResults.PNG)
Output for Example_fitter.xlsx


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
|flowrate                       | The average flow rate through the column. PSDM only considers and average or steady-state condition, not variable flow. | gpm, cm^3/s, m^3/s, ft^3/s, mL/s, L/min, mL/min, mgd |
|diameter                     | The diameter of a cylindrical column.  | cm, mm, m, ft, in |
|tortuosity                   | Parameter of flow between curve and length. | N/A |
|units                        | Influent and Effluent Concentrations | ug, ng, mg |
|time                         | The units for time in the corresponding "data" sheet in the Excel-based files or "Concentration Points" table under Input>Ion's tab. |




### Compounds Tab

The Compounds tab contains information about the compounds to be simulated. The Compounds tab should contain any compounds found in the system. Compound characteristics (molecular weight, molar volume, boiling point, density, solubility, vapor pressure) are stored in the "Properties" tab in Excel-based input files and listed under "Compound List" within PSDM under the Input>Compounds tab. Information on K data within the system over is stored in the "Kdata" tab in Excel-based input files or "K Data" within PSDM under the Input>Compounds tab.


### Compound List
|  Input        	                |Column Name   | Description                                                                      |
|---            				    |---        |---                                                                               |
|Molecular Weight               |MW         |Molecular weight of ionic species. Gram per mol.                                          |           |
|Molar Volume                   |MolarVol       |(L/mol) Volume occupied by one mole of ion substance. |                     
|Boiling Point                  | BP          |(degC) Temperature in Celsius that the given chemical boils. |
|Density      					|Density    	| (g/mL) Mass per unit volume of ion.     |
|Solubility						| Solubility  		| (g/L) Amount of substance that will disolve in water. |
|Vapor Pressure  				| VaporPress       | (mmHg) Pressure of vapor form of substance. |      


### K Data

|Name                   |Variable      | Description            |
|---                    |---           |---                     |
|Freundlich K           | K            | (ug/g)(L/ug)^(1/n) Freundlich isotherm parameter: q = K * Ce^(1/n)                       |
|Freundlich 1/n         | 1/n          | (unitless) Exponent in Freundlich isotherm                         |
|Solid-phase capacity   | q            | (ug/g) Solid-phase capacity of media at equilibrium with liquid phase concentration                      |
|Time to breakthrough    | brk          | Days to breakthrough, provided by fitting function                       |
|Average influent concentration    | AveC         | Average influent provided by fitting function, can provide C0 if needed.                 |

### Concentration Points / Influent and Effluent Data

The Concentration Points table (stored in Excel-based files in the "data" sheet) represents time series of concentrations in the system over time (conc_units for an ion defines the concentration units for each column). NOTE: The duration of a simulation is specified by the largest time/last row in the Time column. Times in this table should be specified in ascending order, and a minimum of two (2) rows are needed to run a simulation (time 0, and run duration). If concentrations are the same in all rows, the simulation has a constant concentration, but variable concentrations can be modeled (where linear interpolation is used between points). Time units for the times in this table is specified on the "Input>Column Parameters" tab as time (day or hr). This data can be plotted as well by simply clicking the "Influent Data" button on the 'Output' page.

The effluent data represents the results from a previous experiment or results from a previous model. This data does not need to be present for the model to run. If the user would like to see this data on the plot, the data must be present with the influent in the 'data' sheet. When the analysis has been ran, the user can toggle on the "effluent data" button on the left hand side to see the data plotted. The units that the effluent data is in is represented by the time units as well in the column parameters. Whatever unit corresponds to the chemical in the 'concentration units' selection will be converted to the 'Output Units' on the output page. Moreover, if the user wants to plot the effluent data then the chemical must exist in both the 'ions' page, the 'cin page', and the 'effluent page'.

### Output
The Output tab provides a graphical output of results after a simulation is completed. The units for time and concentration can be adjusted via the dropdown menus in the left-hand column.

![output](DocumentPics/output.png)

The graphs will dynamically update to reflect the current selection. If c/c0 is selected, the output is scaled relative to c0 (the initial concentrations), which is the first row of the "Concentration Points" table. "Bed volumes x1000" presents time in thousands of bed volumes treated, where bed volume is the empty bed volume of the media in the system (relating to bed length if velocity is specified, or bed length and diameter if flow rate is specified).

The data can be exported by clicking the save button. This saves the data points to an excel file where the chemicals inputted into the analysis are the header and the columns are the corresponding concentration points with the first column being the corresponding time.

To export the graph, the user can hover over the graph with their cursor which will display the graph settings in the top right. Clicking the camera icon will open up a file explorer where the user can save the graphic.

![savegraph](DocumentPics/saveplot.png)

## Notes to the User

#### Selection of Collocation Points (nr and nz)
The parameters nr and nz control the size of the grid used to numerically solve the underlying differential equations during the simulation. Increasing nr and nz may increase the accuracy of simulations but doing so also makes them take longer to run. No analytical expression has been found for determining optimal grid dimensions for this class of problems, so selecting nr and nz may take some experimentation. Generally, the sharper the GAC zone is relative to the column length, the higher nz will need to be and the sharper the diffusion gradient in the resin beads becomes, the higher nr will need to be. In practice, it is rare for nr to be the controlling parameter for grid size, with nr=7 being accurate enough for most cases without unduly increasing computational cost. The parameter nz is more likely to need attention. Setting nz too low will often produce erroneous oscillations in the breakthrough curves. The illustration below shows simulations with (a) and without (b) these erroneous oscillations.

![nrnz](DocumentPics/Picture1.png)

If increasing nz does not smooth out the erroneous oscillations, there may be other problems with the simulation. In this case, the user is advised to double check the inputs for errors. If there are no errors in the inputs, it is possible the GAC zones in the requested simulations are simply too sharp for this numerical approximation to handle. Faced with this problem, the user may wish to consider reducing the empty bed contact time of the simulation or seek out an alternate method of solution such as an equilibrium-based column model.

#### Specifying Column Size and Flow Rate

The underlying model equations in this code use column length (bed depth), L, to define filter size and superficial (linear) flow velocity, v, to define flow rate of simulated systems. If both parameters are readily available to the user, they can be input directly selecting the “linear” radio button on the left side of the column specification section. In practice, flow in adsorption systems is often specified as a hydraulic loading rate (or “surface loading rate”) given in units of volumetric flow rate divided by area (for instance, with units of gpm/ft2 or gallons per minute per square foot). This specification is ultimately equivalent to specifying a superficial flow velocity and can also be entered directly for v provided appropriate units are selected from the corresponding drop-down menu. 

Occasionally design specification may include bed dimensions and empty bed contact time (EBCT) but omit flow information. In this case, the user can obtain a superficial flow velocity from the following formula:

![eq3](DocumentPics/eq3.png)

The column size and flow rate may also be defined in terms of L, bed diameter (d), and volumetric flow rate (flrt) by selecting the “volumetric” radio button.
A note on selection of flow convention: The entry field for the two conventions are independent. The values shown in disabled fields (gray backgrounds) are not updated to correspond to values entered using the other convention. Thus, switching between the radio buttons usually results in switching between two different systems.

**Note:** If flow rate is provided, diameter is required


## References


## Development Team
David Colantonio

Levi Haupert

Jonathan Burkhardt

Cole Sandlin
