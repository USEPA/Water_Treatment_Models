# Ion Exchange Models

Preliminary versions of ion exchange models.

For an overview of ion exchange column modeling, consult:
* Slater, M.J., 2013. Principles of ion exchange technology. Butterworth-Heinemann.
* Helfferich, F. G. (1995). Ion exchange. Courier Corporation.

For details on the numerical method of solution (Orthogonal Collocation), consult Crittenden, J. C., Hutzler, N. J., Geyer, D. G., Oravitz, J. L., & Friedman, G. (1986). Transport of organic compounds with saturated groundwater flow: Model development and parameter sensitivity. Water Resources Research, 22(3), 271-284.

## IXPY
Contains a Homogeneous Surface Diffusion Model (HSDM) and a Pore Surface Diffusion Model (PSDM) for ion exchange (IX) in fixed-bed columns. An example input file and use case can be found in the test folder.

Provided the dependencies in [requirements.txt](requirements.txt) are satisfied, the models can be run from the command line without being installed:
```
python -m /path/to/ixpy {hsdm or psdm} input_file output_file
```

After installation, the models can be run without specifying the path to the hsdmix module:
```
python -m ixpy {hsdm or psdm} input_file output_file
```

# Testing

The unit tests for this package can be run from the command line:
```
python -m unittest discover -s /path/to/test -v
```

# Installation

Provided the dependencies in [requirements.txt](requirements.txt) are satisfied, this package can be installed from the command line:
```
python -m pip install /path/to/IonExchangeModel
```

The package can also be uninstalled from the command line:
```
python -m pip uninstall ion-exchange-models
```

# Status 
All code in this repository is being provided in a "draft" state and has not been reviewed or cleared by US EPA. This status will be updated as models are reviewed.

EPA Disclaimer
==============
The United States Environmental Protection Agency (EPA) GitHub project code is provided on an "as is" basis and the user assumes responsibility for its use. EPA has relinquished control of the information and no longer has responsibility to protect the integrity , confidentiality, or availability of the information. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recomendation or favoring by EPA. The EPA seal and logo shall not be used in any manner to imply endorsement of any commercial product or activity by EPA or the United States Government.

By submitting a pull request, you make an agreement with EPA that you will not submit a claim of compensation for services rendered to EPA or any other federal agency. Further, you agree not to charge the time you spend developing software code related to this project to any federal grant or cooperative agreement.
