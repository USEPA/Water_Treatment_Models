# PSDM

The Pore and Surface Diffusion Model (PSDM) was converted into Python from the Fortran version used in AdDesignS(TM) from the **Environmental Technologies Design Option Tool** (ETDOT) from Michigan Technological University (MTU). (v1.0 - D.R. Hokanson, D.W. Hand, J.C. Crittenden, T.N. Rogers, E.J. Oman)

3 Python files are used in this tool.

* PSDM.py
* PSDM_functions.py (automatically imported by PSDM.py)
* PSDM_tools.py (must import to access additional capabilities)

See [Example_TCE.py](Example_TCE.py) for a few basic examples of how to run the PSDM code. (Example_TCE.xlsx is the input file for this example.)

See [Example_Multi.py](Example_Multi.py) for using the multi-component competitive capability (Example_Multi.xlsx is the example input.)

See [Example_isotherm.py](Example_isotherm.py) for using the isotherm fitting function.

Additional functions are available - Examples/documentation will be available in the future.




# Additional Information
See also tools found at https://github.com/USEPA/Environmental-Technologies-Design-Option-Tool

This repository is released under the [MIT License](../LICENSE.md).

EPA Disclaimer
==============
The United States Environmental Protection Agency (EPA) GitHub project code is provided on an "as is" basis and the user assumes responsibility for its use. EPA has relinquished control of the information and no longer has responsibility to protect the integrity , confidentiality, or availability of the information. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recomendation or favoring by EPA. The EPA seal and logo shall not be used in any manner to imply endorsement of any commercial product or activity by EPA or the United States Government.

By submitting a pull request, you make an agreement with EPA that you will not submit a claim of compensation for services rendered to EPA or any other federal agency. Further, you agree not to charge the time you spend developing software code related to this project to any federal grant or cooperative agreement.