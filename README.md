# SWB2 water budget model for Tutuila American Samoa


See page at: 
https://uh-wrrc-swb-model.github.io/SWB2-Tutuila/

### Abstract
A water budget approach using SWB2, a soil water-balance model, was applied to the island of Tutuila in American Samoa. The primary objective for the model was to calculate spatially and temporally distributed net-infiltration, which directly controls groundwater recharge rate. This information is essential for assessing water resources availability on islands such as Tutuila where groundwater makes up the majority of drinking and municipal water supplies. Other water budget components such as evapotranspiration, canopy interception, runoff, and mountain front recharge were also quantified with the SWB2 model for average present-day climate conditions. Additionally, the potential effects of future climate change on water resources availability were simulated by integrating dynamically downscaled climate predictions for 2080 to 2099 derived from externally supplied global climate model results. Notable improvements in this model over previously developed Tutuila water budget models include flow-routing based on land topography, inclusion of the mountain front recharge process, and consideration of direct net infiltration from anthropogenic sources such as on-site wastewater units and leaking water delivery lines. Model results indicated approximately 54% of Tutuila’s rainfall infiltrates as groundwater recharge, 8% is lost to canopy evaporation, another 15% is lost to evapotranspiration from soils, and 21% is removed through surface water features as stormflow-runoff.  The model was able to simulate these processes with a high-spatial and temporal resolution with a 20 by 20 m grid-cell size, and a daily-resolution output time step. Climate scenarios suggested an increase in net-infiltration of 17 to 27% may be expected by the end of the century depending on the emissions scenario used.


### Model run instructions

To run this code, the user should copy the entire repository, keeping the directory structure intact. GITHUB prevents large file uploads and also binary files to be included in the repository structure, and thus there are two files that the user will need to manually download from the github "releases" tab within this repository. These are: 

1) "Land_use_wRO_codes.shp" needs to be downloaded and then physically moved into to the "Raw_GIS_Data//Land_use" folder in the users copy of the project directory. 

2) "swb2.exe" needs to be downloaded and then physically moved into to the "Run" folder in the users copy of the project directory. .


To run the model:
The model code is contained in a Jupyter Notebook ("SWB2_Tutuila_model.ipynb") contained in the "SWB2-Tutuila/Model_workspace/Run/" folder. If the user wishes to run the model outside of an .ipynb format, the notebook is formatted as a single cell for easy copying to a .py file or other format. 

Before running, ensure relative paths to the "Raw_GIS_Data" and "Std_input" folders are correct ensure the swb.exe executable file is installed in the same directory that contains this notebook ("SWB2-Tutuila/Model_workspace/Run/") ensure the control file "Tutuila200_controlFile.ctl" is also located in the "Run" folder.  When executed, the model will create the appropriate directory structure and will create output folders with model results in gridded and tabular summarized format, as well as figures in .tif format.

The user needs to have all of the modules specified in the first code block installed within the active python environment. If any modules are not installed, the can be either downloaded from https://anaconda.org/conda-forge, or if using anaconda, the syntax to have anaconda perform the install can be typed into the conda command prompt.

&nbsp;

&nbsp;

If you wish to contact the author I can be reached at cshuler@hawaii.edu

&nbsp;

&nbsp;


# Disclaimer
This script is provided as open-source software on the condition that neither Chris Shuler Hydrologic LLC nor the UH Water Resource Research Center shall be held liable for any damages resulting from the authorized or unauthorized use of the information. No warranty, expressed or implied, is made by the authors as to the accuracy and functioning of the program and related program material nor shall the fact of distribution constitute any such warranty and no responsibility is assumed by the authors in connection therewith. This information is preliminary or provisional and is subject to revision. This software is provided "AS IS." Note that sensitive information, or datasets that are not publically available, are not posted in raw forms. The model code is licensed under the GNU General Public License v3.0 which is an open-access license designed to explicitly affirm any user’s unlimited permission to run, copy, and use the unmodified code from this repository. Please note that some raw datasets used in this work are not owned by the authors and may be subject to other licenses or conditions.
