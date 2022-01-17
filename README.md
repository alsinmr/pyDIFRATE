# pyDIFRATE
Dynamics analysis using MD and NMR


There is NO INSTALLATION required for the code. Just place everything in a folder, navigate there, and run. However, python3 and the following modules must be installed from other sources (these are the tested versions, although other versions may work).

Python v. 3.7.3
numpy v. 1.17.2,
scipy v. 1.3.0,
pandas v. 0.25.1,
MDAnalysis v. 0.19.2,
matplotlib v. 3.0.3,
cv2 v. 4.1.0 (only for Landscape plotting)
ChimeraX v. 1.0

We recommend installing Anaconda: https://docs.continuum.io/anaconda/install/
The Anaconda installation includes numpy, scipy, pandas, and matplotlib. 

MDAnalysis is installed by running:
conda config --add channels conda-forge
conda install mdanalysis
(https://www.mdanalysis.org/pages/installation_quick_start/)

cv2 can be installed by running 
conda install -c conda-forge opencv
(https://anaconda.org/conda-forge/opencv)
(cv2 is not required for the main functionality and can be omitted for most users)


Additionally, for 3D visualization, we use ChimeraX. In order to use the 3D visualization, this must also be installed (https://www.rbvi.ucsf.edu/chimerax/download.html). The path to the executable must then be provided to pyDIFRATE(DR.chimeraX.set_chimera_path(path), replace path with the full path to the executable)


All files are copyrighted under the GNU General Public License. A copy of the license has been provided in the file LICENSE

Funding for this project provided by:

Deutsche Forschungsgemeinschaft (DFG) grant SM 576/1-1

European Social Funds (ESF) and the Free State of Saxony (Junior Research Group UniDyn, Project No. SAB 100382164)



Copyright 2022 Albert Smith-Penzel, Kai Zumpfe