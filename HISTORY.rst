History
=======

v1.0.0
-------
* Project export from Blender (handled by the Blender plugin `Mesh2Input/exportMesh2HRTF.py`)
	* Upgraded to Blender >= 2.80.0
	* Better organization and modularization of the Blender plugin
	* Re-design the parameters, appearance, and in-app documentation of the export menu to be less prone to erroneous input
		* The source type is now defined by a single drop down menu
		* Source properties (point source position and plane wave direction) are now obtained from Objects in the Blender Scene and do not have to be entered manually anymore.
		* Added an option for referencing the HRTF to a source in the origin of coordinates according to the HRTF definition
		* Referencing is now also supported for point sources
		* Add an option to calculate HRIRs from single sided HRTF spectra
		* More flexible selection of the frequencies that are simulated
		* Multiple evaluation grids can be entered in a single text field
		* Added handling for custom definitions of materials (boundary conditions)
		* Moved options for parallelization from the export menu to NumCalc (see below)
		* Removed undocumented and unfinished options fro near-field calculation and frequency dependent meshes
	* Boundary conditions can now be frequency dependent
	* Files containing custom evaluation grids and material data can be located outside the Mesh2HRTF repository
* NumCalc (contained in NumCalc/Source)
	* Introduced command line parameters `istart` and `iend` to select a range of frequencies for simulation and to ease parallelization
	* Introduce command line parameter `nitermax` to control the maximum number of iterations
	* Reduce the default number of maximum iterations to 250
	* Minor bug fixes and stability improvements
* Added a Python API for processing NumCalc output and save HRTF and HRIR data (contained in `Output2HRTF/Python`)
	* Installable via pip
	* Full online documentation
	* Added function to generate a project report and notify in case of issues and/or errors that occurred during the NumCalc simulation
	* Added flexible plot function for quick inspection of the results
	* Added Python tools to read and generate custom evaluation grids
	* Added function to merge results from multiple sources (e.g. left and right ear) into a single SOFA file
* Improved structure of the output data (Sofa files, project reports, exports, plots)
	* Data is now stored in a separate folder `Output2HRTF`
	* Data is named according to the evaluation girds and object meshes
	* Data for multiple evaluation grids is now stored in separate files
	* Frequencies in SOFA files now contain decimal values
* Added testing
	* All three arts of Mesh2HRTF (Project export, NumCalc, and Output2HRTF) are tested using pytest to improve and monitor the code quality
	* The Matlab/Octave API is not tested. New users are recommended to use the Python API


v0.4.0
------
* this will be the last release supporting Blender < 3.8
* new directory structure

v0.3.2
------
* big fix and improvements in PreProcessing/MeshGrading (ticket #25, commit r38)
* bug fix in ExportMesh2HRTF.py (Tickets #13, #15, and #23 on sourcefourge)
* fixed a bug Output2HRTF_Main.m. Now the low frequency HRTF should at 0 dB if using reciprocal simulation and if setting reference = true in Output2HRTF.m, which is auto-generated when exporting from blender.

v0.3.1
------
* bug fix in NumCalc

v0.3.0
------
* New directory structure
* Pascal-case naming of the files
* Small bugfixes in the scripts

v0.2.0
------
* Mesh2Input:
	* MaterialAssignment.py: A Python script that can be loaded into Blender to center the head mesh in the coordinate system
	* MeshCentering.py: A Python script that can be loaded into Blender to automatically assign the materials 'Skin', 'Left ear', and 'Right ear'
	* export_mesh2hrtf.py: Bug fix to correctly export data for calculating the left ear, right ear and both ears.
	* EvaluationGrids (Matlab):
		* Arbitrary user defined spatial grids can now be generated (see the code in demo.m)
		* 'User' in 'Mesh2Input/Data/Evaluation Grids' renamed to 'Custom' because 'User' is a reserved variable in Blender/Python
		* Evaluation grids can be plotted with the Matlab code
* NumCalc: MS VS Solution added to compile NumCalc on Windows.
* Output2HRTF:
	* Output2HRTF_Main.m: Added optional referencing of HRTFs if calculated reciprocally to achieve that the low frequency magnitude of the HRTFs is 0 dB. This is done by dividing the complex pressure by the area of the ear elements (radiating element), compensating for the velocity of the radiating element, and by a division of the complex pressure with the pressure of a point source in the origin of coordinates. (export_mesh2hrtf.py writes the area of the radiating elements, and the flag for referencing to Output2HRTF.m)
	* Output2HRTF_Main.m: Big fix to correctly export SOFA files with data for the left ear, right ear, and both ears.
* Mesh-grading tool moved to Mesh2Input

v0.1.3 (1.6.2018)
-----------------
* mesh-grading plugin for Open Flipper added
* Output: various bug fixes
* Output: Paraview scripts added
* NumCalc: Dissolve tiny clusters and add their elements to next bigger cluster. This seems to enhance the stability of the MLFMM.

v1.0.2 (18.6.2016)
------------------
* initial commit

Contributors
============

* Fabian Brinkmann, AC-TUB
* Slim Ghorbal, IETR/CentraleSupélec, Mimi
* Junaid Khan, FH Technikum Wien
* Michael Kalcher, ARI
* Wolfgang Kreuzer, ARI
* Piotr Majdak, ARI
* Robert Pelzer, AC-TUB
* Johan Pauwels, ICL
* Jeffrey Thomsen, AC-TUB
* Filip Tsai, KTH
* Oliver Weissbarth, CG-TUB
* Harald Ziegelwanger, ARI

**Institutions**

* ARI: Austrian Research Institute, Austian Academy of Sciences, Vienna
* AC-TUB: Audio Communication Group, Technical University of Berlin
* CG-TUB: Computer Graphics Group, Technical University of Berlin
* FH Technikum Wien: University of Applied Sciences, Technikum Wien
* ICL: Imperial College London, UK
* KTH: Royal Institute of Technology, Stockholm, Sweden
* Mimi: Mimi Hearing Technologies, Berlin, Germany