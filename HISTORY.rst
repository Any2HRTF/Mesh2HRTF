History
=======

v0.0.1
------
* Fork from Mesh2HRTF
* Remove blender and create projects from stl meshes in Python
* change project structure
* remove direct sound from NumCalc
* read simulation results and export to Sofa
* calculate scattering coefficients and save in 

Mesh2HRTF
=========

v1.0.0
-------
* Project export from Blender (handled by the Blender plugin `Mesh2Input/mesh2input.py`)
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
	* Added NumCalc manager for automatic parallelization of frequency steps
	* Bugfix: when NC.inp does not contain an "END" of file, NumCalc now throws one error and aborts programme execution (was infinite loop resulting in large output files)
	* Introduced command line parameters `istart` and `iend` to select a range of frequencies for simulation and to ease parallelization
	* Introduce command line parameter `nitermax` to control the maximum number of iterations
	* Introduce command line parameter `estimate_ram` for an a priori estimate of the RAM consumption per frequency step
	* Reduce the default number of maximum iterations to 250
	* Minor bug fixes and stability improvements
* Added a Python API for processing NumCalc output and save HRTF and HRIR data
	* Installable via pip
	* Full online documentation
	* Added function to generate a project report and notify in case of issues and/or errors that occurred during the NumCalc simulation
	* Added flexible plot function for quick inspection of the results
	* Added Python tools to read and generate custom evaluation grids
	* Added function to merge results from multiple sources (e.g. left and right ear) into a single SOFA file
	* Added function to write boundary conditions to material files
* Improved structure of the output data (Sofa files, project reports, exports, plots)
	* Data is now stored in a separate folder `Output2HRTF`
	* Data is named according to the evaluation girds and object meshes
	* Data for multiple evaluation grids is now stored in separate files
	* Frequencies in SOFA files now contain decimal values
* Added testing
	* All three parts of Mesh2HRTF (Project export, NumCalc, and Output2HRTF) are tested using pytest to improve and monitor the code quality
	* The Matlab/Octave API is not tested. New users are recommended to use the Python API
* General
	* Unified names of functions across programming languages
	* Updated project wiki and moved to github


v0.4.0
------
* this will be the last release supporting Blender < 2.8
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

Mesh2scattering Developers
==========================

Mesh2scattering is currently maintained and developed by
Anne Heimes

Mesh2HRTF Developers
====================

Mesh2HRTF is currently maintained and developed by
Piotr Majdak (Conceptualization),
Fabian Brinkmann (Python & Matlab API, Blender Export, Testing, Documentation),
Wolfang Kreuzer (NumCalc, Documentation),
Katharina Pollack (Matlab API, Documentation)

Contributors
============

The following persons contributed to Mesh2HRTF, named in reverse chronological
order

* Jeffrey Thomsen (2022): Testing and documentation
* Sergejs Dombrovskis (2022): Initial NumCalc manager version, documentation and tutorials
* Timon Palm (2021): Hybrid mesh grading tool
* Sebastian Koch (2021): Hybrid meh grading tool
* Junaid Khan (2020): Bugfixes and restructuring
* Oliver Weissbarth (2020): Update of the OpenFlipper mesh grading plug-in
* Slim Ghorbal (2019): Improved Blender export
* Robert Pelzer (2018): Blender AddOns for head centering and material assignment
* Harald Ziegelwanger (2016-2018): Initial development of Mesh2HRTF
* Z. S. Chen (2016): Initial development of NumCalc
* Michael Kalcher
* Johan Pauwels

**Institutions**

* Austrian Research Institute, Austian Academy of Sciences, Vienna
* Audio Communication Group, Technical University of Berlin
* Computer Graphics Group, Technical University of Berlin
* University of Applied Sciences, Technikum Wien
* Imperial College London, UK
* Royal Institute of Technology, Stockholm, Sweden
* Mimi Hearing Technologies, Berlin, Germany
