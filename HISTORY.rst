History
=======

v1.1.2 (21 May 2024)
--------------------
* Fixed a bug when exporting a Mesh2HRTF Project from Blender that uses a plane wave as sound source (PR #108)
* Update testing to work with all supported Blender versions (PR #189, #129)
* Update `inspect_sofa_files` from the Python API to work with the latest pyfar version (PR #125)
* Improve tooltip help in Mesh2HRTF Export GUI in Blender (PR #124)
* Add `__version__` attribute to Python API (PR #111)

v1.1.1 (09 June 2023)
---------------------
* When creating a Mesh2HRTF project from Blender, the Blender project is now saved to the project folder as a compressed copy. Previously it was not compressed and saved directly. The latter came with the danger of making changes to Blender project after creating the Mesh2HRTF project without noting. In this case the Blender file might not reflect the settings of the rest of the project folder (PR #105).
* Fix a bug in detecting issues in NumCalc simulations during post-processing using `ouput2hrtf` from the Python API. Errors were not reported in case `NumCalc -estimate_ram` was called before using `manage_numcalc` from the Python API.
* Improve the output of `manage_numcalc` from the Python API.
* Fix a bug in detecting the free RAM in `manage_numcalc` from the Python API.

v1.1.0 (11 Mai 2023)
--------------------
* Add dockerfile running the Mesh2HRTF Python API, `NumCalc`, and `hrtf_mesh_grading`
* Fixed a bug in installing the Mesh2HRTF Python API
* Fixed a bug in `manage_numcalc` if running multiple Mesh2HRTF projects
* Fix order of HRIRs/HRTFs plotted by `inspect_sofa_files` when choosing the median or frontal plane

v1.0.0 (28 April 2023)
----------------------
* Mesh2Input (Project export from Blender handled by the Blender plugin `Mesh2Input/mesh2input.py`):
	* Upgraded to support Blender versions >= 2.80.0
	* Improved organization and modularization of the Blender plugin
	* Re-design the parameters, appearance, and in-app documentation of the Export menu to be less prone to erroneous input
		* Source type can be defined by a single drop down menu
		* Source properties (point source position or plane wave direction) obtained from Objects in the Blender Scene (separate manual input not required).
		* Support for referencing the HRTFs when calculated for point sources
		* Option for referencing the HRTF to a source placed in the origin of coordinates
		* Option to calculate HRIRs from single-sided HRTF spectra
		* More flexible selection of the simulated frequencies
		* Support of multiple evaluation grids in a single text field
		* Support of custom definitions of materials (defining the boundary conditions)
		* Options for parallelization moved from the Export menu to NumCalc (see below)
		* Clean up: Undocumented and unfinished options for the near-field calculation and frequency-dependent meshes removed
	* Support of frequency-dependent boundary conditions
	* Support for the custom evaluation grids and material data being located outside the Mesh2HRTF directory
	* Improved detection and display of errors in the user input.
* NumCalc:
	* `NumCalc/manage_numcalc.py` added: a NumCalc manager for automatic parallelization of frequency steps
	* Bugfix: for `NC.inp` without "END", NumCalc throws an error and aborts the execution preventing an infinite loop
	* New command line parameters: `istart` and `iend` select the range of frequencies for the simulation to ease the parallelization
	* New command line parameter: `nitermax` controls the maximum number of iterations
	* New command line parameter: `estimate_ram` provides an a-priori estimation of the RAM consumption per frequency step
	* Default number of maximum iterations reduced to 250
	* Minor bug fixes and stability improvements
* Output2HRTF:
	* New Python API for processing NumCalc output and save HRTF and HRIR as SOFA files
		* Installable via pip
		* Full online documentation
		* Added function to generate a project report and notify in case of issues and/or errors that occurred during the NumCalc simulation
		* Added flexible plot function for quick inspection of the results
		* Added Python tools to read and generate custom evaluation grids
		* Added function to merge results from multiple sources (e.g. left and right ear) into a single SOFA file
		* Added function to write boundary conditions to material files
	* Improved structure of the output data (SOFA files, project reports, exports, and plots):
		* Data stored in a separate folder `Output2HRTF`
		* Data named according to the evaluation grids and object meshes
		* Data for multiple evaluation grids stored in separate files
		* Frequencies in SOFA files stored as decimal numbers
* General:
	* Testing: NumCalc and the Python-based parts of Mesh2HRTF (Project export, and the Python-based part of Output2HRTF) are automatically tested using pytest to improve and monitor the code quality. The Matlab/Octave API is tested manually.
	* Unified names of functions across the programming languages
	* Project Wiki migrated to Github and updated

v0.5.0 (July 2022)
------------------
* preparation for the upgrade to Mesh2HRTF 1.x
* license changed to the EUPL 1.2
* this is the last Mesh2HRTF version supporting Blender versions up to 2.79

v0.4.0
------
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

v0.2.0 (2018)
-------------
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

v0.1.3 (2015)
-----------------
* mesh-grading plugin for Open Flipper added
* Output: various bug fixes
* Output: Paraview scripts added
* NumCalc: Dissolve tiny clusters and add their elements to next bigger cluster. This seems to enhance the stability of the MLFMM.

v0.1.2 (2015)
------------------
* initial commit and release via SourceForge

v0.1.1 (2014)
* initial version by Harald Ziegelwanger, Piotr Majdak, and Wolfgang Kreuzer

Mesh2HRTF Developers
====================

Mesh2HRTF is currently maintained and developed by
* Piotr Majdak (Conceptualization, Maintainence),
* Fabian Brinkmann (Python & Matlab API, Blender Export, Testing, Documentation),
* Wolfang Kreuzer (NumCalc, Documentation),
* Katharina Pollack (Matlab API, Documentation)

Contributors
============

The following persons contributed to Mesh2HRTF (in reverse chronological order):

* Tim Wennemann (2023): Update for CenterHead and AssignMaterial scripts
* Jeffrey Thomsen (2022): Testing and documentation
* Sergejs Dombrovskis (2022): Initial NumCalc manager version, documentation and tutorials
* Johan Pauwels (2022): various
* Timon Palm (2021): Hybrid mesh grading tool
* Sebastian Koch (2021): Hybrid meh grading tool
* Junaid Khan (2020): Bugfixes and restructuring
* Oliver Weissbarth (2020): Update of the OpenFlipper mesh grading plug-in
* Slim Ghorbal (2019): Improved Blender export
* Robert Pelzer (2018): Blender AddOns for head centering and material assignment
* Michael Kalcher (2016): various
* Harald Ziegelwanger (2013-2015): Initial development of Mesh2HRTF
* Z. S. Chen (until 2012): Initial development of NumCalc

**Involved Institutions**

* Acoustics Research Institute, Austrian Academy of Sciences, Vienna, Austria
* Audio Communication Group, Technical University of Berlin, Germany.
* Computer Graphics Group, Technical University of Berlin, Germany.
* University of Applied Sciences, Technikum Wien, Austria.
* Imperial College London, United Kingdom.
* Royal Institute of Technology, Stockholm, Sweden.
* Mimi Hearing Technologies, Berlin, Germany.
