==========================================================================
Mesh2Input: Tools and data to prepare the project folder
==========================================================================

mesh2input.py
=============

This file contains the main functionality to create the project folder which stores all the input for the numerical HRTF calculations by NumCalc. To this end,

1. Define this file as an AddOn in Blender
2. Prepare your mesh for the calculations
3. Start this AddOn (in Blender) to create the project folder.

EvaluationGrids
===============

Evaluation grids are used to define the sound-source positions for which the HRTFs are
calculated. The following is available:

- **Data**: pre-defined evaluation grids
- Scripts for reading, writing, and exporting from Blender custom evaluation grids

Materials
=========

Materials are used to define acoustic boundary conditions for a mesh. The following is available:

- **Data**: pre-defined materials
- Scripts for writing custom materials

Meshes
======

Here, we provide various tools to manipulate meshes. The following tools are available:

- **AssignmenMaterials**: A Blender AddOn to automatically assign the default materials *Skin*, *Left ear*, and *Right ear* to a mesh. Requires a centered head mesh.
- **CenterHead**: A Blender AddOn to semi-automatic center the mesh in the origin of coordinates and orient it along the interaural axis.
- **Data**: An example of a head mesh that can be used for testing Mesh2HRTF. This is the mesh also being used in the online tutorials.
- **Grading, Distance based** (discontinued): A C++ tools for a distance-based grading of a mesh in order to reduce the number of faces and thus to decrease the duration of the numerical calculations [1].
- **Grading, Hybrid** (recommended): A C++ tool for a distance-based grading considering the local curvature [2].

Tutorials
=========
Here, we provide python scripts to generate the projects folders for the `online tutorials <https://github.com/Any2HRTF/Mesh2HRTF/wiki/Unix_Tutorials>`_.
To generate the project folders, execute the scripts in Blender as described in the `scripting Section <https://github.com/Any2HRTF/Mesh2HRTF/wiki/Scripting>`_.


References
==========

1 H. Ziegelwanger, W. Kreuzer, and P. Majdak, "A-priori mesh grading for the numerical calculation of the head-related transfer functions," Appl. Acoust. 114, 99–110 (2016). doi:`10.1016/j.apacoust.2016.07.005 <https://doi.org/10.1016/j.apacoust.2016.07.005>`_

2 T. Palm, S. Koch, F. Brinkmann, and M. Alexa, "`Curvature-adaptive mesh grading for numerical approximation of head-related transfer functions <https://www.researchgate.net/publication/280007918_MESH2HRTF_AN_OPEN-SOURCE_SOFTWARE_PACKAGE_FOR_THE_NUMERICAL_CALCULATION_OF_HEAD-RELATED_TRANFER_FUNCTIONS>`_", in Fortschritte der Akustik – DAGA 2021 (Vienna, Austria, 2021) pp. 1111–1114.

