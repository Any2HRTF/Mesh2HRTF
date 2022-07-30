======
Readme
======

Collection of tools for Mesh2HRTF.

BlenderExport
=============

The main Blender AddOn that is used to generate Mesh2HRTF projects for
numerical HRTF calculation.

EvaluationGrid
==============

Evaluation grids are used to define the positions for which the HRTFs are
calculated.

- Mesh2HRTF default evaluation grids
- Python, Matlab/Octave, and Blender Tools for generating custom evaluation grids
  for Mesh2HRTF

Materials
=========

Materials are used to define acoustic boundary conditions for a mesh.

- Mesh2HRTF default materials
- Python, Matlab/Octave, and Blender Tools for generating custom materials

MeshManipulation/AssignmenMaterials
===================================

Blender AddOn for automatically assigning the default materials *Skin*,
*Left ear*, and *Right ear* to a mesh. Requires a centered head mesh.

MeshManipulation/CenterHead
===========================

Blender AddOn for the semi automatic centering and orientation of a head mesh
in the origin of coordinates.


MeshManipulation/Grading*
=========================

C++ tools for grading (remeshing) a mesh before using it for numerical
calculations. Grading reduces the number of faces in a mesh and thus decreases
the time of the numerical calculations.

Hybrid (recommended)
  Mesh grading based on the local curvature and distance from the ear as described by Palm et al. [1].

Distance Based (not maintained)
  Distance based mesh grading as described by Ziegelwanger et al. [2].

References
==========

1 T. Palm, S. Koch, F. Brinkmann, and M. Alexa, "`Curvature-adaptive mesh grading for numerical approximation of head-related transfer functions <https://www.researchgate.net/publication/280007918_MESH2HRTF_AN_OPEN-SOURCE_SOFTWARE_PACKAGE_FOR_THE_NUMERICAL_CALCULATION_OF_HEAD-RELATED_TRANFER_FUNCTIONS>`_", in Fortschritte der Akustik – DAGA 2021 (Vienna, Austria, 2021) pp. 1111–1114.

2 H. Ziegelwanger, W. Kreuzer, and P. Majdak, "A-priori mesh grading for the numerical calculation of the head-related transfer functions," Appl. Acoust. 114, 99–110 (2016). doi:`10.1016/j.apacoust.2016.07.005 <https://doi.org/10.1016/j.apacoust.2016.07.005>`_
