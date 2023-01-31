Mesh2scattering
==============

Mesh2scattering is based on Mesh2HRTF and is an open-source project aiming an easy-to-use software package for the numerical calculation of scattering pattern and scattering and diffusion coeffients of any surface. In a nutshell, Mesh2scattering consists of three parts:

* Mesh2Input: prepares geometrical data and acoustic parameters,
* NumCalc: based on the input from Mesh2Input, it calculates the corresponding sound field
* Output2scattering: processes the output from NumCalc to scattering and/or diffusion coeffients.

Please notice that this project does not support HRTF post processing, use Mesh2HRTF instead.

Installation
============

Using mesh2scattering with Python
---------------------------------

* install miniconda, Visual Studio Codes, git
The Conda API needs to installed in the terminal, before it can be used.

* ``git clone git@github.com:ahms5/Mesh2scattering.git`` in your folder.

* ``conda create --name mesh2scattering python`` creates an virtual environment for mesh2scattering

* ``conda activate mesh2scattering`` actaives the environment

* ``pip install -e .`` install mesh2scattering in pip, ``-e`` option will make changes in mesh2scattering immediately available. This is helpful if you adapt parts of the code to your needs. If you do not want this, omit this option.

* ``pip install -r requirements_python.txt`` install all packages for develop

* ``pip install ipykernel`` if you want to use the interactive window

* go to mesh2scattering/Mesh2Input/mesh2input.py and replace the ``/path/to/mesh2scattering/mesh2scattering`` by the path to this directory containing Mesh2Input, NumCalc, Output2scattering folders.

setting up blender
------------------

* install and open blender

* go to Edit -> Preferences -> Add-ons -> install

* choose ``/path/to/mesh2scattering/mesh2scattering``

* activate the add on.

NumCalc
-------

for Linux:

* Install the C++ build essentials by running ``sudo apt-get install build-essential``
* Go into the NumCalc directory by running ``cd path/to/your/Mesh2HRTF/mesh2hrtf/NumCalc/src``
* Compile NumCalc by running make. It is now located in the folder ``mesh2hrtf/NumCalc/bin``
* Copy NumCalc to a folder in your program path: in the same directory run ``sudo cp NumCalc /usr/local/bin/``
* Now NumCalc can be used by running NumCalc (don't do this yet).

for Windows:

download the executable from the release.

Usage
=====

Mesh2Input - Set up the simualtion
----------------------------------

* open 'mesh2scattering/Mesh2Input/Tutorials/mesh2scattering.py'

* change the parameters as decribe

* now we need to run this script, therefore open 'mesh2scattering/Mesh2Input/create_projects.py'

* change the ``blender_executable`` to the path of the blener installation, on windows it would end with blender.blender_executable

* run the script and in the ``project_path`` folder will be 2 new folders calles ``sample`` and ``reference`` these folders contain two BEM projects which we can execute now.

NumCalc - run the simualtion
----------------------------

* open 'mesh2scattering/NumCalc/run_local.py'

* change the ``project_path`` and the ``numcalc_path`` according to the locations. 

* now we can run the file and wait ... the script will automaticly parallize the task as much as possible

Output2scattering - read and export the data
--------------------------------------------

* open the 'mesh2scattering/Output2scattering/export_results.py'

* change the variables on the top and run it.

done. now you can postprocess the data as you like.


Further Documentation
=====================
For further Documentation please checkout the original project to calculate HRTFs.

Documentation
    https://github.com/Any2HRTF/Mesh2HRTF/wiki

Installation
    https://github.com/Any2HRTF/Mesh2HRTF/wiki/Installation

General Tutorials
    https://github.com/Any2HRTF/Mesh2HRTF/wiki/Tutorials


References Mesh2HRTF
====================

For information about the BEM algorithm please refer to
    H\. Ziegelwanger, W. Kreuzer, and P. Majdak, `''Mesh2HRTF: An open-source software package for the numerical calculation of head-related transfer functions,'' <https://www.researchgate.net/publication/280007918_MESH2HRTF_AN_OPEN-SOURCE_SOFTWARE_PACKAGE_FOR_THE_NUMERICAL_CALCULATION_OF_HEAD-RELATED_TRANFER_FUNCTIONS>`_ Florence, Italy, Jul. 2015.

For general recomendations on the mesh and microphone area see
    H\. Ziegelwanger, P. Majdak, and W. Kreuzer, `''Numerical calculation of listener-specific head-related transfer functions and sound localization: Microphone model and mesh discretization,'' <https://doi.org/10.1121/1.4922518>`_ J. Acoust. Soc. Am., vol. 138, no. 1, pp. 208–222, Jul. 2015, doi: 10.1121/1.4922518.
