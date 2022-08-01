Mesh2HRTF
=========

Mesh2HRTF is an open-source project aiming at providing an easy-to-use software package for the numerical calculation of head-related transfer functions (HRTFs). It targets researchers in the field of binaural audio. In a nutshell, Mesh2HRTF consists of three parts:
- Mesh2Input: prepares geometrical data and acoustic parameters,
- NumCalc: based on the input from Mesh2Input, it calculates the corresponding sound field
- Output2HRTFs: processes the output from NumCalc to HRTFs saved in the SOFA format.

Documentation
=============

Documentation
    https://github.com/Any2HRTF/Mesh2HRTF/wiki

Installation
    https://github.com/Any2HRTF/Mesh2HRTF/wiki/Installation

General Tutorials
    https://github.com/Any2HRTF/Mesh2HRTF/wiki/Tutorials

Releases and Development
========================

Releases are available at https://github.com/Any2HRTF/Mesh2HRTF/releases. 

The ``master`` branch contains the latest stable code basis and is available for cloning, forking, or downloads. 

The ``develop`` branch contains the latest stable developer version. Other branches track the more recent development, which might be unstable. If you would like to contribute to Mesh2HRTF please check `here <https://github.com/Any2HRTF/Mesh2HRTF/blob/develop/CONTRIBUTING.rst>`_.

References
==========

For information about the BEM algorithm please refer to
    H\. Ziegelwanger, W. Kreuzer, and P. Majdak, `''Mesh2HRTF: An open-source software package for the numerical calculation of head-related transfer functions,'' <https://www.researchgate.net/publication/280007918_MESH2HRTF_AN_OPEN-SOURCE_SOFTWARE_PACKAGE_FOR_THE_NUMERICAL_CALCULATION_OF_HEAD-RELATED_TRANFER_FUNCTIONS>`_ Florence, Italy, Jul. 2015.

For general recomendations on the mesh and microphone area see
    H\. Ziegelwanger, P. Majdak, and W. Kreuzer, `''Numerical calculation of listener-specific head-related transfer functions and sound localization: Microphone model and mesh discretization,'' <https://doi.org/10.1121/1.4922518>`_ J. Acoust. Soc. Am., vol. 138, no. 1, pp. 208–222, Jul. 2015, doi: 10.1121/1.4922518.
