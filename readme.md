#General information

Mesh2HRTF was developed at the Acoustics Research Institute. Mesh2HRTF is an open-source project aiming at providing an easy-to-use software package for the numerical calculation of HRTFs. It targets researchers in the field of binaural audio. In a nutshell, Mesh2HRTF simply reads geometrical data, calculates the corresponding sound field and outputs HRTFs. To support multiple computer platforms, the concept of Mesh2HRTF is to focus on a command-line tool, which forms the numerical core, i.e., an implementation of the 3-dimensional Burton-Miller collocation BEM coupled with the multi-level fast multipole method (ML-FMM), and to provide add-ons for existing cross-platform applications for the preprocessing of geometrical data and for the visualization of results.

#Relases
Mesh2HRTF releases are also contained in the `releases` folder.

#Installation
https://sourceforge.net/p/mesh2hrtf/wiki/Installation/

#Tutorials
https://sourceforge.net/p/mesh2hrtf/wiki/Tutorials/

#Development
The `master` branch always contains the latest release, while the current work is contained in `develop`. For bug-fixes, enhancements, and new ideas please use `feature_*` branches.

#References

For information about the BEM algorithm please refer to:

[H. Ziegelwanger, W. Kreuzer, and P. Majdak, “Mesh2HRTF: An open-source software package for the numerical calculation of head-related transfer functions,” Florence, Italy, Jul. 2015.](https://www.researchgate.net/publication/280007918_MESH2HRTF_AN_OPEN-SOURCE_SOFTWARE_PACKAGE_FOR_THE_NUMERICAL_CALCULATION_OF_HEAD-RELATED_TRANFER_FUNCTIONS)

The mesh-grading is described in:

[H. Ziegelwanger, W. Kreuzer, and P. Majdak, “A-priori mesh grading for the numerical calculation of the head-related transfer functions,” Appl. Acoust., vol. 114, pp. 99–110, Dec. 2016, doi: 10.1016/j.apacoust.2016.07.005.](https://doi.org/10.1016/j.apacoust.2016.07.005)

For general recomendations on the mesh and microphone area see:

[H. Ziegelwanger, P. Majdak, and W. Kreuzer, “Numerical calculation of listener-specific head-related transfer functions and sound localization: Microphone model and mesh discretization,” J. Acoust. Soc. Am., vol. 138, no. 1, pp. 208–222, Jul. 2015, doi: 10.1121/1.4922518.](https://doi.org/10.1121/1.4922518)
