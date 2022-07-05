.. highlight:: shell

Contributing
------------

Contributions are welcome, and they are greatly appreciated! Every little bit
helps, and credit will always be given.

You can contribute in many ways:

Types of Contributions
----------------------

Report Bugs and Submit Feedback
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The best way to report bugs of send feedback is to open an issue at https://github.com/Any2HRTF/mesh2hrtf/issues.

If you are reporting a bug, please include:

* Your operating system name and version.
* Any details about your local setup that might be helpful in troubleshooting.
* Detailed steps to reproduce the bug.

If you are proposing a feature:

* Explain in detail how it would work.
* Keep the scope as narrow as possible, to make it easier to implement.
* Remember that this is a volunteer-driven project, and that contributions
  are welcome :)

Fix Bugs or Implement Features
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Look through the GitHub issues for bugs. Anything tagged with "bug" or
"enhancement" is open to whoever wants to implement it. It might be good to
contact us first, to see if anyone is already working on it.

Write Documentation
~~~~~~~~~~~~~~~~~~~

Mesh2HRTF could always use more documentation, whether as part of the
official Mesh2HRTF docs, in docstrings, or even on the web in blog posts,
articles, and such.

Start Contributing
------------------

Ready to contribute? Here's how to set up `Mesh2HRTF` for local development.
Work on the numerical core `NumCalc` and the Matlab/Octave API can directly be
done. Working on the Python API requires a local copy of mesh2hrtf to install
the API for local development

1. Fork the `Mesh2HRTF` repo on GitHub.
2. Clone your fork locally and cd into the Mesh2HRTF directory::

    $ git clone https://github.com/Any2HRTF/Mesh2HRTF.git
    $ cd Mesh2HRTF

3. Install your local copy into a virtualenv. Assuming you have Anaconda or Miniconda installed, this is how you set up your fork for local development::

    $ conda create --name mesh2hrtf python
    $ conda activate mesh2hrtf
    $ conda install pip
    $ pip install -e .
    $ pip install -r requirements_python.txt


The latest work is contained in `develop`. For bug-fixes, enhancements, and new
ideas please create a new branch based on develop.

Testing
-------

Mesh2HRTF uses py ``pytest`` for testing. For all tests to work, you must
configure your blender path in `test_export` (variable ``blender_paths``).

- All tests are ran by

    $ pytest

- Run a single test with

    $ pytest tests/test_plot.py::test_line_plots

- Exclude tests (for example the time consuming test of plot) with

    $ pytest -k 'not plot'

- Create an html report on the test `coverage <https://coverage.readthedocs.io/en/coverage-5.5/>`_ with

    $ pytest --cov=. --cov-report=html

Releasing
---------

To release a new Mesh2HRTF version do the following

- Write the new version to the file VERSION
- Update HISTORY.rst (also include new contributors)
- Commit all changes to develop
- add a tag with the version number ``git tag <tagname>``, e.g. ``git tag v1.0.0``
- push the tag using ``git push origin --follow-tags``
- merge develop into main