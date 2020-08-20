Installation
============

Requirements
-------------

You need to have the following software properly installed in order to
build *SeismicMesh*:

* Python >= 3.0
.. note ::
    The file ``requirements.txt`` in the main directory indicates all the Python dependencies and their respective version numbers. These packages should be installed at compile time by setuptools
.. note ::
    Some users have experienced  problems with the `skfmm-fmm` Python package not being found. If this occurs, try uninstalling and then reinstalling this package after attempting installation of SeismicMesh.

* Pybind11 >= 2.5
.. note ::
    Pybind11 is shipped as a submodule thus you do not need to install manually. You will activate it below when performing the installation.

* C++ compiler (GNU or Intel) with support for std++14 or newer.

* Cmake >= 3.0

* CGAL >= 5.0.0

    * MPFR

    * GMP

    * Boost > 1.4.8


.. warning ::
    Be aware that apt-get does not install the necessary version >=5.0.0 of CGAL for this project. Thus, you will need to build CGAL manually by source using CMake instead. To do this, please follow these instructions https://doc.cgal.org/latest/Manual/installation.html

.. note ::
    CGAL requires Boost, MPFR and GMP. These packages may already be installed on your standard Linux box.




Clusters
-------------

If installing on a cluster with a local installation of ``CGAL`` and ``Boost``, you'll need to edit ``setup.py`` with the CMake arguments so to point the installation to the correct directories. Namely, in ``setup.py`` you'll have to edit the list called ``cmake_args`` to include ::

  -DCMAKE_CXX_COMPILER=+/PATH/TO/CPPCOMPILER

  -DBoost_INCLUDE_DIR=+/PATH/TO/BOOST/

  -DMPFR_LIBRARIES=+/PATH/TO/libmpfr.a

  -DMPFR_INCLUDE_DIR=+/PATH/TO/MPFR/include


Compilation
-------------

After installing all dependencies, perform the two steps in the main directory in the order that they appear below::

$  git submodule update --init --recursive
$  pip install -e .

.. note ::
    If you do not have administrative rights on your system, add the flag ``--user`` to end of the second command

Testing
-------

Testing is accomplished with `pytest`. The `pytest` package can be installed like so::

    pip install pytest

.. note ::
    If you do not have administrative rights on your system, add the flag ``--user`` to end of the pip command.

To test the installation, serial and parallel capabilites, you can use `pytest` from the `test/` directory::

$ cd tests/
$ pytest -m "serial" .
$ mpirun -np 2 pytest -m "parallel" .
