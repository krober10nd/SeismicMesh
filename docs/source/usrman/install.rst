Installation
============

Requirements
------------

You need to have the following software properly installed in order to
build *SeismicMesh*:

* Python >= 3.0
.. note ::
    The file ``requirements.txt`` in the main directory indicates all the Python dependencies and their respective version numbers.

* Pybind11 >= 2.5
.. note ::
    Pybind11 is shipped as a submodule thus you do not need to install manually. You will activate it below when performing the installation.

* C++ compiler (GNU or Intel) with support for std++14 or newer.

* Cmake >= 3.0

* MPFR

* GMP
.. note ::
    CGAL requires MPFR and GMP. These packages may already be installed on your standard Linux box.

* CGAL >= 5.0
.. warning ::
    Be aware that apt-get does not install the necessary version >=5.0 of CGAL. It is strongly recommended to build CGAL manually by source using Cmake instead. Please follow these instructions https://doc.cgal.org/latest/Manual/installation.html

* Boost > 1.4.8

Clusters
-------------

If installing on a cluster with a local  installation of ``CGAL`` and ``Boost``, you'll need to edit ``setup.py`` with the CMake arguments so to point the installation to the correct directories. Namely, in ``setup.py`` you'll have to edit the list called ``cmake_args`` to include ::

  -DCMAKE_CXX_COMPILER=+/PATH/TO/CPPCOMPILER

  -DBoost_INCLUDE_DIR=+/PATH/TO/BOOST/

  -DMPFR_LIBRARIES=+/PATH/TO/libmpfr.a

  -DMPFR_INCLUDE_DIR=+/PATH/TO/MPFR/include


Compilation
-------------

After installing all dependencies, perform the two steps in the main directory in the order they appear::

$  git submodule update --init --recursive
$  python setup.py develop

.. note ::
    If you do not have administrative rights on your system, add the flag ``--user`` to end of the second command

Testing
-------

Testing is accomplished with `pytest`. The `pytest` package can be installed like so::

    pip install pytest

To quickly test the installation, serial and parallel capabilites, you can use `pytest` from the `test/` directory::

$ cd tests/
$ pytest -m "serial" .
$ mpirun -np 2 -m "parallel" .
