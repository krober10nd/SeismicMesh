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
    On some Linux systems, users may have to resort to `apt install python3-segyio` to installing segyio on their systems.

* pybind11 >= 2.5

* C++ compiler (GNU or Intel) with support for std++14 or newer.

* cmake >= 3.0

* CGAL >= 5.0.0 which requires:

    * MPFR

    * GMP

    * Boost > 1.4.8

.. note ::
    CGAL requires Boost, MPFR and GMP. These packages may already be installed on your standard Linux box.

.. warning ::
    Make sure your package manager is downloading CGAL >= 5.0 otherwise you will not be able to install *SeismicMesh*!




Compilation by source
----------------------

After installing all dependencies, perform ::

$  pip install -e .

.. note ::
    If you do not have administrative rights on your system, add the flag ``--user`` to the command above.

.. warning ::
    With this said, the preferred method of installation using pypi: pip install SeismicMesh

Testing
-------

Testing is accomplished with `tox`. The `tox` package can be installed like so::

    pip install tox

To test the installation, serial and parallel capabilities, you can use `tox` from the top directory of the package::

$ tox

Installation on Clusters
--------------------------

.. note::
    Make sure the CXX environment variable points to your intended compiler!

If installing on a cluster by source with a local installation of ``CGAL`` and ``Boost``, you'll need to edit ``setup.cfg`` with the CMake arguments so to point the installation to the correct directories. Namely, in ``setup.py`` you'll have to edit the list called ``cmake_args`` to include ::

  -DCMAKE_CXX_COMPILER=+/PATH/TO/CPPCOMPILER

  -DBoost_INCLUDE_DIR=+/PATH/TO/BOOST/

  -DMPFR_LIBRARIES=+/PATH/TO/libmpfr.a

  -DMPFR_INCLUDE_DIR=+/PATH/TO/MPFR/include
