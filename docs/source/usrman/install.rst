Installation
============

Requirements
------------

You need to have the following software properly installed in order to
build *SeismicMesh*:

* Python 3.5 or above.

* C++ compiler (GNU or Intel) with support std++14 or newer.

* Cmake >= 3.0

* CGAL >= 5.0

.. note ::
    CGAL requires MPFR and GMP which may or may not already be installed on your standard Linux box.

* Boost > 1.4.8

.. warning ::

    Notice the file ``requirements.txt`` which indicates all the dependencies and their respective version numbers. If installing on a cluster with a local  installation of ``CGAL`` and ``Boost``, you'll need to edit ``setup.py`` with the CMake arguments so to point the installation to the correct directories. Namely, in ``setup.py`` you'll have to edit the list called ``cmake_args`` to include ::

  -DCMAKE_CXX_COMPILER=+/PATH/TO/CPPCOMPILER

  -DBoost_INCLUDE_DIR=+/PATH/TO/BOOST/

  -DMPFR_LIBRARIES=+/PATH/TO/libmpfr.a

  -DMPFR_INCLUDE_DIR=+/PATH/TO/MPFR/include


Installation :
==============================================

Perform the steps in the order they appear::

$  git submodule update --init --recursive
$  python setup.py develop

.. note ::
    If you do not have administrative rights on your system, add the flag ``--user``

Testing
-------

To quickly test the installation, you can use `pytest` from the `test/` directory::

$ cd tests/
$ pytest .
