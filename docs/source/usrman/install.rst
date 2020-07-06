Installation
============

Requirements
------------

You need to have the following software properly installed in order to
build *MPI for Python*:

* A working MPI implementation, preferably supporting MPI-3 and built
  with shared/dynamic libraries.

  .. note::

     If you want to build some MPI implementation from sources,
     check the instructions at :ref:`building-mpi` in the appendix.

* Python 2.7, 3.5 or above.

  .. note::

     Some MPI-1 implementations **do require** the actual
     command line arguments to be passed in :c:func:`MPI_Init()`. In
     this case, you will need to use a rebuilt, MPI-enabled, Python
     interpreter executable. *MPI for Python* has some support for
     alleviating you from this task. Check the instructions at
     :ref:`python-mpi` in the appendix.


Using **pip** or **easy_install**
---------------------------------

If you already have a working MPI (either if you installed it from
sources or by using a pre-built package from your favourite GNU/Linux
distribution) and the :program:`mpicc` compiler wrapper is on your
search path, you can use :program:`pip`::

  $ [sudo] pip install mpi4py

or alternatively *setuptools* :program:`easy_install` (deprecated)::

  $ [sudo] easy_install mpi4py

.. note::

   If the :program:`mpicc` compiler wrapper is not on your
   search path (or if it has a different name) you can use
   :program:`env` to pass the environment variable :envvar:`MPICC`
   providing the full path to the MPI compiler wrapper executable::

     $ [sudo] env MPICC=/path/to/mpicc pip install mpi4py

     $ [sudo] env MPICC=/path/to/mpicc easy_install mpi4py



Testing
-------

To quickly test the installation::

  $ mpiexec -n 5 python -m mpi4py.bench helloworld
  Hello, World! I am process 0 of 5 on localhost.
  Hello, World! I am process 1 of 5 on localhost.
  Hello, World! I am process 2 of 5 on localhost.
  Hello, World! I am process 3 of 5 on localhost.
  Hello, World! I am process 4 of 5 on localhost.

If you installed from source, issuing at the command line::

  $ mpiexec -n 5 python demo/helloworld.py

or (in the case of ancient MPI-1 implementations)::

  $ mpirun -np 5 python `pwd`/demo/helloworld.py

will launch a five-process run of the Python interpreter and run the
test script :file:`demo/helloworld.py` from the source distribution.

You can also run all the *unittest* scripts::

  $ mpiexec -n 5 python test/runtests.py

or, if you have nose_ unit testing framework installed::

  $ mpiexec -n 5 nosetests -w test

.. _nose: http://nose.readthedocs.io/

or, if you have `py.test`_ unit testing framework installed::

  $ mpiexec -n 5 py.test test/

.. _py.test: http://docs.pytest.org/
