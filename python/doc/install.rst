..
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   You should have received a copy of the GNU General Public License
   along with this program; if not, see <http://www.gnu.org/licenses/>.
   Copyright (C) Peter Benner, Martin Koehler, Jens Saak and others
                 2009-2018

.. include:: global.rst

.. _installation:


Installation
============

In order to install |pymess| you have to peform the following steps:

    1. :ref:`Compile C-M.E.S.S.  <installation_compile>`
    2. :ref:`Install Py-M.E.S.S. <installation_install>`
    3. :ref:`Start Py-M.E.S.S.   <installation_start>`


.. _installation_compile:

Compile C-M.E.S.S.
------------------------------

Perform an out of source build of |cmess| using |cmake|.

.. code-block:: bash

    mkdir   <MY BUILD>
    cd      <MY BUILD>
    cmake   <PATH TO SOURCE> -DPYTHON=ON [FURTHER CMAKE OPTIONS]

|cmake| will print out a :code:`Configuration Summary`.
|numpy| and |scipy| are necessary for |pymess|.

.. code-block:: bash
    :emphasize-lines: 9
    :linenos:

    -- --------------------- Configuration Summary -------------------
    -- ----- Enabled Features -----
    -- OpenMP support:                OPENMP          = FALSE
    -- SuiteSparse and Components:    SUITESPARSE     = TRUE
    -- SUPERLU support:               SUPERLU         = ON
    -- ARPACK support                 ARPACK          = FALSE
    -- X11 Plot support               X11             = TRUE
    -- MATIO support                  MATIO           = TRUE
    -- Python support (Py-M.E.S.S.)   PYTHON          = ON (2.7;3.5)
    -- Documentation:                 DOC             = ON
    ...

|cmake| reported in line `9` that |pymess| can be build for :code:`Python 2.7` and  :code:`Python 3.5`.

.. _installation_install:

Install Py-M.E.S.S.
-------------------

Now simply start the compilation process and change directory into  `python`.

.. code-block:: bash

    make
    cd python/python_2.7
    python2.7 setup.py build && sudo python2.7 setup.py install

|pymess| is now installed on your machine.
You can also use the custom |cmake| targets:

.. code-block:: bash

    make pymess-2.7-build
    make pymess-2.7-install
    make pymess-2.7-run-examples    # optional, run examples
    make pymess-2.7-run-tests       # optional, run tests

.. _installation_start:

Start Py-M.E.S.S.
-----------------

Try to import the  |pymess| module using:

.. code-block:: bash

    python -c "import pymess"


..
    Problems
    --------

..
    .. topic:: Cannot find `libmess.so` libray

..
    .. code-block:: python

 ..
        ImportError: libmess.so.0: cannot open shared object file: No such file or directory
..
    You have to append the path of ``libmess.so`` dynamic library to the ``LD_LIBRARY_PATH``
    variable
..
    .. code-block:: bash

 ..
        export LD_LIBRARY_PATH=<MY BUILD>/lib:$LD_LIBRARY_PATH
..
    If you use a Mac OS X system, then you have to adapt the ``DYLD_LIBRARY_PATH`` variable

 ..
    .. code-block:: bash

..
        export DYLD_LIBRARY_PATH=<MY BUILD>/lib:$DYLD_LIBRARY_PATH

