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

Py-M.E.S.S.
===========

:Release: |release|
:Date: |today|


This is the documentation of |pymess| :cite:`BarKPetal13` the |python| interface for the |cmess| :cite:`messweb`  library.


Features
---------

* Interface to the |cmess| library.
* Solving large scale algebraic (generalized) Lyapunov equations:

.. math::

        A^T X   +     X A &= -C^T C,  \\
        A^T XE  + E^T X A &= -C^T C.

* Solving large scale algebraic (generalized) Riccati equations:

.. math::

        A^T X   +     X A   -    XBB^T X  &= -C^T C, \\
        A^T XE  + E^T X A   - E^TXBB^T XE &= -C^T C.

* |blas| Level 3 solver for (generalized) Lyapunov / Stein equations:

.. math::

        A^T X   +     X A &= -Y,  \\
        A^T XE  + E^T X A &= -Y.

.. math::

        A^T X A  -      X   &= -Y,  \\
        A^T X A  -  E^T X E &= -Y.

* Solving dense positive algebraic (generalized) Riccati equations:

.. math::

        A^T X   +     X A   +    XG X  + Q &= 0, \\
        A^T XE  + E^T X A   + E^TXGXE  + Q &= 0.


.. currentmodule:: pymess.easy

* Algebraic Lyapunov Equation Solver: :func:`lyap`.

.. currentmodule:: pymess.glyap3

* Algebraic |blas| Level 3 Lyapunov Equation Solver: :func:`glyap`.

* Algebraic |blas| Level 3 Stein Equation Solver: :func:`gstein`.

.. currentmodule:: pymess.easy

* Algebraic Sylvester Equation Solver: :func:`sylvester_sparsedense`.

.. currentmodule:: pymess.easy

* Algebraic Riccati Equation Solver: :func:`care`.

.. currentmodule:: pymess.dense_nm_gmpare

* Dense Newton Positive Algebraic Riccati Equation Solver: :func:`dense_nm_gmpare`.

.. currentmodule:: pymess.lradi

* Low Rank ADI Method: :func:`lradi`.

.. currentmodule:: pymess.lrnm

* Low Rank Newton Method: :func:`lrnm`.

.. currentmodule:: pymess.equations.std

* First Order Systems:

    * :class:`EquationGLyap`,
    * :class:`EquationGRiccati`.

.. currentmodule:: pymess.equations.so1

* Second Order Systems (SO1):

    * :class:`EquationGLyapSO1`,
    * :class:`EquationGRiccatiSO1`.

.. currentmodule:: pymess.equations.so2

* Second Order Systems (SO2):

    * :class:`EquationGLyapSO2`,
    * :class:`EquationGRiccatiSO2`.

.. currentmodule:: pymess.equations.dae1

* Index-1 Differential Algebraic Systems:

    * :class:`EquationGLyapDAE1`,
    * :class:`EquationGRiccatiDAE1`.

.. currentmodule:: pymess.equations.dae2

* Index-2 Differential Algebraic Systems:

    * :class:`EquationGLyapDAE2`,
    * :class:`EquationGRiccatiDAE2`.

* Callback functionality:

    * :ref:`Unittest framework <unittest>`.



Topics
--------

.. toctree::
    :maxdepth: 2

    Installation <install>
    Tutorial <tutorial/index>
    unittest
    release/release


API Reference
--------------

The API of all available functions and classes.

.. toctree::
    :maxdepth: 3

    api/equations
    api/options
    api/status
    api/methods
    api/enum
    api/equation
    api/exceptions


Indices and Tables
-------------------

* :ref:`genindex`
* :ref:`search`


Copyright
---------

Copyright 2009 - 2018 by

* |bennerp|
* |koehlerm|
* |saakj|

The software is licensed under GPLv2 or later. See :download:`COPYING <../../COPYING>` for details.


Contact
--------
More Information about the |mess| project can be found under  http://www.mpi-magdeburg.mpg.de/projects/mess.


References
----------

.. bibliography:: ../../documents/csc-bibfiles/csc.bib
    :labelprefix: A
    :style: plain
.. bibliography:: ../../documents/csc-bibfiles/mor.bib
    :labelprefix: B
    :style: plain
.. bibliography:: ../../documents/csc-bibfiles/software.bib
    :labelprefix: C
    :style: plain
