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

Equations
==========


.. contents::

Lyapunov Equations
..................

.. currentmodule:: pymess.equations
.. autosummary::
    :nosignatures:

    std.EquationGLyap
    so1.EquationGLyapSO1
    so2.EquationGLyapSO2
    dae1.EquationGLyapDAE1
    dae2.EquationGLyapDAE2


Riccati Equations
...................

.. currentmodule:: pymess.equations
.. autosummary::
    :nosignatures:

    std.EquationGRiccati
    so1.EquationGRiccatiSO1
    so2.EquationGRiccatiSO2
    dae1.EquationGRiccatiDAE1
    dae2.EquationGRiccatiDAE2


First Order Systems
......................

.. currentmodule:: pymess.equations.std
.. autosummary::
    :nosignatures:

    EquationGLyap
    EquationGRiccati

.. autoclass:: EquationGLyap
.. autoclass:: EquationGRiccati

Second Order Systems
.......................

Second Order (SO 1)
--------------------

.. currentmodule:: pymess.equations.so1
.. autosummary::
    :nosignatures:

    EquationGLyapSO1
    EquationGRiccatiSO1

----

.. autoclass:: EquationGLyapSO1
.. autoclass:: EquationGRiccatiSO1


Second Order (SO 2)
--------------------

.. currentmodule:: pymess.equations.so2
.. autosummary::
    :nosignatures:

    EquationGLyapSO2
    EquationGRiccatiSO2

----

.. autoclass:: EquationGLyapSO2
.. autoclass:: EquationGRiccatiSO2


Differential Algebraic Systems
...............................

Index 1 Differential Algebraic System (DAE 1)
---------------------------------------------

.. currentmodule:: pymess.equations.dae1
.. autosummary::
    :nosignatures:

     EquationGLyapDAE1
     EquationGRiccatiDAE1

----

.. autoclass:: EquationGLyapDAE1
.. autoclass:: EquationGRiccatiDAE1


Index 2 Differential Algebraic System (DAE 2)
---------------------------------------------

.. currentmodule:: pymess.equations.dae2
.. autosummary::
    :nosignatures:

    EquationGLyapDAE2
    EquationGRiccatiDAE2

----

.. autoclass:: EquationGLyapDAE2
.. autoclass:: EquationGRiccatiDAE2


