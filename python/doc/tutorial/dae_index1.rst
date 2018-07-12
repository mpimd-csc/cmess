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

.. include:: ../global.rst

.. _tutorial_dae1:

Index-1 Differential Algebraic Systems
======================================

.. contents::
.. sectionauthor:: |mbehr|

Index-1 Systems can be transformed to a standard or generalized
LTI System using a projection.
You can use  :class:`~pymess.equations.dae1.EquationGLyapDAE1` and :class:`~pymess.equations.dae1.EquationGRiccatiDAE1`
classes to create instances of a Lyapunov or Riccati Equation of a Index-1 DAE System.

Example - Lyapunov Equation
---------------------------

.. literalinclude:: ../../example/pymess_lradi_dae1.py
    :language: python
    :pyobject: main


Example: :download:`pymess_lradi_dae1.py <../../example/pymess_lradi_dae1.py>`


Example - Riccati Equation
--------------------------

.. literalinclude:: ../../example/pymess_lrnm_dae1.py
    :language: python
    :pyobject: main


Example: :download:`pymess_lrnm_dae1.py <../../example/pymess_lrnm_dae1.py>`



