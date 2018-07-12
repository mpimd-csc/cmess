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

.. _tutorial_secondorder:


Second Order Systems
====================

.. contents::
.. sectionauthor:: |mbehr|


.. currentmodule: pymess.equations


We transform Second Order Systems to First Order System
using :class:`~pymess.equations.so1.EquationGLyapSO1` and :class:`~pymess.equations.so1.EquationGRiccatiSO1`.
As the transformation is not unique there is also a second
variant implemented
:class:`~pymess.equations.so2.EquationGLyapSO2` and :class:`~pymess.equations.so2.EquationGRiccatiSO2`.

The remaining steps are the same.

Example - Lyapunov Equation
---------------------------

.. literalinclude:: ../../example/pymess_lradi_so1.py
    :language: python
    :pyobject: main


Example: :download:`pymess_lradi_so1.py <../../example/pymess_lradi_so1.py>`

Example: :download:`pymess_lradi_so2.py <../../example/pymess_lradi_so2.py>`


Example - Riccati Equation
--------------------------

.. literalinclude:: ../../example/pymess_lrnm_so1.py
    :language: python
    :pyobject: main


Example: :download:`pymess_lrnm_so1.py <../../example/pymess_lrnm_so1.py>`

Example: :download:`pymess_lrnm_so2.py <../../example/pymess_lrnm_so2.py>`


