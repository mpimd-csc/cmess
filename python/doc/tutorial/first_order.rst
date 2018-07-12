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

.. _tutorial_firstorder:

First Order Systems LTI Systems
===============================

.. contents::
.. sectionauthor:: |mbehr|


.. currentmodule: pymess.equations


A more advanced way to solve large scale Lyapunov and Riccati equation
is to use the  `equation` objects from |pymess|.
The idea behind the `equation` objects is, that the Lyapunov and
Riccati equation arises in general from a linear dynamical system of the form

.. math::

        \mathcal{E} \dot{x} (t) &= \mathcal{A} x  + \mathcal{B} u        \\
                    y(t)        &= \mathcal{C} x


Such systems are also called linear time invariant systems (LTI).
:math:`x` is in general referred as the state,
:math:`u` is the control and :math:`y` is the output.
|pymess| as well as |cmess| focus on the case where
:math:`\mathcal{A}` and :math:`\mathcal{E}` are large and sparse.
This is for example the case, when the matrices come from a semidiscretization of
a PDE in space. The control to state matrix :math:`\mathcal{B}`
has often in that case only a few columns. So the number of
control is small compared to the dimension of the state space.
The state to output matrix :math:`\mathcal{C}`
has also only a few rows, so the number of "measured" outputs
from the state :math:`x` is also small.
The fact that the number of inputs and outputs is quiet important,
because we can not handle large scale matrix equation with many
inputs/outputs or rows of :math:`\mathcal{B}` / columns of :math:`\mathcal{C}`
efficiently in our algorithms.


Equation Objects
--------------------

In order to exploit the structure of the matrices :math:`\mathcal{A}`
and :math:`\mathcal{E}`, we have `equation` objects.
Consider the case that :math:`\mathcal{E}` or :math:`\mathcal{A}`
have a special block structure which allows efficient multiplication and
solving with these matrices.


The simplest `equation` objects are the

    - :py:class:`~pymess.equations.std.EquationGLyap`
    - :py:class:`~pymess.equations.std.EquationGRiccati`


They store the necassary data for a standard LIT system like above.

We consider more difficult `equation` objects in  :ref:`tutorial_secondorder`, :ref:`tutorial_dae2` and  :ref:`tutorial_dae1`.

It is also possible to define your own equation objects using `callback` functionality,
see here :ref:`tutorial_callback`.


Example First Order System - Lyapunov Equation
----------------------------------------------

.. literalinclude:: ../../example/pymess_lradi.py
    :language: python
    :pyobject: main


Example: :download:`pymess_lradi.py <../../example/pymess_lradi.py>`




Example First Order System - Riccati Equation
---------------------------------------------

.. literalinclude:: ../../example/pymess_lrnm.py
    :language: python
    :pyobject: main


Example: :download:`pymess_lrnm.py <../../example/pymess_lrnm.py>`



































































