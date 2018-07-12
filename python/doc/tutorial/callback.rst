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



.. _tutorial_callback:

Callback Interface
=====================

.. contents::
.. sectionauthor:: |mbehr|

|pymess| allow also define own equations.
You have to inherit from :class:`~pymess.equation.Equation`.
It is necessary that you have to implement the methods:

 * :meth:`~pymess.equation.Equation.ax_apply`
 * :meth:`~pymess.equation.Equation.ex_apply`
 * :meth:`~pymess.equation.Equation.ainv_apply`
 * :meth:`~pymess.equation.Equation.einv_apply`
 * :meth:`~pymess.equation.Equation.apex_apply`
 * :meth:`~pymess.equation.Equation.apeinv_apply`
 * :meth:`~pymess.equation.Equation.parameter`


Depending on your matrices :math:`A` and :math:`E` you can implement
the necessary methods above in a efficient way.
Please not that the input arguments and return values are in each call copied from |cmess| to
|pymess| and vice versa.
Internally the defined methods are wrapped around :code:`C` functions and we setup a
equation structure in :code:`C` and call the implemented methods from :code:`Python`.

Here you can see a basic example:


.. literalinclude:: ../../example/pymess_callback.py
    :language: python
    :pyobject: MyEquation


Now we can use our derived class to solve a algebraic Riccati Equation:

.. literalinclude:: ../../example/pymess_callback.py
    :language: python
    :pyobject: main





Example: :download:`pymess_callback.py <../../example/pymess_callback.py>`


