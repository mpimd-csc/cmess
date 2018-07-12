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


Introduction
============

.. contents::



Py-M.E.S.S. modules
--------------------

After successful installation of |pymess| you can import the `pymess` module.
Let us call the :py:func:`~pymess.misc.mess_version`  function it will print some information about your
|cmess| build.


.. doctest::

    >>> from pymess import *
    >>> mess_version()


You should obtain a similiar like this:

.. code-block:: python

    MESS - The Matrix Equations Sparse Solver Library  [ 32 bit integer, UMFPACK, AMD, COLAMD, CSPARSE, ARPACK, X11_XPM, SSE2, ]


Now we call the :py:func:`~pymess.misc.eps`  function, which is actually an interface to `mess_eps`.
It returns the machine epsilon:

.. doctest::

    >>> from pymess import *
    >>> eps()   # doctest: +SKIP


This was the first part of the tutorial, go to :ref:`tutorial_easy` and learn how to use  :py:func:`~pymess.easy.lyap`
and  :py:func:`~pymess.easy.care` to solve algebraic Lyapunov and Riccati equations.



