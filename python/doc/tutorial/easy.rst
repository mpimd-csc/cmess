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

.. _tutorial_easy:


Easy Interface
===============================

.. contents::
.. sectionauthor:: |mbehr|


Solve algebraic Lyapunov equation
---------------------------------

In this section we want to solve the algebraic Lyapunov equation
using |pymess|. The function :py:func:`~pymess.easy.lyap` can solve
the standard Lyapunov Equation :math:`AX + XA^T = -B B^T`
and the generalized Lyapunov Equation :math:`AXE^T + E X A^T = -BB^T`.
:py:func:`~pymess.easy.lyap` automatically switches to a dense/sparse solver if your
matrix :math:`A` and :math:`E` are dense/sparse.
The solution :math:`X` is given in factored form, e.g.
:math:`X\approx ZZ^T`.
The matrices :math:`A` and :math:`E` must be real, quadratic and
of the same size
We assume that the eigenvalues of :math:`A` or :math:`(A,E)` lie in the left open halfplane,
otherwise such a factorization of :math:`X` is not feasible. We also
assume that the number of colums of :math:`B` is small compared to the dimension of the
size of :math:`A` and :math:`E`  for large sparse matrices.
The ADI methods becomes inefficient, if :math:`B` has many columns.
In all cases :math:`B` has to be dense.

After solving the Lyapunov we compute the absolute and relative residual using
:py:func:`~pymess.residual.res2_lyap`.


The first example solve a dense Lyapunov equation:

.. doctest::

    >>> from __future__ import print_function
    >>> from pymess import *
    >>> from numpy.random import rand
    >>> from numpy import eye
    >>> n = 50;
    >>> p = 40;
    >>> A = rand(n,n)-2*n*eye(n);                           # create random matrix A, B
    >>> B = rand(n,p);
    >>> Z = lyap(A, None, B)                                # solve Lyapunov equation
    >>> res,rhs,rel = res2_lyap(A,None, B, Z,MESS_OP_NONE); # compute and print the residual
    >>> print("relative residual =", rel)                   # doctest: +SKIP



The steps for large sparse matrices are exactly the same.

.. literalinclude:: ../../example/pymess_lyap.py
    :language: python
    :pyobject: main


Example: :download:`pymess_lyap.py <../../example/pymess_lyap.py>`



Solve Sylvester equation
------------------------

Here we solve the standard/semi-generalized/generalized Sylvester equation
:math:`A X + X H + M = 0`,
:math:`A X + E X H + M = 0` or
:math:`A X F + E X H + M = 0`.
Therefore we use the function :py:func:`~pymess.easy.sylvester_sparsedense`.
This function can handle all cases, the standard, the semi-generalized and the generalized one.
The residual can be computed using :py:func:`~pymess.residual.res2_sylv_sd`.
The documentation of :py:func:`~pymess.easy.sylvester_sparsedense` contains small examples.


.. literalinclude:: ../../example/pymess_sylvester_sparsedense.py
    :language: python
    :pyobject: main


Example: :download:`pymess_sylvester_sparsedense.py <../../example/pymess_sylvester_sparsedense.py>`



Solve algebraic Riccati equation
---------------------------------

Now let us consider the standard/generalized algebraic Riccati equation
:math:`A^T X + X A  - XBB^TX + C^T C = 0`  or
:math:`A^T X E     + E^TX A    - E^TXBB^TXE    + C^T C = 0`.
In order to find a symmetric positive semidefinite solution
:math:`X=\approx ZZ^T` we use the function :py:func:`~pymess.easy.care`.
This function can handle both cases the standard and the generalized one.
We assume also that all eigenvalues :math:`A` or :math:`(A,E)` lie in the left open halfplane.
The residual can be computed using :py:func:`~pymess.residual.res2_ric`.
The documentation of :py:func:`~pymess.easy.care` contains small examples.


.. literalinclude:: ../../example/pymess_care.py
    :language: python
    :pyobject: main


Example: :download:`pymess_care.py <../../example/pymess_care.py>`


