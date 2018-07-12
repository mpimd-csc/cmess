#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>.
#
# Copyright (C) Peter Benner, Martin Koehler, Jens Saak and others
#               2009-2018
#

r"""Interface to the mess_glyap, mess_gstein function."""
from ._c_interface import _glyap, _gstein
from .enum import MESS_OP_TRANSPOSE


def glyap(a, e, y, op=MESS_OP_TRANSPOSE, ahat=None, ehat=None, qa=None, qe=None):
    r"""Interface to the mess_glyap function.

    This functions solves the Lyapunov equations

    .. math::

        A^T X       +     X A       + Y     &= 0          \\
        A^T X E     + E^T X A       + Y     &= 0          \\
        A   X       +     X A^T     + Y     &= 0          \\
        A   X E^T   + E   X A^T     + Y     &= 0


   The eigenvalues of :math:`A` / :math:`(A, E)` must fullfill the condition
   :math:`\lambda_i + \lambda_j \neq 0`.

   In contrast to :func:`pymess.easy.lyap` no stability is requiered
   and the only condition on the right hand side :math:`Y` is symmetry.
   Please not that :func:`pymess.glyap3.glyap` does not check the symmetry of :math:`Y`.

   :func:`pymess.glyap3.glyap` uses internally a BLAS level 3 block algorithm.
   :func:`pymess.glyap3.glyap` computes the real Schur decomposition of :math:`A` / :math:`(A, E)`
   in order to solve the equation.
   The factorization is returned too.
   You can reuse the real Schur decomposition in subsequent calls of :func:`pymess.glyap3.glyap`.
   This is usefull if you want to solve several standard / generalized Lyapunov
   equation with different right hand sides.

   The real Schur decomposition of :math:`A` fullfills: :math:`A = Q_A \hat{A} Q_A^T`.

   The generalized real Schur decomposition of :math:`(A,E)` fullfills:
   :math:`A = Q_A \hat{A} Q_E^T` and :math:`E = Q_A \hat{E} Q_E^T`.

   If you call :func:`pymess.glyap3.glyap` with the real Schur decomposition
   only the solution :math:`X` is returned.

   If you do not provide the Schur decomposition then the Schur decomposition
   is returned.

    Parameters
    ----------
    a: (n,n) array_like real matrix
        Matrix :math:`A` of Lyapunov equation.

    e: (n,n) array_like real matrix, optional
        Matrix :math:`E` of Lyapunov equation.
        If e is None, then :math:`E` of Lyapunov equation will be the identity matrix.

    y: (n,n)  array_like real symmetric matrix
        The matrix :math:`Y` of Lyapunov equation.

    op: int 0 or 1
        0 refers to the equation where the first matrix of the equation is not transposed.
        1 refers to the equation where the first matrix of the equation is transposed.
        Use :data:`pymess.enum.MESS_OP_NONE` for 0  or
        :data:`pymess.enum.MESS_OP_TRANSPOSE` for 1.

    ahat: (n,n)  array_like real matrix, optional
        The matrix :math:`\hat{A}` of real Schur decomposition equation.

    ehat: (n,n)  array_like real matrix, optional
        The matrix :math:`\hat{E}` of real Schur decomposition equation.

    qa: (n,n)  array_like real matrix, optional
        The matrix :math:`Q_A` of real Schur decomposition equation.

    qe: (n,n)  array_like real matrix, optional
        The matrix :math:`Q_E` of real Schur decomposition equation.


    Returns
    -------
    x: (n,n) array_like real matrix
        Solution :math:`X` of Lyapunov equation.

    ahat: (n,n)  array_like real matrix, optional
        The matrix :math:`\hat{A}` of real Schur decomposition equation.

    ehat: (n,n)  array_like real matrix, optional
        The matrix :math:`\hat{E}` of real Schur decomposition equation.

    qa: (n,n)  array_like real matrix, optional
        The matrix :math:`Q_A` of real Schur decomposition equation.

    qe: (n,n)  array_like real matrix, optional
        The matrix :math:`Q_E` of real Schur decomposition equation.

    Raises
    ------

    See Also
    --------

    Notes
    -----
    Supported types are ndarray and matrix of numpy.

    References
    ----------

    Examples
    --------

    Further examples :download:`pymess_glyap.py <../../example/pymess_glyap.py>`

    """

    return _glyap(a, e, y, op, ahat, ehat, qa, qe)

def gstein(a, e, y, op=MESS_OP_TRANSPOSE, ahat=None, ehat=None, qa=None, qe=None):
    r"""Interface to the mess_gstein function.

    This functions solves the Stein equations

    .. math::

        A^T X A     -     X         + Y     &= 0          \\
        A^T X A     - E^T X E       + Y     &= 0          \\
        A   X A^T   -     X         + Y     &= 0          \\
        A   X A^T   + E   X E^T     + Y     &= 0


   The eigenvalues of :math:`A` / :math:`(A, E)` must fullfill the condition
   :math:`\lambda_i  \lambda_j \neq 1`.

   Please not that :func:`pymess.glyap3.gstein` does not check the symmetry of :math:`Y`.

   :func:`pymess.glyap3.gstein` uses internally a BLAS level 3 block algorithm.
   :func:`pymess.glyap3.gstein` computes the real Schur decomposition of :math:`A` / :math:`(A, E)` in order
   to solve the equation.
   The factorization is returned too.
   You can reuse the real Schur decomposition in subsequent calls of :func:`pymess.glyap3.gstein`.
   This is usefull if you want to solve several standard / generalized Stein
   equation with different right hand sides.

   The real Schur decomposition of :math:`A` fullfills: :math:`A = Q_A \hat{A} Q_A^T`.

   The generalized real Schur decomposition of :math:`(A,E)` fullfills:
   :math:`A = Q_A \hat{A} Q_E^T` and :math:`E = Q_A \hat{E} Q_E^T`.

   If you call :func:`pymess.glyap3.gstein` with the real Schur decomposition
   only the solution :math:`X` is returned.

   If you do not provide the Schur decomposition then the Schur decomposition
   is returned.

    Parameters
    ----------
    a: (n,n) array_like real matrix
        Matrix :math:`A` of Stein equation.

    e: (n,n) array_like real matrix, optional
        Matrix :math:`E` of Stein equation.
        If e is None, then :math:`E` of Stein equation will be the identity matrix.

    y: (n,n)  array_like real symmetric matrix
        The matrix :math:`Y` of Stein equation.

    op: int 0 or 1
        0 refers to the equation where the first matrix of the equation is not transposed.
        1 refers to the equation where the first matrix of the equation is transposed.
        Use :data:`pymess.enum.MESS_OP_NONE` for 0  or
        :data:`pymess.enum.MESS_OP_TRANSPOSE` for 1.

    ahat: (n,n)  array_like real matrix, optional
        The matrix :math:`\hat{A}` of real Schur decomposition equation.

    ehat: (n,n)  array_like real matrix, optional
        The matrix :math:`\hat{E}` of real Schur decomposition equation.

    qa: (n,n)  array_like real matrix, optional
        The matrix :math:`Q_A` of real Schur decomposition equation.

    qe: (n,n)  array_like real matrix, optional
        The matrix :math:`Q_E` of real Schur decomposition equation.

    Returns
    -------
    x: (n,n) array_like real matrix
        Solution :math:`X` of Stein equation.

    ahat: (n,n)  array_like real matrix, optional
        The matrix :math:`\hat{A}` of real Schur decomposition equation.

    ehat: (n,n)  array_like real matrix, optional
        The matrix :math:`\hat{E}` of real Schur decomposition equation.

    qa: (n,n)  array_like real matrix, optional
        The matrix :math:`Q_A` of real Schur decomposition equation.

    qe: (n,n)  array_like real matrix, optional
        The matrix :math:`Q_E` of real Schur decomposition equation.

    Raises
    ------

    See Also
    --------

    Notes
    -----
    Supported types are ndarray and matrix of numpy.

    References
    ----------

    Examples
    --------

    Further examples :download:`pymess_gstein.py <../../example/pymess_gstein.py>`

    """

    return _gstein(a, e, y, op, ahat, ehat, qa, qe)
