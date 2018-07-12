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

r"""Interface to mess_dense_nm_gmpare."""
from pymess._c_interface import _dense_nm_gmpare
from .enum import MESS_OP_TRANSPOSE, MESS_2_NORM

# pylint: disable=line-too-long
def dense_nm_gmpare(x0, a, e, q, g, plus=False, linesearch=False, trans=MESS_OP_TRANSPOSE, \
                    maxit=50, absres_tol=1e-11, relres_tol=1e-12, nrm=MESS_2_NORM):
    r"""This class represents the interface to the mess_dense_nm_gmpare function.

    This functions solves several instances of Riccati Equations for dense matrices:

    .. math::

        A^T X       +     X A     + / -      X G X    + Q  &= 0         \\
        A^T X E     + E^T X A     + / -  E^T X G X E  + Q  &= 0         \\

        A   X E^T   + E   X A^T   + / -  E   X G X E^T  + Q  &= 0       \\
        A   X E     + E^T X A     + / -  E^T X G X E    + Q  &= 0       \\



    Parameters
    ----------
    x0: (n,n) dense symmetric matrix (optional)
        Initial guess matrix :math:`X_0` for Newton Iteration.

    a: (n,n) dense real matrix
        Matrix :math:`A` of Riccati Equation.

    e: (n,n) dense real matrix (optional)
        Matrix :math:`E` of Riccati Equation.

    q: (n,n) dense real symmetric matrix
        Matrix :math:`Q` of Riccati Equation.

    g: (n,n) dense real symmetric matrix
        Matrix :math:`G` of Riccati Equation.

    plus: int
        Determines whether quadratic term :math:`XGX` should be added or subtracted.

    linesearch: int
        Determines whether linesearch should be used in computation.

    trans: int
        Determines whether :math:`A^TX E +  E^T XA +/- E^T XGX E + Q = 0` or
        :math:`AXE^T + EXA^T +/- EXGXE^T + Q = 0` should be computed.
        Use :data:`pymess.enum.MESS_OP_NONE` for 0  or
        :data:`pymess.enum.MESS_OP_TRANSPOSE` for 1.

    maxit: int
        Maximum number of Newton iterations.

    absres_tol: real and positive, tolerance for absolute residual

    relres_tol: real and positive, tolerance for relative residual

    nrm: int
        Determines which norm is used for residual computations
        Use :data:`pymess.enum.MESS_2_NORM`, :data:`pymess.enum.MESS_FROBENIUS_NORM`, :data:`pymess.enum.MESS_1_NORM`, or :data:`pymess.enum.MESS_INF_NORM`

    Returns
    -------
    x: (n,n) symmetric dense solution matrix
        :math:`X`.

    absres: real and nonnegative
        achieved absolute residual

    relres: real and nonnegative
        achieved relative residual

    Raises
    ------

    See Also
    ---------
    :func:`pymess.easy.care`

    Notes
    -----

    References
    ----------
    :cite:`Ben97b`, :cite:`BenB98`, :cite:`JonV06`

    Examples
    ---------

    Further examples :download:`pymess_gmpare.py <../../example/pymess_gmpare.py>`



    """

    return _dense_nm_gmpare(x0, a, e, q, g, plus, linesearch, trans, maxit, absres_tol, relres_tol, nrm)
