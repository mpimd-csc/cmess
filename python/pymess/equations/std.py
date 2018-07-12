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

r"""Interface to mess_equation_lyap, mess_equation_glyap, mess_equation_riccati, mess_equation_griccati."""
from pymess.equation import Equation


class EquationGLyap(Equation):
    r"""Interface to mess_equation_lyap and mess_equation_glyap.

    Represents a standard or generalized Lyapunov equation

    .. math::

        A^T X       + X A           &= -B^T B,          \\
        A X         + X A^T         &= -C C^T,          \\
        A^T X E     + E^T X A       &= -B^T B,          \\
        A X E^T     + E X A^T       &= -C C^T.


    arising from a first order system

    .. math::

            E\dot{x} &= Ax + B u,\\
                y    &= C x.

    Parameters
    ----------
    opt: options instance
        :class:`~pymess.options.Options` instance

    a: (n,n) sparse real matrix
        Matrix :math:`A` of Lyapunov equation.

    e: (n,n) sparse real matrix or None
        Matrix :math:`E` of Lyapunov equation.
        If e is None, then :math:`E` of Lyapunov equation will be the identity matrix.

    rhs: (n,p) array_like real matrix or (p,n) array_like real matrix
        The right hand side of Lyapunov equation either :math:`B` or :math:`C`.
        We assume that :math:`p \ll n`.

    Raises
    ------


    See Also
    --------
    :class:`.EquationGRiccati`
    :func:`pymess.lradi.lradi`

    Notes
    -----
    Supported types are ndarray of numpy and csc, coo or csr matrix of scipy.

    References
    ----------


    Examples
    --------
    Further examples :download:`pymess_lradi.py <../../example/pymess_lradi.py>`


    """

    def __init__(self, opt, a, e, rhs):
        r"""Create a :class:`.EquationGLyap` instance"""
        shapea = a.shape
        dim = shapea[0]
        super(EquationGLyap, self).__init__("EquationGLyap", opt, dim)
        self.a = a
        self.rhs = rhs
        self.e = e

    def __delattr__(self, key):
        raise TypeError(
            "Trying to delete attributes of %r is not allowed" % self)


class EquationGRiccati(Equation):
    r"""Interface to mess_equation_riccati and mess_equation_griccati.

    Represents a standard or generalized Riccati equation

    .. math::

        A^T X   + X A       - XBB^T X       &= -C^T C,      \\
        A X     + X A^T     - XC^TC X       &= -BB^T,       \\
        A^T XE  + E^T X A   - E^TXBB^T XE   &= -C^T C,      \\
        A XE^T  + E X A^T   - EXC^TC XE^T   &= -BB^T.

    arising from a first order system

    .. math::

            E\dot{x}&= Ax + B u,   \\
                y   &= C x.

    Parameters
    ----------
    opt: options instance
        :class:`~pymess.options.Options` instance

    a: (n,n) sparse real matrix
        Matrix :math:`A` of Riccati equation.

    e: (n,n) sparse real matrix or None
        Matrix :math:`E` of Riccati equation.
        If e is None, then :math:`E` of Riccati equation will be the identity matrix.


    b: (n,p) array_like real matrix
        Matrix :math:`B` of Riccati equation.
        We assume that  :math:`p \ll n`.

    c: (q,n) array_like real matrix
        Matrix :math:`C` of Riccati equation.
        We assume that  :math:`q \ll n`.


    Raises
    ------


    See Also
    --------
    :class:`pymess.equations.std.EquationGLyap`
    :func:`pymess.lrnm.lrnm`

    Notes
    -----
    Supported types are ndarray of numpy and csc, coo or csr matrix of scipy.

    References
    ----------

    Examples
    --------
    Further examples :download:`pymess_lrnm.py <../../example/pymess_lrnm.py>`

    """

    def __init__(self, opt, a, e, b, c):
        r"""Create a :class:`.EquationGRiccati` instance"""
        shapea = a.shape
        dim = shapea[0]
        super(EquationGRiccati, self).__init__("EquationGRiccati", opt, dim)
        self.a = a
        self.b = b
        self.c = c
        self.e = e
        self.rhs = None

    def __delattr__(self, key):
        raise TypeError(
            "Trying to delete attributes of %r is not allowed" % self)
