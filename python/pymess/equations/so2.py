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

r"""Interface to mess_equation_glyap_so2 and mess_equation_griccati_so2."""
from pymess.equation import Equation


class EquationGLyapSO2(Equation):
    r"""Interface to mess_equation_glyap_so2.

    Represents a standard or generalized Lyapunov equation for a second order system

    .. math::

        M\ddot{x} + D \dot{x} + K x = B u,

    where :math:`M`, :math:`D`, :math:`K` are real matrices of dimension :math:`n`.

    The system is implicitly transformed to its first order representation:

    .. math::

        \underbrace{\begin{bmatrix} I & 0 \\ 0 & M \end{bmatrix}}_{\mathcal{E}}
        \begin{bmatrix} \dot{x} \\ \ddot{x} \end{bmatrix}
        =
        \underbrace{\begin{bmatrix} 0 & I \\ -K & -D  \end{bmatrix}}_{\mathcal{A}}
        \begin{bmatrix} x \\ \dot{x} \end{bmatrix}
        +
        \underbrace{\begin{bmatrix} 0 \\ B \end{bmatrix}}_{\mathcal{B}}

    Therefore an instance of :class:`.EquationGLyapSO2` represents
    the generalized Lyapunov equation

    .. math::

            \mathcal{A} X \mathcal{E}^T + \mathcal{E} X\mathcal{A}^T + \mathcal{B}\mathcal{B}^T &= 0 \\
            \mathcal{E}^T X \mathcal{A} + \mathcal{A}^T X\mathcal{E} + \mathcal{B}^T\mathcal{B} &= 0

    The linearization of the second order system leads to quadratic shifts.
    In order to avoid numerical problems one can define a lower and upper bounds
    for the absolute value of the shift parameter :math:`p` used in the
    :func:`pymess.lradi.lradi` method. Shifts that not fullfill the inequality
    :math:`lowerbound < |p|<upperbound` are automatically sorted out.
    Common choices are :math:`lowerbound=10^{-5}` and :math:`upperbound=10^5`.

    Parameters
    ----------
    opt: options instance
        :class:`~pymess.options.Options` instance

    m: (n,n) sparse real matrix
        Matrix :math:`M` of second order system.

    d: (n,n) sparse real matrix
        Matrix :math:`D` of second order system.

    k: (n,n) sparse real matrix
        Matrix :math:`K` of second order system.


    b: (2n,p) array_like real matrix or (p,2n) array_like real matrix
        The right hand side of Lyapunov equation :math:`\mathcal{B}`
        We assume that :math:`p \ll n`.

    lowerbound: positive real
        Lower bound for shift parameter criterion.

    upperbound: positive real
        Upper bound for shift parameter criterion.


    Raises
    ------


    See Also
    --------
    :class:`pymess.equations.so1.EquationGLyapSO1`
    :func:`pymess.lradi.lradi`


    Notes
    -----
    Supported types are ndarray of numpy and csc, coo or csr matrix of scipy.

    References
    ----------

    Examples
    --------
    Further examples :download:`pymess_lradi_so2.py <../../example/pymess_lradi_so2.py>`

    """

    def __init__(self, opt, m, d, k, b, lowerbound, upperbound):

        super(EquationGLyapSO2, self).__init__("EquationGLyapSO2", opt, m.shape[0])
        self.m = m
        self.d = d
        self.k = k
        self.b = b
        self.lowerbound = lowerbound
        self.upperbound = upperbound

    def __delattr__(self, key):
        raise TypeError(
            "Trying to delete attributes of %r is not allowed" % self)


class EquationGRiccatiSO2(Equation):
    r"""Interface to mess_equation_griccati_so2.

    Represents a standard or generalized Riccati equation for a second order system

    .. math::
            M\ddot{x} + D \dot{x} + K x &= B u,    \\
                y &= Cx

    where :math:`M`, :math:`D`, :math:`K` are real matrices of dimension :math:`n`.

    The system is implicitly transformed to its first order representation:

    .. math::

        \underbrace{\begin{bmatrix} I & 0 \\ 0 & M \end{bmatrix}}_{\mathcal{E}}
        \begin{bmatrix} \dot{x} \\ \ddot{x} \end{bmatrix}
        &=
        \underbrace{\begin{bmatrix} 0 & I \\ -K & -D  \end{bmatrix}}_{\mathcal{A}}
        \begin{bmatrix} x \\ \dot{x} \end{bmatrix}
        +
        \underbrace{\begin{bmatrix} 0 \\ B \end{bmatrix}}_{\mathcal{B}}            \\
                y &= \underbrace{C}_{\mathcal{C}}x

    The Riccati solver  will then solve

    .. math::

        \mathcal{A}^T X \mathcal{E} + \mathcal{E}^T X \mathcal{A} -
        \mathcal{E}X\mathcal{B}\mathcal{B}^T X \mathcal{E}^T +
        \mathcal{C}^T \mathcal{C}=0

    or

    .. math::

        \mathcal{A}X \mathcal{E}^T + \mathcal{E} X\mathcal{A}^T -
        \mathcal{E}X \mathcal{C}^T \mathcal{C}X \mathcal{E}^T  +
        \mathcal{B}\mathcal{B}^T = 0


    The linearization of the second order system leads to quadratic shifts.
    In order to avoid numerical problems one can define a lower and upper bounds
    for the absolute value of the shift parameter :math:`p` used in the
    :func:`pymess.lradi.lradi` method. Shifts that not fullfill the inequality
    :math:`lowerbound < |p|<upperbound` are automatically sorted out.
    Common choices are :math:`lowerbound=10^{-5}` and :math:`upperbound=10^5`.

    Parameters
    ----------
    opt: options instance
        :class:`~pymess.options.Options` instance

    M: (n,n) sparse real matrix
        Matrix :math:`M` of second order system.

    D: (n,n) sparse real matrix
        Matrix :math:`D` of second order system.

    K: (n,n) sparse real matrix
        Matrix :math:`K` of second order system.

    B: (n,p) or (2n,p) array_like real matrix
        Matrix :math:`\mathcal{B}` of Riccati equation.

    C: (2q,n) or (q,n) array_like real matrix
        Matrix :math:`\mathcal{C}` of Riccati equation.

    lowerbound: positive real
        Lower bound for shift parameter criterion.

    upperbound: positive real
        Upper bound for shift parameter criterion.


    Raises
    ------


    See Also
    --------
    :class:`pymess.equations.so2.EquationGRiccatiSO2`
    :func:`pymess.lrnm.lrnm`


    Notes
    -----
    Supported types are ndarray of numpy and csc, coo or csr matrix of scipy.

    References
    ----------


    Examples
    --------
    Further examples :download:`pymess_lrnm_so2.py <../../example/pymess_lrnm_so2.py>`

    """

    def __init__(self, opt, m, d, k, b, c, lowerbound, upperbound):

        super(EquationGRiccatiSO2, self).__init__("EquationGRiccatiSO2", opt, m.shape[0])
        self.m = m
        self.d = d
        self.k = k
        self.b = b
        self.c = c
        self.lowerbound = lowerbound
        self.upperbound = upperbound

    def __delattr__(self, key):
        raise TypeError(
            "Trying to delete attributes of %r is not allowed" % self)
