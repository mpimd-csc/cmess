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

r"""Interface to mess_equation_glyap_dae2 and mess_equation_griccati_dae2."""
from pymess.equation import Equation

#
class EquationGLyapDAE2(Equation):
    r"""Interface to mess_equation_glyap_dae2.

    Represents a generalized Lyapunov equation for a Index 2 system

    .. math::
        :nowrap:

        \begin{eqnarray*}
            \begin{bmatrix} M & 0  \\ 0  & 0 \end{bmatrix}
            \begin{bmatrix} \dot{z} \\ 0 \end{bmatrix}
            =
            \begin{bmatrix} A & G  \\ G^T & 0  \end{bmatrix}
            \begin{bmatrix}  z \\ p \end{bmatrix}
            +
            \begin{bmatrix} B \\ 0 \end{bmatrix} u.
        \end{eqnarray*}

    Preliminaries:
        * :math:`n_p<n_v` and :math:`\delta \in \mathbb{R}` and :math:`\delta <0`
        * :math:`M \in \mathbb{R}^{n_v \times n_v}` symmetric and positive definite
        * :math:`G \in \mathbb{R}^{n_v \times n_p}` has full rank
        * :math:`B \in \mathbb{R}^{n_v \times n_p}`


    We solve the following Lyapunov equations:

    .. math::
        :nowrap:

        \begin{eqnarray*}
        (\Pi A) X  M^T  +  M  X  (\Pi A)^T &=& -(\Pi B)(\Pi B)^T          \\
                (\Pi A)^T X  M  +  M^T  X  (\Pi A) &=& -(C\Pi^T )^T(C\Pi^T )
        \end{eqnarray*}

    Parameters
    ----------
    opt: options instance
        :class:`~pymess.options.Options` instance

    m: (nv,nv) sparse real symmetric positive definite matrix
        Matrix :math:`M` of Index 2 system

    a: (nv,nv) sparse real matrix
        Matrix :math:`A` of Index 2 system

    g: (nv,np) sparse real full rank matrix
        Matrix :math:`G` of Index 2 system.

    b: (nv,p)  or (p,nv) array_like real matrix
        Right Hand Side of Lyapunov equation.

    delta: real and negativ, standard value: -0.02

    k0: (p,nv)  or (nv,p) array_like real matrix or None
        Initial Feedback for Stabilized Matrix Pencil

    Returns
    -------

    Raises
    ------

    See Also
    --------
    :class:`.EquationGRiccatiDAE2`
    :func:`pymess.lradi.lradi`

    Notes
    -----
    Supported types are ndarray of numpy and csc, coo or csr matrix of scipy.

    References
    ----------
    :cite:`BenSSetal13`, :cite:`BaeBSetal15`

    Examples
    --------
    Further examples :download:`pymess_lradi_dae2.py <../../example/pymess_lradi_dae2.py>`

    """

    def __init__(self, opt, m, a, g, b, delta, k0):

        super(EquationGLyapDAE2, self).__init__("EquationGLyapDAE2", opt, m.shape[0])
        self.m = m
        self.a = a
        self.g = g
        self.b = b
        self.delta = delta
        self.k0 = k0

    def __delattr__(self, key):
        raise TypeError(
            "Trying to delete attributes of %r is not allowed" % self)


class EquationGRiccatiDAE2(Equation):
    r"""Interface to mess_equation_griccati_dae2.

    Represents a generalized Riccati equation for a Index 2 system

    .. math::
        :nowrap:

        \begin{eqnarray*}
            \begin{bmatrix} M & 0  \\ 0  & 0 \end{bmatrix}
            \begin{bmatrix} \dot{z} \\ 0 \end{bmatrix}
            &=&
            \begin{bmatrix} A & G  \\ G^T & 0  \end{bmatrix}
            \begin{bmatrix}  z \\ p \end{bmatrix}
            +
            \begin{bmatrix} B \\ 0 \end{bmatrix} u \\
                    y &=& C z
        \end{eqnarray*}

    Preliminaries:
        * :math:`n_p<n_v` and :math:`\delta \in \mathbb{R}` and :math: `\delta <0`
        * :math:`M \in \mathbb{R}^{n_v \times n_v}` symmetric and positive definite
        * :math:`G \in \mathbb{R}^{n_v \times n_p}` has full rank
        * :math:`B \in \mathbb{R}^{n_v \times n_p}`
        * :math:`C \in \mathbb{R}^{q \times n_v}`

    We solve the following Riccati equations:

    with the default `options`:

    .. math::
        :nowrap:

        \begin{eqnarray*}
            \Pi A X M^T + M  X A^T \Pi^ T - M^T X C^T C X M^T &=& - \Pi B B^T \Pi^T       \\
                    \Pi^T Z &=& Z
        \end{eqnarray*}

    or, with `options.type = pymess.MESS_OP_TRANSPOSE`,

    .. math::
        :nowrap:

        \begin{eqnarray*}
            \Pi A^T X M  + M^T X A \Pi^T - M^T X B B^T X M  &=& - \Pi C^T C \Pi^T       \\
                    \Pi^T Z &=& Z
        \end{eqnarray*}


    Parameters
    ----------
    opt: options instance
        :class:`~pymess.options.Options` instance

    m: (nv,nv) sparse real symmetric positive definite matrix
        Matrix :math:`M` of Index 2 system

    a: (nv,nv) sparse real matrix
        Matrix :math:`A` of Index 2 system

    g: (nv,np) sparse real full rank matrix
        Matrix :math:`G` of Index 2 system.

    b: (nv,q)  or (q,nv) array_like real matrix
        Matrix :math:`B` of Index 2 system.

    c: (p,nv)  or (nv,p) array_like real matrix
        Matrix :math:`C` of Index 2 system.


    delta: real and negativ, common choice -0.02


    Raises
    ------

    See Also
    --------
    :class:`.EquationGLyapDAE2`
    :func:`pymess.lrnm.lrnm`

    Notes
    -----
    Supported types are ndarray of numpy and csc, coo or csr matrix of scipy.

    References
    ----------


    Examples
    --------
    Further examples :download:`pymess_lrnm_dae2.py <../../example/pymess_lrnm_dae2.py>`

    """

    def __init__(self, opt, m, a, g, b, c, delta):

        super(EquationGRiccatiDAE2, self).__init__("EquationGRiccatiDAE2", opt, m.shape[0])
        self.m = m
        self.a = a
        self.g = g
        self.b = b
        self.c = c
        self.delta = delta

    def __delattr__(self, key):
        raise TypeError(
            "Trying to delete attributes of %r is not allowed" % self)
