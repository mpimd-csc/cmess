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

r"""Interface to mess_equation_glyap_dae1 and mess_equation_griccati_dae1."""
from pymess.equation import Equation

class EquationGLyapDAE1(Equation):
    r"""
    Interface to mess_equation_glyap_dae1.

    Represents a generalized Lyapunov equation for a Index 1 system

    .. math::

            \begin{bmatrix} E_{11} & 0  \\ 0  & 0 \end{bmatrix}
            \begin{bmatrix} \dot{x}_1 \\ \dot{x}_2 \end{bmatrix}
            =
            \begin{bmatrix} A_{11} & A_{12}  \\ A_{21} & A_{22}  \end{bmatrix}
            \begin{bmatrix}  x_1 \\ x_2 \end{bmatrix}
            +
            \begin{bmatrix} B \\ 0 \end{bmatrix} u.

    Preliminaries:

        * :math:`E_{11},A_{11} \in \mathbb{R}^{n_1 \times n_1}`
        * :math:`A_{12} \in \mathbb{R}^{n_1 \times n_2}`
        * :math:`A_{21} \in \mathbb{R}^{n_2 \times n_1}`
        * :math:`A_{22} \in \mathbb{R}^{n_2 \times n_2}` and :math:`A_{22}` is invertible.

    Using that :math:`A_{22}` is invertible leads to the ODE system

    .. math::

        E_{11} \dot{x}_1 =
        (A_{11} - A_{12}A_{22}^{-1}A_{21})x_1 + (B(1:n_1,:)- A_{12}A_{22}^{-1}B(n_1+1:end,:)) u

    Therefore we solve the following Lyapunov equations:

    .. math::

        E_{11} X (A_{11}-A_{12}A_{22}^{-1}A_{21})^T + (A_{11}-A_{12}A_{22}^{-1}A_{21}) X E_{11}^T =
        -\tilde{B} \tilde{B}^T

        \tilde{B} := (B(1:n_1,:) - A_{12}A_{22}^{-1}B(n_1+1:end,:))

        E_{11}^T X (A_{11}-A_{12}A_{22}^{-1}A_{21}) + (A_{11}-A_{12}A_{22}^{-1}A_{21})^T X E_{11} =
        -\tilde{B}^T\tilde{B}

        \tilde{B} := (B(:,1:n_1) - B(:,n_1+1:end) A_{22}^{-1}A_{21})

    Parameters
    ----------
    opt: options instance
        :class:`~pymess.options.Options` instance

    e11: (n1,n1) sparse real matrix
        Matrix :math:`E_{11}` of Index 1 system

    a11: (n1,n1) sparse real matrix
        Matrix :math:`A_{11}` of Index 1 system.

    a12: (n1,n2) sparse real matrix
        Matrix :math:`A_{12}` of Index 1 system.

    a21: (n2,n1) sparse real matrix
        Matrix :math:`A_{21}` of Index 1 system.

    a22: (n2,n2) sparse real matrix
        Matrix :math:`A_{22}` of Index 1 system.

    b: (n1+n2,p) or (p,n1+n2) array_like real matrix
        The right hand side of Lyapunov equation.


    Raises
    ------


    See Also
    --------
    :class:`.EquationGRiccatiDAE1`
    :func:`pymess.lradi.lradi`

    Notes
    -----
    Supported types are ndarray of numpy and csc, coo or csr matrix of scipy.


    References
    ----------


    Examples
    --------
    Further examples :download:`pymess_lradi_dae1.py <../../example/pymess_lradi_dae1.py>`

    """

    def __init__(self, opt, e11, a11, a12, a21, a22, b):

        super(EquationGLyapDAE1, self).__init__("EquationGLyapDAE1", opt, e11.shape[0])
        self.e11 = e11
        self.a11 = a11
        self.a12 = a12
        self.a21 = a21
        self.a22 = a22
        self.b = b

    def __delattr__(self, key):
        raise TypeError(
            "Trying to delete attributes of %r is not allowed" % self)


class EquationGRiccatiDAE1(Equation):
    r"""Interface to mess_equation_griccati_dae1.

    Represents a standard or generalized Riccati equation for a Index 1 system

    .. math::

            \begin{bmatrix} E_{11} & 0  \\ 0  & 0 \end{bmatrix}
            \begin{bmatrix} \dot{x}_1 \\ \dot{x}_2 \end{bmatrix}
            &=
            \begin{bmatrix} A_{11} & A_{12}  \\ A_{21} & A_{22}  \end{bmatrix}
            \begin{bmatrix}  x_1 \\ x_2 \end{bmatrix}
            +
            \begin{bmatrix} B \\ 0 \end{bmatrix} u          \\
                    y &= C\begin{bmatrix} x_1 & x_2 \end{bmatrix}

    Preliminaries:

        * :math:`E_{11},A_{11} \in \mathbb{R}^{n_1 \times n_1}`
        * :math:`A_{12} \in \mathbb{R}^{n_1 \times n_2}`
        * :math:`A_{21} \in \mathbb{R}^{n_2 \times n_1}`
        * :math:`A_{22} \in \mathbb{R}^{n_2 \times n_2}` and :math:`A_{22}` is invertible.

    Using that :math:`A_{22}` is invertible leads to the ODE system

    .. math::

        E_{11} \dot{x}_1 &=
        (A_{11} - A_{12}A_{22}^{-1}A_{21})x_1 + (B(1:n_1,:)- A_{12}A_{22}^{-1}B(n_1+1:end,:)) u   \\
                y &= C \begin{bmatrix} x_1 \\ -A_{22}^{-1}A_{21} x_1 \end{bmatrix} -  C \begin{bmatrix} 0 \\
                                        A_{22}^{-1}B(n_1+1:end,:) u \end{bmatrix}

    Therefore we solve the following Riccati equations:

    .. math::

            E_{11} X (A_{11}-A_{12}A_{22}^{-1}A_{21})^T + (A_{11}-A_{12}A_{22}^{-1}A_{21}) X E_{11}^T -
                    E_{11}X C^T C X E_{11}^T + \tilde{B}\tilde{B}^T = 0 \\
                    \tilde{B} := (B(1:n_1,:) - A_{12}A_{22}^{-1}B(n_1+1:end,:))


    .. math::

            E_{11}^T X (A_{11}-A_{12}A_{22}^{-1}A_{21}) + (A_{11}-A_{12}A_{22}^{-1}A_{21})^T X E_{11} =
                    -\tilde{B}^T\tilde{B} \\
                    \tilde{B} := (B(:,1:n_1) - B(:,n_1+1:end) A_{22}^{-1}A_{21})


    Parameters
    ----------
    opt: options instance
        :class:`~pymess.options.Options` instance

    e11: (n1,n1) sparse real matrix
        Matrix :math:`E_{11}` of Index 1 system

    a11: (n1,n1) sparse real matrix
        Matrix :math:`A_{11}` of Index 1 system.

    a12: (n1,n2) sparse real matrix
        Matrix :math:`A_{12}` of Index 1 system.

    a21: (n2,n1) sparse real matrix
        Matrix :math:`A_{21}` of Index 1 system.

    a22: (n2,n2) sparse real matrix
        Matrix :math:`A_{22}` of Index 1 system.

    b: (n1,p) or (n1+n2,p) array_like real matrix

    c: (q,n1+2) or (q,n1) array_like real matrix


    Raises
    ------

    See Also
    --------
    :class:`.EquationGLyapDAE1`
    :func:`pymess.lrnm.lrnm`


    Notes
    -----
    Supported types are ndarray of numpy and csc, coo or csr matrix of scipy.

    References
    ----------


    Examples
    --------
    Further examples :download:`pymess_lrnm_dae1.py <../../example/pymess_lrnm_dae1.py>`


    """

    def __init__(self, opt, e11, a11, a12, a21, a22, b, c):

        super(EquationGRiccatiDAE1, self).__init__("EquationGRiccatiDAE1", opt, e11.shape[0])
        self.e11 = e11
        self.a11 = a11
        self.a12 = a12
        self.a21 = a21
        self.a22 = a22
        self.b = b
        self.c = c

    def __delattr__(self, key):
        raise TypeError(
            "Trying to delete attributes of %r is not allowed" % self)
