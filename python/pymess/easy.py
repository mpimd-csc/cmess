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

r"""Interface to the mess_lyap, mess_direct_create_sylvester_sparsedense, mess_care function."""
import numpy
from pymess._c_interface import _care, _lyap, _sylvester_sparsedense


def care(a, e, b, c):
    r"""Interface to the mess_care function.

    This functions solves the Riccati Equations

    .. math::

        A^T X       + X A       - XBB^TX        + C^T C      &= 0          \\
        A^T X E     + E^TX A    - E^TXBB^TXE    + C^T C      &= 0


    Parameters
    ----------
    a: (n,n) array_like or sparse real matrix
        Matrix :math:`A` of Riccati Equation.

    e: (n,n) array_like or sparse real matrix, optional
        Matrix :math:`E` of Riccati equation.
        If e is None, then :math:`E` of Riccati equation will be the identity matrix.

    b: (n,p)  array_like real matrix
        The matrix :math:`B` of Riccati Equation.
        We assume that :math:`p\ll n`.

    c: (q,n)  array_like real matrix
        The matrix :math:`C` of Riccati Equation.
        We assume that :math:`q\ll n`.

    Returns
    -------
    Z: (n,s) array_like, Low Rank solution factor.
        :math:`X\approx ZZ^T`.

    Raises
    ------

    See Also
    --------

    Notes
    -----
    Supported types are ndarray of numpy and csc, coo or csr matrix of scipy.

    References
    ----------

    Examples
    --------

    Solve standard Riccati Equation :math:`A^T X + X A - XBB^TX + C^T C  = 0`.

    .. doctest::

        >>> from pymess import *
        >>> from scipy.linalg import norm
        >>> import numpy as np
        >>> a = np.array([[-1., 0.1],[ 0., -2.]])
        >>> b = np.array([[1.], [2.]])
        >>> c = np.array([[1., 0.]])
        >>> z = care(a,None,b,c)
        >>> x = z.dot(z.T)
        >>> x1 = a.T.dot(x)
        >>> print(norm(x1 + x1.T - x.dot(b).dot(b.T).dot(x) + c.T.dot(c))<1e-8)
        True

    Solve generalized Riccati Equation :math:`A^T X E + E^T X A - E^T XBB^TX E + C^T C = 0`.

    .. doctest::

        >>> from pymess import *
        >>> from scipy.linalg import norm
        >>> import numpy as np
        >>> a = np.array([[-1., 0.1], [0., -2.]])
        >>> b = np.array([[1.], [2.]])
        >>> c = np.array([[1., 0.]])
        >>> e = np.array([[0.5, 0.1], [0., 0.5]])
        >>> z = care(a,e,b,c)
        >>> x = z.dot(z.T)
        >>> x1 = a.T.dot(x).dot(e)
        >>> print(norm(x1 + x1.T - e.T.dot(x).dot(b).dot(b.T).dot(x).dot(e) + c.T.dot(c))<1e-8)
        True

    Further examples :download:`pymess_care.py <../../example/pymess_care.py>`

    """

    z = _care(a, e, b, c)

    if isinstance(b, numpy.ndarray) and isinstance(c, numpy.ndarray):
        return numpy.array(z)

    if isinstance(b, numpy.matrix) and isinstance(c, numpy.ndarray) or \
            isinstance(b, numpy.ndarray) and isinstance(c, numpy.matrix):
        return numpy.matrix(z)

    return z


def sylvester_sparsedense(a, f, e, b, m):
    r"""Interface to the mess_direct_create_sylvester_sparsedense function.

    This functions solves the Sylvester Equations

    .. math::

        A X       + X B        + M      &= 0          \\
        A X       + E X B        + M      &= 0        \\
        A X F       + E X B        + M      &= 0

    with sparse A, E and small dense B, F.



    Parameters
    ----------
    a: (n,n) sparse matrix
        Matrix :math:`A` of Sylvester Equation.

    f: (m,m) dense matrix, optional
        Matrix :math:`F` of Sylvester equation.
        If f is None, then :math:`F` of Sylvester equation will be the identity matrix.

    e: (n,n) sparse matrix, optional
        Matrix :math:`E` of Sylvester equation.
        If e is None, then :math:`E` of Sylvester equation will be the identity matrix.

    b: (m,m) dense matrix
        The matrix :math:`B` of Sylvester Equation.

    m: (n,m) dense matrix
        The matrix :math:`M` of Sylvester Equation.

    Returns
    -------
    x: (n,m) dense solution matrix :math:`X`.

    Raises
    ------

    See Also
    --------


    Notes
    -----


    References
    ----------
    :cite:`morBenKS11`


    Examples
    --------

    Further examples :download:`pymess_sylvester_sparsedense.py <../../example/pymess_sylvester_sparsedense.py>`

    """

    x = _sylvester_sparsedense(a, f, e, b, m)

    return x


def lyap(a, e, b):
    r"""Interface to the mess_lyap function.

    This functions solves the Lyapunov Equations

    .. math::

        A X       + X A^T       &= -B B^T,          \\
        A X E^T   + E X A^T     &= -B B^T.

    It is import that the eigenvalues of :math:`A` or :math:`(A,E)` lie
    in the left open half plane.

    Parameters
    ----------
    a: (n,n) array_like or sparse real matrix
        Matrix :math:`A` of Lyapunov Equation.

    e: (n,n) array_like or sparse real matrix, optional
        Matrix :math:`E` of Lyapunov equation.
        If e is None, then :math:`E` of Lyapunov equation will be the identity matrix.

    rhs: (n,p)  array_like real matrix
        The right hand side :math:`B` of Lyapunov Equation.
        We assume that :math:`p\ll n`.

    Returns
    -------
    z: (n,s) array_like, Low Rank solution factor.
        :math:`X\approx ZZ^T`.


    Raises
    ------

    See Also
    --------

    Notes
    -----
    Supported types are ndarray of numpy and csc, coo or csr matrix of scipy.

    References
    ----------

    Examples
    --------

    Solve standard Lyapunov Equation :math:`A X + X A^ T = -BB^T`

    .. doctest::

        >>> from pymess import *
        >>> from scipy.linalg import norm
        >>> import numpy as np
        >>> a = np.array([[-1., 0.1], [0., -2.]])
        >>> b = np.array([[1.], [2.]])
        >>> z = lyap(a, None, b)
        >>> x = z.dot(z.T)
        >>> x1 = a.dot(x)
        >>> print(norm(x1 + x1.T + b.dot(b.T))<1e-8)
        True

    Solve generalized Lyapunov Equation :math:`A X E^T + E X A^ T = -BB^T`

    .. doctest::

        >>> from pymess import *
        >>> from scipy.linalg import norm
        >>> import numpy as np
        >>> a = np.array([[-1., 0.1], [0., -2.]])
        >>> b = np.array([[1.], [2.]])
        >>> e = np.array([[0.5, 0.1], [0., 0.5]])
        >>> z = lyap(a, e, b)
        >>> x = z.dot(z.T)
        >>> x1 = a.dot(x).dot(e.T)
        >>> print(norm(x1 + x1.T + b.dot(b.T))<1e-8)
        True

    Further examples :download:`pymess_lyap.py <../../example/pymess_lyap.py>`

    """

    z = _lyap(a, e, b)

    if isinstance(b, numpy.ndarray):
        return numpy.array(z)

    if isinstance(b, numpy.matrix):
        return numpy.matrix(z)
    return z
