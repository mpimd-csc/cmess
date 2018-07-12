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

r"""This module contains methods for residual computation of Lyapunov and Riccati equations."""

from scipy.sparse.linalg import eigsh, LinearOperator
from scipy.sparse import issparse
from scipy.linalg import norm
import numpy as np
from .enum import MESS_OP_NONE, MESS_OP_TRANSPOSE
from .misc import MESSErrorArguments


def res2_lyap(a, e, rhs, z, trans):
    r"""This function computes the residual of a standard/generalized Lyapunov equation.

    This function computes the absolute and relative residual of a Lyapunov equation


    standard Lyapunov equation:

    .. math::

        ||A ZZ^T + ZZ^T A^T   + BB^T||_2

        \frac{||A ZZ^T  + ZZ^T A^T  + BB^T||_2} {||BB^T||_2}


    standard (dual) Lyapunov equation:

    .. math::

        ||A^T ZZ^T + ZZ^T A   + C^TC||_2

        \frac{||A^T ZZ^T  + ZZ^T A  + C^TC||_2} {||C^TC||_2 }


    generalized Lyapunov equation:

    .. math::

        ||A ZZ^T E^T + E ZZ^T A^T   + BB^T||_2

        \frac{||A ZZ^T E^T  + E ZZ^T A^T  + BB^T||_2} {||BB^T||_2 }


    generalized (dual) Lyapunov equation:

    .. math::

        ||A^T ZZ^T E + E^T ZZ^T A   + C^TC||_2

        \frac{||A^T ZZ^T E  + E^T ZZ^T A  + C^TC||_2} {||C^T C||_2 }

    Parameters
    ----------
    a: (n,n) array_like or sparse real matrix
        Matrix :math:`A` of Lyapunov Equation.

    e: (n,n) array_like or sparse real matrix, optional
        Matrix :math:`E` of Lyapunov Equation. Set to None if not wanted.

    rhs: (n,p) array_like real matrix or (p,n) array_like real matrix
        The right hand side of Lyapunov equation either :math:`B` or :math:`C`.

    z: (n,s) array_like, Low Rank solution factor.
        The solution of the Lyapunov equation :math:`X\approx ZZ^T`

    trans: int 0 or 1
        0 refers to the equation where the first matrix of the equation is not transposed.
        1 refers to the equation where the first matrix of the equation is transposed.
        Use :data:`pymess.enum.MESS_OP_NONE` for 0  or
        :data:`pymess.enum.MESS_OP_TRANSPOSE` for 1.


    Returns
    -------

    nrm2: real, absolute residual.

        * :math:`||A ZZ^T + ZZ^T A^T + BB^T||_2`
        * :math:`||A ZZ^T E^T + E ZZ^T A^T + BB^T||_2`
        * :math:`||A^T ZZ^T + ZZ^T A + C^TC||_2`
        * :math:`||A^T ZZ^T E + E^T ZZ^T A + C^TC||_2`

    nrmRHS: real, residual of right hand side.
        * :math:`||BB^T||_2`
        * :math:`||C^TC||_2`

    relnrm: real, relative residual.
        * :math:`\frac{ ||AZZ^T + ZZ^T A^T + BB^T ||_2 }{ ||BB^T||_2 }`
        * :math:`\frac{ ||AZZ^T E^T + E ZZ^T A^T + BB^T ||_2 }{ ||BB^T||_2 }`
        * :math:`\frac{ ||A^TZZ^T + ZZ^T A + C^TC ||_2 }{ ||C^TC||_2 }`
        * :math:`\frac{ ||A^TZZ^T E + E^T ZZ^T A + C^TC ||_2 }{ ||C^TC||_2 }`

    Raises
    ------


    See Also
    --------

    Notes
    -----

    References
    ----------

    Examples
    --------
    Further examples :download:`easy_lyap.py <../../unittests/easy/lyap.py>`

    """

    if (trans not in [MESS_OP_NONE, MESS_OP_TRANSPOSE]):
        raise MESSErrorArguments("trans has to be 0 or 1", trans)

    if a.shape[0] >= 100:
        # setup linear operator for eigenvalue computation
        if trans == MESS_OP_NONE:
            nrmrhs = norm(rhs, 2)**2
            if e is None:

                def mv_lyap(v):
                    """ Lyapunov operator for standard Lyapunov equation."""
                    v = np.matrix(v).T
                    return  a.dot(z.dot(z.T.dot(v)))    +   \
                            z.dot(z.T.dot(a.T.dot(v)))  +   \
                            rhs.dot(rhs.T.dot(v))
            else:

                def mv_lyap(v):
                    """ Lyapunov operator for generalized Lyapunov equation."""
                    v = np.matrix(v).T
                    return  a.dot(z.dot(z.T.dot(e.T.dot(v)))) +     \
                            e.dot(z.dot(z.T.dot(a.T.dot(v)))) +     \
                            rhs.dot(rhs.T.dot(v))
        else:
            nrmrhs = norm(rhs, 2)**2
            if e is None:

                def mv_lyap(v):
                    """ Lyapunov operator for dual Lyapunov equation."""
                    v = np.matrix(v).T
                    return  a.t.dot(z.dot(z.T.dot(v)))  +   \
                            z.dot(z.T.dot(a.dot(v)))    +   \
                            rhs.T.dot(rhs.dot(v))
            else:

                def mv_lyap(v):
                    """ Lyapunov operator for generalized dual Lyapunov equation."""
                    v = np.matrix(v).T
                    return  a.T.dot(z.dot(z.T.dot(e.dot(v))))   +   \
                            e.T.dot(z.dot(z.T.dot(a.dot(v))))   +   \
                            rhs.T.dot(rhs.dot(v))

        lyap = LinearOperator(a.shape, matvec=mv_lyap, rmatvec=mv_lyap, matmat=mv_lyap, dtype=a.dtype)
        nrm2 = abs(eigsh(lyap, 1, None, which='LM', return_eigenvectors=False))
        nrm2 = nrm2[0]
    else:
        # matrix is dense so compute residual dense
        if trans == MESS_OP_NONE:
            nrmrhs = norm(rhs, 2)**2
            if e is None:
                nrm2 = norm(a.dot(z.dot(z.T)) + z.dot(z.T.dot(a.T)) + rhs.dot(rhs.T), 2)
            else:
                nrm2 = norm(a.dot(z.dot(z.T.dot(e.T))) + e.dot(z.dot(z.T.dot(a.T))) + rhs.dot(rhs.T), 2)
        else:
            nrmrhs = norm(rhs, 2)**2
            if e is None:
                nrm2 = norm(a.T.dot(z.dot(z.T)) + z.dot(z.T.dot(a)) + rhs.T.dot(rhs), 2)
            else:
                nrm2 = norm(a.T.dot(z.dot(z.T.dot(e))) + e.T.dot(z.dot(z.T.dot(a))) + rhs.T.dot(rhs), 2)

    return nrm2, nrmrhs, nrm2 / nrmrhs


def res2_glyap(a, e, y, x, trans):
    r"""This function computes the residual of a standard/generalized Lyapunov equation.

    This function computes the absolute and relative residual of a Lyapunov equation


    standard Lyapunov equation:

    .. math::

        ||A X + X A^T + Y||_2

        \frac{||A X  + X A^T  + Y||_2} {||Y||_2}


    standard (dual) Lyapunov equation:

    .. math::

        ||A^T X + X A + Y||_2

        \frac{||A^T X  + X A  + Y||_2} {||Y||_2 }


    generalized Lyapunov equation:

    .. math::

        ||A X E^T + E X A^T   + Y||_2

        \frac{||A X E^T  + E X A^T  + Y||_2} {||Y||_2 }


    generalized (dual) Lyapunov equation:

    .. math::

        ||A^T X E + E^T X A   + Y||_2

        \frac{||A^T X E  + E^T X A  + Y||_2} {||Y||_2 }

    Parameters
    ----------
    a: (n,n) array_like real matrix
        Matrix :math:`A` of Lyapunov equation.

    e: (n,n) array_like real matrix, optional
        Matrix :math:`E` of Lyapunov equation. Set to None if not wanted.

    y: (n,n) array_like real matrix
        The right hand side of Lyapunov equation :math:`Y`.

    x: (n,n) array_like, numerical solution.
        The numerical solution of the Lyapunov equation :math:`X`

    trans: int 0 or 1
        0 refers to the equation where the first matrix of the equation is not transposed.
        1 refers to the equation where the first matrix of the equation is transposed.
        Use :data:`pymess.enum.MESS_OP_NONE` for 0  or
        :data:`pymess.enum.MESS_OP_TRANSPOSE` for 1.


    Returns
    -------

    nrm2: real, absolute residual.

        * :math:`||A X + X A^T + Y||_2`
        * :math:`||A X E^T + E X A^T + Y||_2`
        * :math:`||A^T X + X A + Y||_2`
        * :math:`||A^T X E + E^T X A + Y||_2`

    nrmRHS: real, residual of right hand side.
        * :math:`||Y||_2`

    relnrm: real, relative residual.
        * :math:`\frac{ ||A X + X A^T + Y ||_2 }{ ||Y||_2 }`
        * :math:`\frac{ ||A X E^T + E X A^T + Y ||_2 }{ ||Y||_2 }`
        * :math:`\frac{ ||A^T X + X A + Y ||_2 }{ ||Y||_2 }`
        * :math:`\frac{ ||A^T X E + E^T X A + Y ||_2 }{ ||Y||_2 }`

    Raises
    ------


    See Also
    --------

    Notes
    -----

    References
    ----------

    Examples
    --------
    Further examples :download:`glyap3.py <../../unittests/glyap3/glyap3.py>`

    """

    if (not trans in [MESS_OP_NONE, MESS_OP_TRANSPOSE]):
        raise MESSErrorArguments("trans has to be 0 or 1", trans)

    if issparse(a) or issparse(y) or e is not None and issparse(e):
        raise MESSErrorArguments("Matrices must be dense.")

    nrmrhs = norm(y, 2)

    # matrix is dense so compute residual dense
    if trans == MESS_OP_NONE:
        if e is None:
            nrm2 = norm(a.dot(x) + x.dot(a.T) + y, 2)
        else:
            nrm2 = norm(a.dot(x.dot(e.T)) + e.dot(x.dot(a.T)) + y, 2)
    else:
        if e is None:
            nrm2 = norm(a.T.dot(x) + x.dot(a) + y, 2)
        else:
            nrm2 = norm(a.T.dot(x.dot(e)) + e.T.dot(x.dot(a)) + y, 2)

    return nrm2, nrmrhs, nrm2 / nrmrhs


def res2_gstein(a, e, y, x, trans):
    r"""This function computes the residual of a standard/generalized Stein equation.

    This function computes the absolute and relative residual of a Stein equation


    standard Stein equation:

    .. math::

        ||A X A^T - X + Y||_2

        \frac{||A X A^T  - X  + Y||_2} {||Y||_2}


    standard (dual) Stein equation:

    .. math::

        ||A^T X A - X + Y||_2

        \frac{||A^T X A - X  + Y||_2} {||Y||_2 }


    generalized Stein equation:

    .. math::

        ||A X A^T - E X E^T + Y||_2

        \frac{||A X A^T  - E X E^T  + Y||_2} {||Y||_2 }


    generalized (dual) Stein equation:

    .. math::

        ||A^T X A - E^T X E   + Y||_2

        \frac{||A^T X A  - E^T X E  + Y||_2} {||Y||_2 }

    Parameters
    ----------
    a: (n,n) array_like real matrix
        Matrix :math:`A` of Stein equation.

    e: (n,n) array_like real matrix, optional
        Matrix :math:`E` of Stein equation. Set to None if not wanted.

    y: (n,n) array_like real matrix
        The right hand side of Stein equation :math:`Y`.

    x: (n,n) array_like, numerical solution.
        The numerical solution of the Stein equation :math:`X`

    trans: int 0 or 1
        0 refers to the equation where the first matrix of the equation is not transposed.
        1 refers to the equation where the first matrix of the equation is transposed.
        Use :data:`pymess.enum.MESS_OP_NONE` for 0  or
        :data:`pymess.enum.MESS_OP_TRANSPOSE` for 1.


    Returns
    -------

    nrm2: real, absolute residual.

        * :math:`||A X A^T - X + Y||_2`
        * :math:`||A X A^T - E X E^T + Y||_2`
        * :math:`||A^T X A - X  + Y||_2`
        * :math:`||A^T X A - E^T X E + Y||_2`

    nrmRHS: real, residual of right hand side.
        * :math:`||Y||_2`

    relnrm: real, relative residual.
        * :math:`\frac{ ||A X A^T - X  + Y ||_2 }{ ||Y||_2 }`
        * :math:`\frac{ ||A X A^T - E X E^T + Y ||_2 }{ ||Y||_2 }`
        * :math:`\frac{ ||A^T X A - X  + Y ||_2 }{ ||Y||_2 }`
        * :math:`\frac{ ||A^T X A - E^T X E + Y ||_2 }{ ||Y||_2 }`

    Raises
    ------


    See Also
    --------

    Notes
    -----

    References
    ----------

    Examples
    --------
    Further examples :download:`glyap3.py <../../unittests/glyap3/glyap3.py>`

    """

    if (trans not in [MESS_OP_NONE, MESS_OP_TRANSPOSE]):
        raise MESSErrorArguments("trans has to be 0 or 1", trans)

    if issparse(a) or issparse(y) or e is not None and issparse(e):
        raise MESSErrorArguments("Matrices must be dense.")

    nrmrhs = norm(y, 2)

    # matrix is dense so compute residual dense
    if trans == MESS_OP_NONE:
        if e is None:
            nrm2 = norm(a.dot(x).dot(a.T) - x + y, 2)
        else:
            nrm2 = norm(a.dot(x).dot(a.T) - e.dot(x.dot(e.T)) + y, 2)
    else:
        if e is None:
            nrm2 = norm(a.T.dot(x).dot(a) - x + y, 2)
        else:
            nrm2 = norm(a.T.dot(x).dot(a) - e.T.dot(x.dot(e)) + y, 2)

    return nrm2, nrmrhs, nrm2 / nrmrhs


def res2_sylv_sd(a, f, e, h, rhs, x):
    r"""This function computes the residual of a standard/generalized Sylvester equation.

    This function computes the absolute and relative residual of a Sylvester equation


    standard Sylvester equation:

    .. math::

        ||A X + X H   + M||_2

        \frac{||A X  + X H  + M||_2} {||M||_2}


    semi-generalized Sylvester equation:

    .. math::

        ||A X + E X H   + M||_2

        \frac{||A X  + E X H  + M||_2} {||M||_2 }


    generalized Sylvester equation:

    .. math::

        ||A X F + E X H   + M||_2

        \frac{||A X F  + E X H  + M||_2} {||M||_2 }

    Parameters
    ----------
    a: (n,n) sparse matrix
        Matrix :math:`A` of Sylvester Equation.

    f: (n,n) dense matrix, optional
        Matrix :math:`F` of Sylvester Equation. Set to None if not wanted.

    e: (n,n) sparse matrix, optional
        Matrix :math:`E` of Sylvester Equation. Set to None if not wanted.

    h: (n,n) dense matrix
        Matrix :math:`H` of Sylvester Equation.

    rhs: (n,p) dense matrix
        The right hand side of Sylvester equation :math:`M`

    x: (n,s) dense solution matrix.
        The solution :math:`X` of the Sylvester equation


    Returns
    -------

    nrm2: real, absolute residual.

        * :math:`||A X + X H + M||_2`
        * :math:`||A X + E X H + M||_2`
        * :math:`||A X F + E X H + M||_2`
        * :math:`||A^T X + X H^T + M||_2`
        * :math:`||A^T X + E^T X H^T + M||_2`
        * :math:`||A^T X F^T + E^T X H^T + M||_2`

    nrmRHS: real, residual of right hand side.
        * :math:`||M||_2`

    relnrm: real, relative residual.
        * :math:`\frac{ ||A X + X H + M ||_2 }{ ||M||_2 }`
        * :math:`\frac{ ||A X + E X H + M ||_2 }{ ||M||_2 }`
        * :math:`\frac{ ||A X F + E X H + M ||_2 }{ ||M||_2 }`
        * :math:`\frac{ ||A^T X + X H^T + M ||_2 }{ ||M||_2 }`
        * :math:`\frac{ ||A^T X + E^T X H^T + M ||_2 }{ ||M||_2 }`
        * :math:`\frac{ ||A^T X F^T + E^T X H^T + M ||_2 }{ ||M||_2 }`

    Raises
    ------


    See Also
    --------

    Notes
    -----

    References
    ----------

    Examples
    --------
    Further examples :download:`pymess_sylvester_sparsedense.py <../../example/pymess_sylvester_sparsedense.py>`

    """

    nrmrhs = norm(rhs, 2)

    if e is None:
        nrm2 = norm(a.dot(x) + x.dot(h) + rhs, 2)
    elif f is None:
        nrm2 = norm(a.dot(x) + e.dot(x.dot(h)) + rhs, 2)
    else:
        nrm2 = norm(a.dot(x.dot(f)) + e.dot(x.dot(h)) + rhs, 2)


    return nrm2, nrmrhs, nrm2 / nrmrhs


def res2_ric(a, e, b, c, z, trans):
    r"""This function computes the residual of a standard/generalized Riccati equation.

    This function computes the absolute and relative residual of a Riccati equation

    standard Riccati equation:

    .. math::

        ||A ZZ^T  +  ZZ^T A^T - ZZ^T C^TC ZZ^T + BB^T  ||_2

        \frac{||A ZZ^T  +  ZZ^T A^T - ZZ^T C^TC ZZ^T + BB^T  ||_2} {|| BB^T||_2}


    standard (dual) Riccati equation:

    .. math::

        ||A^T ZZ^T  +  ZZ^T A - ZZ^T BB^T ZZ^T + C^T C  ||_2

        \frac{||A^T ZZ^T  +  ZZ^T A^ - ZZ^T B B^T ZZ^T + C^T C  ||_2} {|| C^T C ||_2}


    generalized Riccati equation:

    .. math::

        ||A ZZ^T E^T  +   E Z Z^T A^T -  E ZZ^T C^TC ZZ^T E^T + BB^T  ||_2

        \frac{||A ZZ^T E^T  +  E  ZZ^T A^T - E ZZ^T C^TC ZZ^T E^T + BB^T  ||_2} {|| BB^T||_2}


    generalized dual Riccati equation:

    .. math::

        ||A^T ZZ^T E  +  E^T  ZZ^T A - E^T ZZ^T BB^T ZZ^T E + C^T C  ||_2

        \frac{||A^T ZZ^T  E + E^T ZZ^T A^ - E^T ZZ^T B B^T ZZ^T + C^T C  ||_2} {|| C^T C ||_2}


    Parameters
    ----------
    A: (n,n) array_like or sparse real matrix
        Matrix :math:`A` of Riccati Equation.

    E: (n,n) array_like or sparse real matrix, optional
        Matrix :math:`E` of Riccati Equation. Set to None if not wanted.

    B: (n,q)  array_like real matrix
        Matrix :math:`B` of Riccati Equation.

    C: (p,n)  array_like real matrix
        Matrix :math:`C` of Riccati Equation.

    Z: (n,s) array_like, Low Rank solution factor.

    trans: int 0 or 1
        0 refers to the equation where the first matrix of the equation is not transposed.
        1 refers to the equation where the first matrix of the equation is transposed.
        You can use :data:`~pymess.enum.MESS_OP_NONE` for 0 or
        :data:`~pymess.enum.MESS_OP_TRANSPOSE` for 1.


    Returns
    -------

    nrm2: real, absolute residual.
        * :math:`||AZZ^T + ZZ^T A^T - ZZ^TC^TCZZ^T + BB^T||_2`
        * :math:`||AZZ^T E^T + E ZZ^T A^T - E ZZ^T C^TC ZZ^T E^T + BB^T||_2`
        * :math:`||A^TZZ^T + ZZ^T A - ZZ^T B B^T ZZ^T + C^T C||_2`
        * :math:`||A^TZZ^TE + E^T ZZ^T A - E^T ZZ^T B B^T ZZ^T E + C^T C||_2`


    nrmRHS: real, residual of right hand side.
        * :math:`||BB^T||_2`
        * :math:`||C^T C||_2`

    relnrm: real, relative residual.
        * :math:`\frac{||AZZ^T + ZZ^T A^T - ZZ^T C^TC ZZ^T + BB^T||_2} {||BB^T||_2}`
        * :math:`\frac{||AZZ^T E^T + E ZZ^T A^T - E ZZ^T C^TC ZZ^T E^T + BB^T||_2} {||BB^T||_2}`
        * :math:`\frac{||A^TZZ^T + ZZ^T A - ZZ^T BB^T ZZ^T + C^T C  ||_2} {||C^T C||_2}`
        * :math:`\frac{||A^TZZ^TE + E^TZZ^T A - E^T ZZ^T B B^T ZZ^T E + C^T C||_2}{||C^T C||_2}`

    Raises
    ------

    See Also
    --------


    Notes
    -----


    References
    ----------


    Examples
    --------
    Further examples :download:`easy_care.py <../../unittests/easy/care.py>`

    """

    if (not trans in [MESS_OP_NONE, MESS_OP_TRANSPOSE]):
        raise MESSErrorArguments("trans has to be 0 or 1", trans)

    if a.shape[0] >= 100:
        # setup linear operator for eigenvalue computation
        if trans == MESS_OP_NONE:
            nrmrhs = norm(b, 2)**2
            if e is None:

                def mv_ric(v):
                    """ Riccati operator for standard Riccati equation."""
                    v = np.matrix(v).T
                    return  a.dot(z.dot(z.T.dot(v)))                            +   \
                            z.dot(z.T.dot(a.T.dot(v)))                          -   \
                            z.dot(z.T.dot(c.T.dot(c.dot(z.dot(z.T.dot(v))))))   +   \
                            b.dot(b.T.dot(v))
            else:

                def mv_ric(v):
                    """ Riccati operator for generalized Riccati equation."""
                    v = np.matrix(v).T
                    return  a.dot(z.dot(z.T.dot(e.T.dot(v))))                                 + \
                            e.dot(z.dot(z.T.dot(a.T.dot(v))))                                 - \
                            e.dot(z.dot(z.T.dot(c.T.dot(c.dot(z.dot(z.T.dot(e.T.dot(v)))))))) + \
                            b.dot(b.T.dot(v))
        else:
            nrmrhs = norm(c, 2)**2

            if e is None:

                def mv_ric(v):
                    """ Riccati operator for dual standard Riccati equation."""
                    v = np.matrix(v).T
                    return  a.T.dot(z.dot(z.T.dot(v)))                              +   \
                            z.dot(z.T.dot(a.dot(v)))                                -   \
                            z.dot(z.T.dot(b.dot(b.T.dot(z.dot(z.T.dot(v))))))       +   \
                            c.T.dot(c.dot(v))
            else:

                def mv_ric(v):
                    """ Riccati operator for dual generalized Riccati equation."""
                    v = np.matrix(v).T
                    return  a.T.dot(z.dot(z.T.dot(e.dot(v))))                                 + \
                            e.T.dot(z.dot(z.T.dot(a.dot(v))))                                 - \
                            e.T.dot(z.dot(z.T.dot(b.dot(b.T.dot(z.dot(z.T.dot(e.dot(v)))))))) + \
                            c.T.dot(c.dot(v))

        ric = LinearOperator(a.shape, matvec=mv_ric, rmatvec=mv_ric, matmat=mv_ric, dtype=a.dtype)
        nrm2 = abs(eigsh(ric, 1, None, which='LM', return_eigenvectors=False))
        nrm2 = nrm2[0]
    else:
        if trans == MESS_OP_NONE:
            nrmrhs = norm(b, 2)**2
            if e is None:
                nrm2 = norm(a.dot(z.dot(z.T)) + z.dot(z.T.dot(a.T)) -                   \
                        z.dot(z.T.dot(c.T.dot(c.dot(z.dot(z.T))))) + b.dot(b.T), 2)
            else:
                nrm2 = norm(a.dot(z.dot(z.T.dot(e.T))) + e.dot(z.dot(z.T.dot(a.T))) -       \
                        e.dot(z.dot(z.T.dot(c.T.dot(c.dot(z.dot(z.T.dot(e.T))))))) +    \
                        b.dot(b.T), 2)

        else:
            nrmrhs = norm(c, 2)**2
            if e is None:
                nrm2 = norm(a.T.dot(z.dot(z.T)) + z.dot(z.T.dot(a)) -                       \
                        z.dot(z.T.dot(b.dot(b.T.dot(z.dot(z.T))))) + c.T.dot(c), 2)
            else:
                nrm2 = norm(a.T.dot(z.dot(z.T.dot(e))) + e.T.dot(z.dot(z.T.dot(a))) -       \
                        e.T.dot(z.dot(z.T.dot(b.dot(b.T.dot(z.dot(z.T.dot(e))))))) +    \
                        c.T.dot(c), 2)

    return nrm2, nrmrhs, nrm2 / nrmrhs

def res2_gmpare(a, e, q, g, x, plus=False, trans=MESS_OP_TRANSPOSE):
    r"""This function computes the residual of a dense standard/generalized positive/negative Riccati equation.

    This function computes the absolute and relative residual of a Riccati equation

    standard Riccati equation:

    .. math::

        ||A^T X  +  X A +/- X G X + Q  ||_2

        \frac{||A^T X  +  X A +/- X G X + Q  ||_2} {|| Q ||_2}


    standard dual Riccati equation:

    .. math::

        ||A X  +  X A^T +/- X G X + Q  ||_2

        \frac{||A X  +  X A^T +/- X G X + Q  ||_2} {|| Q ||_2}


    generalized Riccati equation:

    .. math::

        ||A^T X E  +  E^T X A +/- E^T X G X E + Q  ||_2

        \frac{||A^T X E  +  E^T X A +/- E^T X G X E + Q  ||_2} {|| Q ||_2}


    generalized dual Riccati equation:

    .. math::

        ||A X E^T  +  E X A^T +/- E X G X E^T + Q  ||_2

        \frac{||A X E^T  +  E X A^T +/- E X G X E^T + Q  ||_2} {|| Q ||_2}



    Parameters
    ----------
    A: (n,n) symmetric dense real matrix
        Matrix :math:`A` of Riccati Equation.

    E: (n,n) symmetric dense real matrix, optional
        Matrix :math:`E` of Riccati Equation. Set to None if not wanted.

    Q: (n,n)  symmetric dense real matrix
        Matrix :math:`Q` of Riccati Equation.

    G: (n,n)  symmetric dense real matrix
        Matrix :math:`G` of Riccati Equation.

    X: (n,n) symmetric dense solution matrix.

    plus: bool
        True refers to the positive equation.
        False refers to the negative equation.

    trans: int 0 or 1
        0 refers to the equation where the first matrix of the equation is not transposed.
        1 refers to the equation where the first matrix of the equation is transposed.
        You can use :data:`~pymess.enum.MESS_OP_NONE` for 0 or
        :data:`~pymess.enum.MESS_OP_TRANSPOSE` for 1.


    Returns
    -------

    nrm2: real, absolute residual.
        * :math:`||A^TX + XA +/- XGX + Q||_2`
        * :math:`||A^TXE + E^TXA +/- E^TXGXE + Q||_2`
        * :math:`||AX + XA^T +/- XGX + Q||_2`
        * :math:`||AXE^T + EXA^T +/- EXGXE^T + Q||_2`


    nrmRHS: real, residual of right hand side.
        * :math:`||Q||_2`

    relnrm: real, relative residual.
        * :math:`\frac{||A^TX + XA +/- XGX + Q||_2} {||Q||_2}`
        * :math:`\frac{||A^TXE + E^TXA +/- E^TXGXE + Q||_2} {||Q||_2}`
        * :math:`\frac{||AX + XA^T +/- XGX + Q||_2} {||Q||_2}`
        * :math:`\frac{||AXE^T + EXA^T +/- EXGXE^T + Q||_2}{||Q||_2}`

    Raises
    ------

    See Also
    --------


    Notes
    -----


    References
    ----------


    Examples
    --------
    Further examples :download:`pymess_gmpare.py <../../example/pymess_gmpare.py>`

    """

    if (trans not in [MESS_OP_NONE, MESS_OP_TRANSPOSE]):
        raise MESSErrorArguments("trans has to be 0 or 1", trans)
    if (plus not in [False, True]):
        raise MESSErrorArguments("plus has to be False or True", plus)

    nrmrhs = norm(q, 2)

    plus_scalar = 1 if plus else -1

    if trans == MESS_OP_TRANSPOSE:
        if e is None:
            nrm2 = norm(a.T.dot(x) + x.dot(a) + plus_scalar*x.dot(g.dot(x)) + q, 2)
        else:
            nrm2 = norm(a.T.dot(x.dot(e)) + e.T.dot(x.dot(a)) + plus_scalar*e.T.dot(x.dot(g.dot(x.dot(e)))) + q, 2)
    else:
        if e is None:
            nrm2 = norm(a.dot(x) + x.dot(a.T) + plus_scalar*x.dot(g.dot(x)) + q, 2)
        else:
            nrm2 = norm(a.dot(x.dot(e.T)) + e.dot(x.dot(a.T)) + plus_scalar*e.dot(x.dot(g.dot(x.dot(e.T)))) + q, 2)

    return nrm2, nrmrhs, nrm2 / nrmrhs
