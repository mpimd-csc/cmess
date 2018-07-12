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

r"""This module contains enum types from |cmess| defined as constants."""

#: equation types
MESS_EQN_NONE = 1
r"""
Represents no equation. Mostly used for initialization.
"""

MESS_EQN_LYAP = 1
r"""
Represents a Lyapunov Equation.

Represents the Lyapunov Equations:

.. math::
    :nowrap:

    \begin{eqnarray*}
    A^TX + XA   &=& -B^TB   \\
    AX   + XA^T &=& -CC^T.
    \end{eqnarray*}
"""

MESS_EQN_GLYAP = 2
r"""
Represents a generalized Lyapunov Equation.

Represents the generalized Lyapunov Equations:

.. math::
    :nowrap:

    \begin{eqnarray*}
    A^TXE + E^TXA   &=& -B^TB   \\
    AXE^T + EXA^T   &=& -CC^T.
    \end{eqnarray*}
"""

MESS_EQN_RICCATI = 3
r"""
Represents a Riccati Equation.

Represents the Riccati Equations:

.. math::
    :nowrap:

    \begin{eqnarray*}
    A^TX + XA       - X BB^TX   + C^T C  &=& 0  \\
    AX   + XA^T     - X C^T CX  + B B^T  &=& 0
    \end{eqnarray*}
"""

MESS_EQN_GRICCATI = 4
r"""
Represents a generalized Riccati Equation.

Represents the generalized Riccati Equations:

.. math::
    :nowrap:

    \begin{eqnarray*}
    A^TXM + M^TXA   - M^TX BB^TXM   + C^T C  &=& 0  \\
    AXM^T + MXA^T   - MX C^T CXM^T  + B B^T  &=& 0
    \end{eqnarray*}
"""

# memory usage type
MESS_MEMORY_LOW = 0
r"""
 Option for low memory usage.
"""

MESS_MEMORY_MID = 1
r"""
 Option for medium memory usage.
"""

MESS_MEMORY_HIGH = 2
r"""
 Option for high memory usage.
"""

# operation types
MESS_OP_NONE = 0
r"""
Represents the non-transposed operation.
"""

MESS_OP_TRANSPOSE = 1
r"""
Represents the transposed operation.
"""

MESS_OP_HERMITIAN = 2
r"""
Represents the hermitian operation.
"""

#Shift Parameter Options

MESS_LRCFADI_PARA_MINMAX = 0
r"""
Heuristic shift parameters.

Compute complex parameters using the heuristic for the min-max problem.
"""

MESS_LRCFADI_PARA_MINMAX_REAL = 1
r"""
Heuristic real shift parameters.

Compute real parameters using the heuristic for the min-max problem.
"""

MESS_LRCFADI_PARA_WACHSPRESS = 2
r"""
Wachspress parameters.

Compute Wachspress shift parameters.
"""

MESS_LRCFADI_PARA_ADAPTIVE_V = 3
r"""
Adaptive shift parameters, using :math:`V`.

Adaptive shift parameters, using :math:`V` as projection space.
"""

MESS_LRCFADI_PARA_ADAPTIVE_Z = 4
r"""
Adaptive shift parameters, using :math:`Z`.

Adaptive shift parameters, using :math:`Z` as projection space.
"""

# LINEAR SOLVERS

MESS_DIRECT_DEFAULT_LU = 0
r"""
Select the automatic solver detection.
"""

MESS_DIRECT_SPARSE_LU = 1
r"""
Select the internal sparse solver.
"""

MESS_DIRECT_LAPACK_LU = 2
r"""
Select the lapack dense solver.
"""

MESS_DIRECT_UMFPACK_LU = 3
r"""
Select the umfpack sparse solver.
"""

MESS_DIRECT_SUPERLU_LU = 4
r"""
Select the superlu sparse solver.
"""

MESS_DIRECT_CSPARSE_LU = 5
r"""
Select the csparse sparse solver.
"""

MESS_DIRECT_BANDED_LU = 6
r"""
Select the banded sparse solver.
"""

MESS_DIRECT_MKLPARDISO_LU = 7
r"""
Select the  mklpardiso solver.
"""

# CHOLESKY BASED SOLVERS
MESS_DIRECT_DEFAULT_CHOLESKY = 0
r"""
Select the default cholesky solver.
"""

MESS_DIRECT_LAPACK_CHOLESKY = 1
r"""
Select the default lapack cholesky solver.
"""

MESS_DIRECT_CSPARSE_CHOLESKY = 2
r"""
Select the default csparse cholesky solver.
"""

MESS_DIRECT_CHOLMOD_CHOLESKY = 3
r"""
Select the default cholmod cholesky solver.
"""


# MULTI DIRECT SOLVERS
MESS_MULTIDIRECT_SPARSE_LU = 0
r"""
Select the internal Single-Pattern Multi-Value multisolver.
"""

MESS_MULTIDIRECT_UMFPACK_LU = 1
r"""
Select the umfpack based multisolver.
"""


# RESIDUAL METHODS
MESS_RESIDUAL_INDEFINITE = 0
r"""
Residual via Residual Factor Updates (fastest variant).
"""

MESS_RESIDUAL_SPECTRAL = 1
r"""
Compute Largest Eigenvalue iteratively.
"""


# NORM TYPES
MESS_2_NORM = 0
r"""
2-Norm.
"""

MESS_FROBENIUS_NORM = 1
r"""
Frobenius-Norm.
"""

MESS_1_NORM = 3
r"""
1-Norm, maximum absolute column sum norm.
"""

MESS_INF_NORM = 4
r"""
Infinity-Norm, maximum absolute row sum norm.
"""
