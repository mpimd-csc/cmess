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


Py-M.E.S.S. Enums
=========================

.. currentmodule:: pymess.enum


Matrix Equation Enums
-----------------------

.. autosummary::

    MESS_EQN_NONE
    MESS_EQN_LYAP
    MESS_EQN_GLYAP
    MESS_EQN_RICCATI
    MESS_EQN_GRICCATI

.. autodata:: MESS_EQN_NONE
    :annotation:
.. autodata:: MESS_EQN_LYAP
    :annotation:
.. autodata:: MESS_EQN_GLYAP
    :annotation:
.. autodata:: MESS_EQN_RICCATI
    :annotation:
.. autodata:: MESS_EQN_GRICCATI
    :annotation:


Memory Usage Enums
--------------------

.. autosummary::

    MESS_MEMORY_LOW
    MESS_MEMORY_MID
    MESS_MEMORY_HIGH

.. autodata:: MESS_MEMORY_LOW
    :annotation:
.. autodata:: MESS_MEMORY_MID
    :annotation:
.. autodata:: MESS_MEMORY_HIGH
    :annotation:


Shift Parameter Enums
------------------------------

.. autosummary::

    MESS_LRCFADI_PARA_MINMAX
    MESS_LRCFADI_PARA_MINMAX_REAL
    MESS_LRCFADI_PARA_WACHSPRESS
    MESS_LRCFADI_PARA_ADAPTIVE_V
    MESS_LRCFADI_PARA_ADAPTIVE_Z

.. autodata:: MESS_LRCFADI_PARA_MINMAX
    :annotation:
.. autodata:: MESS_LRCFADI_PARA_MINMAX_REAL
    :annotation:
.. autodata:: MESS_LRCFADI_PARA_WACHSPRESS
    :annotation:
.. autodata:: MESS_LRCFADI_PARA_ADAPTIVE_V
    :annotation:
.. autodata:: MESS_LRCFADI_PARA_ADAPTIVE_Z
    :annotation:


Operation Type Enums
-----------------------

.. autosummary::

    MESS_OP_NONE
    MESS_OP_TRANSPOSE
    MESS_OP_HERMITIAN

.. autodata:: MESS_OP_NONE
    :annotation:
.. autodata:: MESS_OP_TRANSPOSE
    :annotation:
.. autodata:: MESS_OP_HERMITIAN
    :annotation:

Linear Solver Enums
---------------------

.. autosummary::

    MESS_DIRECT_DEFAULT_LU
    MESS_DIRECT_SPARSE_LU
    MESS_DIRECT_LAPACK_LU
    MESS_DIRECT_UMFPACK_LU
    MESS_DIRECT_SUPERLU_LU
    MESS_DIRECT_CSPARSE_LU
    MESS_DIRECT_BANDED_LU
    MESS_DIRECT_MKLPARDISO_LU

.. autodata:: MESS_DIRECT_DEFAULT_LU
    :annotation:
.. autodata:: MESS_DIRECT_SPARSE_LU
    :annotation:
.. autodata:: MESS_DIRECT_LAPACK_LU
    :annotation:
.. autodata:: MESS_DIRECT_UMFPACK_LU
    :annotation:
.. autodata:: MESS_DIRECT_SUPERLU_LU
    :annotation:
.. autodata:: MESS_DIRECT_CSPARSE_LU
    :annotation:
.. autodata:: MESS_DIRECT_BANDED_LU
    :annotation:
.. autodata:: MESS_DIRECT_MKLPARDISO_LU
    :annotation:


Cholesky Decomposition Enums
----------------------------

.. autosummary::

    MESS_DIRECT_DEFAULT_CHOLESKY
    MESS_DIRECT_LAPACK_CHOLESKY
    MESS_DIRECT_CSPARSE_CHOLESKY
    MESS_DIRECT_CHOLMOD_CHOLESKY

.. autodata:: MESS_DIRECT_DEFAULT_CHOLESKY
    :annotation:
.. autodata:: MESS_DIRECT_LAPACK_CHOLESKY
    :annotation:
.. autodata:: MESS_DIRECT_CSPARSE_CHOLESKY
    :annotation:
.. autodata:: MESS_DIRECT_CHOLMOD_CHOLESKY
    :annotation:


Multidirect Solver Enums
----------------------------

.. autosummary::
    MESS_MULTIDIRECT_SPARSE_LU
    MESS_MULTIDIRECT_UMFPACK_LU


.. autodata:: MESS_MULTIDIRECT_SPARSE_LU
    :annotation:
.. autodata:: MESS_MULTIDIRECT_UMFPACK_LU
    :annotation:

Residual Methods Enums
----------------------------

.. autosummary::
    MESS_RESIDUAL_INDEFINITE
    MESS_RESIDUAL_SPECTRAL

.. autodata:: MESS_RESIDUAL_INDEFINITE
    :annotation:
.. autodata:: MESS_RESIDUAL_SPECTRAL
    :annotation:

Norm Methods Enums
----------------------------

.. autosummary::
    MESS_2_NORM
    MESS_FROBENIUS_NORM
    MESS_1_NORM
    MESS_INF_NORM

.. autodata:: MESS_2_NORM
    :annotation:
.. autodata:: MESS_FROBENIUS_NORM
    :annotation:
.. autodata:: MESS_1_NORM
    :annotation:
.. autodata:: MESS_INF_NORM
    :annotation:
