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

r"""Interface to the mess_direct_select, mess_direct_chol_select, mess_multidirect_select functions."""

from pymess._c_interface import _direct_select, _direct_chol_select, _multidirect_select


def direct_select(lu_package):
    r"""Interface to the mess_direct_select function.


    Parameters
    ----------
    lu package:

        - MESS_DIRECT_DEFAULT_LU
        - MESS_DIRECT_SPARSE_LU
        - MESS_DIRECT_LAPACKLU
        - MESS_DIRECT_UMFPACK_LU
        - MESS_DIRECT_SUPERLU_LU
        - MESS_DIRECT_CSPARSELU
        - MESS_DIRECT_BANDED_LU
        - MESS_DIRECT_MKLPARDISO_LU

    Returns
    -------

    Raises
    ------
    RuntimeError if solver is not available

    See Also
    --------

    Notes
    ------

    References
    ----------

    Examples
    --------

    .. doctest::

        >>> from pymess import *
        >>> direct_select(MESS_DIRECT_DEFAULT_LU)       # doctest: +SKIP
        >>> direct_select(MESS_DIRECT_SPARSE_LU)        # doctest: +SKIP
        >>> direct_select(MESS_DIRECT_LAPACKLU)         # doctest: +SKIP
        >>> direct_select(MESS_DIRECT_UMFPACK_LU)       # doctest: +SKIP
        >>> direct_select(MESS_DIRECT_SUPERLU_LU)       # doctest: +SKIP
        >>> direct_select(MESS_DIRECT_CSPARSELU)        # doctest: +SKIP
        >>> direct_select(MESS_DIRECT_BANDED_LU)        # doctest: +SKIP
        >>> direct_select(MESS_DIRECT_MKLPARDISO_LU)    # doctest: +SKIP
        >>> direct_select(MESS_DIRECT_DEFAULT_LU)       # doctest: +SKIP

    """

    return _direct_select(lu_package)


def direct_chol_select(chol_package):
    r"""Interface to the mess_direct_chol_select function.


    Parameters
    ----------
    chol_package:

        - MESS_DIRECT_DEFAULTCHOL
        - MESS_DIRECT_LAPACKCHOL
        - MESS_DIRECT_CSPARSECHOL
        - MESS_DIRECT_CHOLMODCHOL

    Returns
    -------

    Raises
    ------
    RuntimeError if solver is not available

    See Also
    --------

    Notes
    -----

    References
    ----------

    Examples
    ---------

    .. doctest::

        >>> from pymess import *
        >>> direct_chol_select(MESS_DIRECT_DEFAULTCHOL)     # doctest: +SKIP
        >>> direct_chol_select(MESS_DIRECT_LAPACKCHOL)      # doctest: +SKIP
        >>> direct_chol_select(MESS_DIRECT_CSPARSECHOL)     # doctest: +SKIP
        >>> direct_chol_select(MESS_DIRECT_CHOLMODCHOL)     # doctest: +SKIP
        >>> direct_chol_select(MESS_DIRECT_DEFAULTCHOL)     # doctest: +SKIP
    """

    return _direct_chol_select(chol_package)


def multidirect_select(mlu_package):
    r"""Interface to the mess_direct_select function.


    Parameters
    ----------
    mlu_package:

        - MESS_MULTIDIRECT_SPARSE_LU
        - MESS_MULTIDIRECT_UMFPACK_LU

    Returns
    -------

    Raises
    ------
    RuntimeError if solver is not available

    See Also
    --------

    Notes
    -----

    References
    ----------

    Examples
    --------


    .. doctest::

        >>> from pymess import *
        >>> multidirect_select(MESS_MULTIDIRECT_SPARSE_LU)      # doctest: +SKIP
        >>> multidirect_select(MESS_MULTIDIRECT_UMFPACK_LU)     # doctest: +SKIP

    """

    return _multidirect_select(mlu_package)
