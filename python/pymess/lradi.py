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

r"""Interface to mess_lradi."""
from pymess._c_interface import _lradi


# pylint: disable=line-too-long
def lradi(eqn, opt):
    r"""This class represents the interface to the mess_lradi function.

    This functions solves the Lyapunov Equations for different :class:`pymess.equation.Equation` instances.

    Parameters
    ----------
    eqn: equation instance defining Lyapunov Equation
        :class:`~pymess.equation.Equation` instance

    opt: options instance to control algorithm parameters.
        :class:`~pymess.options.Options` instance


    Returns
    -------
    z: (n,s) array_like, Low Rank solution factor.
        :math:`X\approx ZZ^T`.

    stat: status instance, informations about the iteration.
        :class:`~pymess.status.Status` instance

    Raises
    ------

    See Also
    --------
    :class:`pymess.equations.std.EquationGLyap`
    :class:`pymess.equations.so1.EquationGLyapSO1`
    :class:`pymess.equations.so2.EquationGLyapSO2`
    :class:`pymess.equations.dae1.EquationGLyapDAE1`
    :class:`pymess.equations.dae2.EquationGLyapDAE2`
    :func:`pymess.lradi.lradi`

    Notes
    -----
    Supported types are ndarray of numpy and csc, coo or csr matrix of scipy.

    References
    ----------
    :cite:`BenKS12`, :cite:`BenKS13a`, :cite:`BenKS13`, :cite:`BenKS14b`, :cite:`Kue16`

    Examples
    --------

    =================================================== =====================================================================
    Lyapunov Equation                                   Example
    =================================================== =====================================================================
    :class:`~pymess.equations.std.EquationGLyap`        :download:`pymess_lradi.py <../../example/pymess_lradi.py>`
    :class:`~pymess.equations.so1.EquationGLyapSO1`     :download:`pymess_lradi_so1.py <../../example/pymess_lradi_so1.py>`
    :class:`~pymess.equations.so2.EquationGLyapSO2`     :download:`pymess_lradi_so2.py <../../example/pymess_lradi_so2.py>`
    :class:`~pymess.equations.dae1.EquationGLyapDAE1`   :download:`pymess_lradi_dae1.py <../../example/pymess_lradi_dae1.py>`
    :class:`~pymess.equations.dae2.EquationGLyapDAE2`   :download:`pymess_lradi_dae2.py <../../example/pymess_lradi_dae2.py>`
    =================================================== =====================================================================


    """

    return _lradi(eqn, opt)
