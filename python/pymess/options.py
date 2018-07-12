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

r"""This module contains :class:`Options` class, which represents the optoins structure mess_options."""

from scipy import ndarray
from scipy.sparse import isspmatrix
from .enum import *

class Options(object):
    r"""This class represents the options structure mess_options_t.

    Represents the options structure mess_options_t.

    Parameters
    ----------
    type: operation type

    nm: Newton Options
        :class:`.OptionsNewton` instance

    adi: ADI Options
        :class:`.OptionsAdi` instance

    Raises
    ------

    See Also
    --------
    :class:`.OptionsNewton`
    :class:`.OptionsAdi`
    :class:`.OptionsAdiShifts`

    Notes
    -----

    References
    ----------

    Examples
    --------

    .. doctest::

        >>> from pymess import *
        >>> opt = Options()
        >>> opt.type = MESS_OP_NONE
        >>> print(opt)              # doctest: +SKIP
        >>> print(opt.nm)           # doctest: +SKIP
        >>> print(opt.adi)          # doctest: +SKIP


    """

    def __init__(self):
        self.type = MESS_OP_NONE
        self.residual_method = MESS_RESIDUAL_INDEFINITE
        self.nm = OptionsNewton()
        self.adi = OptionsAdi()

    def __repr__(self):
        s = "options class:\n"
        if self.type == MESS_OP_NONE:
            s = s + " " + "type = MESS_OP_NONE\n"
        elif self.type == MESS_OP_TRANSPOSE:
            s = s + " " + "type = MESS_OP_TRANSPOSE\n"
        elif self.type == MESS_OP_HERMITIAN:
            s = s + " " + "type = MESS_OP_HERMITIAN\n"
        else:
            raise TypeError(
                "Unknown operation type of %r is not allowed" % self)

        if self.residual_method == MESS_RESIDUAL_INDEFINITE:
            s = s + " " + "residual_method = MESS_RESIDUAL_INDEFINITE\n"
        elif self.residual_method == MESS_RESIDUAL_SPECTRAL:
            s = s + " " + "residual_method = MESS_RESIDUAL_SPECTRAL\n"
        else:
            raise TypeError(
                "Unknown residual_method type of %r is not allowed" % self)

        s = s + " " + "type = %s\n" % self.type
        s = s + "\n->" + self.nm.__repr__()
        s = s + "\n->" + self.adi.__repr__()
        return s

    def __delattr__(self, key):
        raise TypeError(
            "Trying to delete attributes of %r is not allowed" % self)


class OptionsNewton(object):
    r"""This class represents the options for the newton method from the structure mess_options_t.

    Represents the options for the newton method fom the structure mess_options_t.

    Parameters
    ----------
    type: operation type

    maxit: int
        Maximal number of iteration in the Newton method

    res2_tol: real
        Stopping criteria, 2-norm of relative residual.

    gpstep: int
        Galerkin Projection frequency.

    output: 0 or 1
        Turn output in Newton method on/off.

    singleshifts: 0 or 1
        Turn single shift strategy on or off.

    linesearch: 0 or 1
        Turn line search on or off.

    k0: array_like real matrix
        Initial Feedback.


    Raises
    ------

    See Also
    --------
    :class:`.Options`
    :class:`.OptionsAdi`
    :class:`.OptionsAdiShifts`


    Notes
    -----

    References
    ----------

    Examples
    --------


    """

    def __init__(self):
        self.maxit = 30
        self.res2_tol = 1e-10
        self.gpstep = 0
        self.output = 1
        self.singleshifts = 0
        self.linesearch = 0
        self.k0 = None

    def __repr__(self):
        s = "options_newton class:\n"
        #compute maximal length of options string
        name_max_len = max(
            [len(a) if not a.startswith('_') else 0 for a in dir(self)])
        for attr in dir(self):
            if eval('self.' + attr) is not None:
                if isspmatrix(eval('self.' + attr)) or isinstance(
                        eval('self.' + attr), ndarray):
                    (n, m) = eval('self.' + attr + '.shape')
                    dtype = eval('self.' + attr + '.dtype')
                    s = s + " " + attr + " \t\t:\t[%d x %d] " % (n, m)
                    if isspmatrix(eval('self.' + attr)):
                        s = s + "sparse %s\n" % dtype
                    else:
                        s = s + "dense %s\n" % dtype

                elif attr[0] != '_':
                    name_len = len(attr)
                    #print attributes
                    s = s + " " + attr + " " * (
                        name_max_len - name_len + 1
                    ) + "= %lg\n" % eval('self.' + attr)
        return s

    def __delattr__(self, key):
        raise TypeError(
            "Trying to delete attributes of %r is not allowed" % self)


class OptionsAdi(object):
    r"""This class represents the options for the adi method from the structure mess_options_t.

    Represents the options for the adi method fom the structure mess_options_t.

    Parameters
    ----------
    type: operation type

    maxit: int
        Maximal number of iteration in ADI method.

    res2_tol: real
        relative 2-norm stopping tolerance in ADI method.

    res2c_tol: real
        relative 2-norm change tolerance in ADI method.

    rel_change_tol: real
        relative Frobenius norm change tolerance in the ADI method.

    output: 0 or 1
        Turn output in ADI method on/off.

    memory_usage: memory usage option

    shifts: adi_shifts instance
        :class:`.OptionsAdiShifts` instance


    Raises
    ------

    See Also
    --------
    :class:`.OptionsNewton`
    :class:`.Options`
    :class:`.OptionsAdiShifts`


    Notes
    -----

    References
    ----------

    Examples
    --------


    """

    def __init__(self):
        self.maxit = 500
        self.res2_tol = 1e-10
        self.res2c_tol = 1e-11
        self.rel_change_tol = 1e-10
        self.output = 1
        self.memory_usage = MESS_MEMORY_MID
        self.shifts = OptionsAdiShifts()

    def __repr__(self):
        s = "options_adi class:\n"
        #compute maximal length of options string
        name_max_len = max(
            [len(a) if not a.startswith('_') else 0 for a in dir(self)])
        for attr in dir(self):
            if eval('self.' + attr) is not None:
                if attr == 'shifts':
                    s = s + "\n->" + self.shifts.__repr__()
                elif attr[0] != '_':
                    name_len = len(attr)
                    #print attributes
                    s = s + " " + attr + " " * (
                        name_max_len - name_len + 1
                    ) + "= %lg\n" % eval('self.' + attr)
        return s

    def __delattr__(self, key):
        raise TypeError(
            "Trying to delete attributes of %r is not allowed" % self)


class OptionsAdiShifts(object):
    r"""This class represents the options for the adi method from the structure mess_options_t.

    Represents the options for the adi method fom the structure mess_options_t.

    Parameters
    ----------
    p: precomputed shifts

    arp_p: Number of steps of Arnoldi Process.
        Number of steps of Arnoldi Process  for  :math:`A` or :math:`A-UV`.

    arp_m: Number of steps of Arnoldi Process
        Number of steps of Arnoldi Process for  :math:`A^{-1}` or :math:`(A-UV)^{-1}`.

    paratype: Shift parameter type

    l0: Number of shifts

    b0: Initial vector for Arnoldi process

    Raises
    ------

    See Also
    --------
    :class:`.OptionsNewton`
    :class:`.Options`
    :class:`.OptionsAdiShifts`

    Notes
    -----

    References
    ----------

    Examples
    --------


    """

    def __init__(self):
        self.p = None
        self.arp_p = 48
        self.arp_m = 32
        self.paratype = MESS_LRCFADI_PARA_MINMAX
        self.l0 = 16
        self.b0 = None

    def __repr__(self):
        s = "options_adi_shifts class:\n"
        #compute maximal length of options string
        name_max_len = max(
            [len(a) if not a.startswith('_') else 0 for a in dir(self)])
        for attr in dir(self):
            if not eval('self.' + attr) is None:
                if isinstance(eval('self.' + attr), ndarray):
                    (n, m) = eval('self.' + attr + '.shape')
                    datatype = eval('self.' + attr + '.dtype')
                    s = s + " " + attr + " \t\t:\t[%d x %d] %s\n" % (n, m,
                                                                     datatype)
                elif attr[0] != '_':
                    name_len = len(attr)
                    #print attributes
                    s = s + " " + attr + " " * (
                        name_max_len - name_len + 1
                    ) + "= %lg\n" % eval('self.' + attr)
        return s

    def __delattr__(self, key):
        raise TypeError(
            "Trying to delete attributes of %r is not allowed" % self)
