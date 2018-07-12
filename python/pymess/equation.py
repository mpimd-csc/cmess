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

""" The equation class is the base class for all equation objects
    used in Py-M.E.S.S.. Each equation, Lyapunov or Riccati, needs to be derived
    from this class to ensure that the operators are properly defined."""

class Equation(object):
    r""" The equation class is the base class for all equation objects
    used in MESS. Each equation, Lyapunov or Riccati, needs to be derived
    from this class to ensure that the operators are properly defined."""
    b = None
    c = None
    k = None
    rhs = None
    name = ""
    dim = 0
    options = None

    def __init__(self, name="", opt=None, dim=0):
        r""" Constructor of the object. It only sets the name and the options
        to internal variables.  It should be called from the derived class
        using the super mechanism in order to include future changes directly. """
        self.name = name
        self.options = opt
        self.dim = dim

    def __repr__(self):
        r""" Return the representation string. """
        return "Equation: " + self.name

    def __str__(self):
        r""" Return the object as string """
        return "Equation: " + self.name

    # Operator A
    def ax_generate(self):
        r"""Prepares the application of :math:`A` to a right hand side vector or
        matrix. It is called onces before :math:`x = Ay` is called the first time.
        """
        pass

    def ax_clear(self):
        r"""Finalizes the application of the :math:`A` operator. It is called
        after :math:`x = Ay` is applied the last time."""
        pass

    def ax_apply(self, op, y):
        r"""Applies :math:`A` to a right hand side. It has to return
        the result of :math:`x = Ay`. """
        raise NotImplementedError("A_apply needs to be implemented")

    # Operator E
    def ex_generate(self):
        r"""Prepares the application of :math:`E` to a right hand side vector or
        matrix. It is called onces before :math:`x = Ey` is called the first time.
        """
        pass

    def ex_clear(self):
        r"""Finalizes the application of the :math:`E` operator. It is called
        after :math:`x = Ey` is applied the last time."""
        pass

    def ex_apply(self, op, y):
        r"""Applies :math:`E` to a right hand side. It has to return
        the result of :math:`x = Ey`."""
        raise NotImplementedError("ex_apply needs to be implemented")

    # Operator inv(A)
    def ainv_generate(self):
        r"""Prepares the application of the inverse of the :math:`A` operator on a
        given right hand side vector or matrix. It is called onces before :math:`Ax = y` is called the first time.
        """
        pass

    def ainv_clear(self):
        r"""Finalizes the application of the inverse of the :math:`A` operator. It is called
        after :math:`Ax = y` is used the last time."""
        pass

    def ainv_apply(self, op, y):
        r"""Applies the inverse of :math:`A` to a right hand side. It has to return
        the solution of :math:`Ax =y`. """
        raise NotImplementedError("ainv_apply needs to be implemented")

    # Operator inv(E)
    def einv_generate(self):
        r"""Prepares the application of the inverse of the :math:`E` operator on a
        given right hand side vector or matrix. It is called onces before :math:`Ex = y` is called the first time.
        """
        pass

    def einv_clear(self):
        r"""Finalizes the application of the inverse of the :math:`E` operator. It is called
        after :math:`Ex=y` is used the last time."""
        pass

    def einv_apply(self, op, y):
        r"""Applies the inverse of :math:`E` to a right hand side. It has to return
        the solution of :math:`Ex = y`."""
        raise NotImplementedError("einv_apply needs to be implemented")

    # Operator ApE
    def apex_generate(self, p):
        r"""Prepares the application of the :math:`A+pE` operator.
        It is called onces before apex_apply is called the first time."""
        pass

    def apex_clear(self):
        r""" Finalizes the application of the inverse of the :math:`A+pE` operator. It is called
        after apex_apply is used the last time."""
        pass

    def apex_apply(self, op, p, idx_p, y):
        r""" Applies function the operator of :math:`A+pE`. It has to return :math:`x=(A+pE)y`."""
        raise NotImplementedError("apex_apply needs to be implemented")

    # Operator ApE
    def apeinv_generate(self, p):
        r"""Prepares the application of the :math:`(A+pE)^{-1}` operator.
        It is called onces before apeinv_apply is called the first time."""
        pass

    def apeinv_clear(self):
        r""" Finalizes the application of the :math:`(A+pE)^{-1}` operator. It is called
        after apex_apply is used the last time."""
        pass

    def apeinv_apply(self, op, p, idx_p, y):
        r""" Applies the inverse of :math:`A+pE` to a right hand side. It has to return
        the solution of :math:`(A+pE)x = y`."""
        raise NotImplementedError("apeinv_apply needs to be implemented")

    # Additional Helpers
    def parameter(self, arp_p, arp_m, b=None, k=None):
        r""" The parmeter function has to return the shift parameters. If None
        is returned shift paramter will be automatically determined.
        The Shift parameter strategy is determined by the options structure.
        """
        raise NotImplementedError("parameter needs to be implemented")

    def _test(self, t, s):
        return "{0:s} - {1:d} - {2:d}".format(self.name, t, s)

    def __delattr__(self, key):
        raise TypeError(
            "Trying to delete attributes of %r is not allowed" % self)
