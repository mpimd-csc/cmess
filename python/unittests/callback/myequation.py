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

"""Implement  a Callback Equation."""
from __future__ import print_function, division
import scipy.sparse as scsp
from scipy.sparse.linalg import spsolve
from scipy.sparse.linalg import LinearOperator
from scipy.sparse.linalg import eigs
import numpy as np
from pymess import Equation, MESS_OP_NONE, MESS_OP_TRANSPOSE, MESSErrorArguments

# pylint: disable=unused-argument,no-self-use
class MyEquation(Equation):
    """docstring for MyEquation."""

    def __init__(self, opt, a, e, rhs, b, c):
        shapea = a.shape
        dim = shapea[0]
        super(MyEquation, self).__init__(name="MyEquation", opt=opt, dim=dim)
        self.rhs = rhs
        self.b = b
        self.c = c
        self.a = scsp.csc_matrix(a)
        self.e = scsp.csc_matrix(e)
        self.p = []

    # operator A
    def ax_apply(self, op, y):
        """to document"""
        if op == MESS_OP_NONE:
            return self.a * y
        return self.a.T * y

    # operator E
    def ex_apply(self, op, y):
        """to document"""
        if op == MESS_OP_NONE:
            return self.e * y
        return self.e.T * y

    # operator A^{-1}
    def ainv_apply(self, op, y):
        """to document"""
        if op == MESS_OP_NONE:
            return spsolve(self.a, np.array(y))
        return spsolve(self.a.T, np.array(y))

    # operator E^{-1}
    def einv_apply(self, op, y):
        """to document"""
        if op == MESS_OP_NONE:
            return spsolve(self.e, np.array(y))
        return spsolve(self.e.T, np.array(y))

    # operator (A+pE)
    def apex_apply(self, op, p, idx_p, y):
        """to document"""
        if op == MESS_OP_NONE:
            return self.a * y + p * (self.e * y)
        return self.a.T * y + p.conjugate() * (self.e.T * y)

    # operator (A+pE)^{-1}
    def apeinv_apply(self, op, p, idx_p, y):
        """to document"""
        if op == MESS_OP_NONE:
            tmp = self.a + p * self.e
        elif op == MESS_OP_TRANSPOSE:
            tmp = self.a + p * self.e
            tmp = tmp.T
        else:
            tmp = self.a + p * self.e
            tmp = tmp.H

        return spsolve(tmp, np.array(y))

    # parameter function
    def parameter(self, arp_p, arp_m, b=None, k=None):
        """to document"""
        return None


class MyEquation2(MyEquation):
    """docstring for MyEquation2."""

    def parameter(self, arp_p, arp_m, b=None, k=None):
        """to document"""
        if (b is not None and k is None) or (b is None and k is not None):
            raise MESSErrorArguments("You have to provide B and K both.")
        if b is None and k is None:
            lm = eigs(self.a, 20, self.e, which='LM', return_eigenvectors=False)
            sm = eigs(self.a, 20, self.e, which='SM', return_eigenvectors=False)
        else:
            #we build a stabilized linear operator and compute the generalized eigenvalues
            def mv_bk(v):
                """to document"""
                v = np.matrix(v).T
                return self.a.dot(v) - b * (k.dot(v))

            def rmv_bk(v):
                """to document"""
                v = np.matrix(v).T
                return self.a.H.dot(v) - k.H * (b.H.dot(v))

            def mm_bk(v):
                """to document"""
                return self.a * v - b * (k * v)

            abk = LinearOperator(self.a.shape, matvec=mv_bk, rmatvec=rmv_bk, matmat=mm_bk, dtype=self.a.dtype)
            lm = eigs(abk, 20, self.e, which='LM', return_eigenvectors=False)
            sm = eigs(abk, 20, self.e, which='SM', return_eigenvectors=False)

        #concatenate both and take real part, and filter positive ones
        ev = np.concatenate((lm, sm))
        ev = ev.real
        ev = ev[ev < 0]
        return ev
