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

"""Implement a Callback Equation."""
from __future__ import print_function, division
import scipy.sparse as scsp
from scipy.io import mmread
from scipy.sparse.linalg import spsolve
import numpy as np
from pymess import Equation, Options, lrnm, MESS_OP_NONE, MESS_OP_TRANSPOSE, res2_ric



# pylint: disable=unused-argument,no-self-use
class MyEquation(Equation):
    """docstring for MyEquation."""

    def __init__(self, opt, a, e, b, c):
        """Constructor for Riccati Equation"""

        # we have to call the constructor from Equation with dimension of statespace and Options instance
        shapea = a.shape
        dim = shapea[0]
        super(MyEquation, self).__init__(name="MyEquation", opt=opt, dim=dim)

        # mandatory, we have to set the right hand side of the equation our self
        if opt.type == MESS_OP_NONE:
            self.rhs = b
        else:
            self.rhs = c

        # mandatory, we have to set the matrix b and c
        self.b = b
        self.c = c

        # optional, we save the matrix which define our equation, that we can use them in our
        # own function handles
        self.a = scsp.csr_matrix(a)
        self.e = scsp.csc_matrix(e)

    # ***THESE OPERATIONS ARE MANDATORY TO IMPLEMENT***
    def ax_apply(self, op, y):
        """Perform A*y."""
        if op == MESS_OP_NONE:
            return self.a * y
        return self.a.T * y

    def ex_apply(self, op, y):
        """Perform E*y."""
        if op == MESS_OP_NONE:
            return self.e * y
        return self.e.T * y

    def ainv_apply(self, op, y):
        """Perform A^{-1}*y."""
        if op == MESS_OP_NONE:
            return spsolve(self.a, np.array(y))
        return spsolve(self.a.T, np.array(y))

    def einv_apply(self, op, y):
        """Perform E^{-1}*y."""
        if op == MESS_OP_NONE:
            x = spsolve(self.e, np.array(y))
        else:
            x = spsolve(self.e.T, np.array(y))
        return x

    def apex_apply(self, op, p, idx_p, y):
        """Perform (A+pE)*y."""
        if op == MESS_OP_NONE:
            return self.a * y + p * (self.e * y)
        return self.a.T * y + p.conjugate() * (self.e.T * y)

    def apeinv_apply(self, op, p, idx_p, y):
        """Perform (A+pE)^{-1}*y."""
        if op == MESS_OP_NONE:
            tmp = self.a + p * self.e
        elif op == MESS_OP_TRANSPOSE:
            tmp = self.a + p * self.e
            tmp = tmp.T
        else:
            tmp = self.a + p * self.e
            tmp = tmp.H

        return spsolve(tmp, np.array(y))

    def parameter(self, arp_p, arp_m, b=None, k=None):
        """Compute shift parameter. If you return None then C-M.E.S.S. will compute shift parameters for you."""
        return None

     # ** OPERATIONS ARE OPTIONAL TO IMPLEMENT***
    def ax_generate(self):
        """We print something for demonstration."""
        print("--- MyEquation: CALL ax_generate ---")



def main():
    """Demonstrate Callback functionality"""

    # read data
    a = mmread('@CMAKE_SOURCE_DIR@/tests/data/Rail/A.mtx')
    b = mmread('@CMAKE_SOURCE_DIR@/tests/data/Rail/B.mtx')
    c = mmread('@CMAKE_SOURCE_DIR@/tests/data/Rail/C.mtx')
    e = mmread('@CMAKE_SOURCE_DIR@/tests/data/Rail/E.mtx')

    # create options
    opt = Options()
    opt.nm.output = 0
    opt.adi.output = 0

    # create instance of MyEquation and solve
    opt.type = MESS_OP_NONE
    eqn1 = MyEquation(opt, a, e, b, c)
    z1, stat1 = lrnm(eqn1, opt)

    # create instance of MyEquation and solve
    opt.type = MESS_OP_TRANSPOSE
    eqn2 = MyEquation(opt, a, e, b, c)
    z2, stat2 = lrnm(eqn2, opt)

    # print information and compute residual again to make sure that we have everything correct
    print("\n")
    print("MyEquation: eqn1")
    print("Size of Low Rank Solution Factor Z1: %d x %d \n"%(z1.shape))
    print(stat1.lrnm_stat())
    nrmr1, nrmrhs1, relnrm1 = res2_ric(a, e, b, c, z1, MESS_OP_NONE)
    print("check for eqn1:\t rel_res2=%e\t res2=%e\t nrmrhs=%e\n" % (relnrm1, nrmr1, nrmrhs1))

    print("\n")
    print("--------------------------------------------------------------------------------------")
    print("\n")

    # print information and compute residual again to make sure that we have everything correct
    print("MyEquation: eqn2")
    print("Size of Low Rank Solution Factor Z2: %d x %d \n"%(z2.shape))
    print(stat2.lrnm_stat())
    nrmr2, nrmrhs2, relnrm2 = res2_ric(a, e, b, c, z2, MESS_OP_TRANSPOSE)
    print("check for eqn1:\t rel_res2=%e\t res2=%e\t nrmrhs=%e\n" % (relnrm2, nrmr2, nrmrhs2))

if __name__ == "__main__":
    main()
