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

"""to document"""
from __future__ import print_function, division
from scipy.io import mmread
from numpy.random import seed, rand
from numpy.linalg import norm
from pymess import EquationGRiccatiDAE1, lrnm, Options, MESS_OP_TRANSPOSE

def main():
    """Solve Riccati Equation DAE1 System."""

    # read data
    e11 = mmread("@CMAKE_SOURCE_DIR@/tests/data/bips98_606/E11.mtx").tocsr()
    a11 = mmread("@CMAKE_SOURCE_DIR@/tests/data/bips98_606/A11.mtx").tocsr()
    a12 = mmread("@CMAKE_SOURCE_DIR@/tests/data/bips98_606/A12.mtx").tocsr()
    a21 = mmread("@CMAKE_SOURCE_DIR@/tests/data/bips98_606/A21.mtx").tocsr()
    a22 = mmread("@CMAKE_SOURCE_DIR@/tests/data/bips98_606/A22.mtx").tocsr()

    #create opt instance
    opt = Options()
    opt.adi.output = 0
    opt.nm.output = 0
    opt.nm.res2_tol = 1e-10

    # generate some matrices b and c
    seed(1111)
    if opt.type == MESS_OP_TRANSPOSE:
        b = rand(e11.shape[0], 2)
        c = rand(2, e11.shape[1] + a22.shape[1])
    else:
        b = rand(e11.shape[0] + a22.shape[0], 2)
        c = rand(2, e11.shape[1])

    b = b / norm(b, 'fro')
    c = c / norm(c, 'fro')

    # create equation
    eqn = EquationGRiccatiDAE1(opt, e11, a11, a12, a21, a22, b, c)

    # solve equation
    z, status = lrnm(eqn, opt)

    # get residual
    res2 = status.res2_norm
    res2_0 = status.res2_0
    it = status.it
    print("Size of Low Rank Solution Factor Z: %d x %d \n"%(z.shape))
    print("it = %d \t rel_res2 = %e\t res2 = %e \n" % (it, res2 / res2_0, res2))
    print(status.lrnm_stat())

if __name__ == "__main__":
    main()
