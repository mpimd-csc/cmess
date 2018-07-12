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
from numpy.linalg import norm
from pymess import EquationGLyapDAE1, lradi, Options

def main():
    """Solve Lyapunov Equation DAE1 System."""

    # read data
    e11 = mmread("@CMAKE_SOURCE_DIR@/tests/data/bips98_606/E11.mtx").tocsr()
    a11 = mmread("@CMAKE_SOURCE_DIR@/tests/data/bips98_606/A11.mtx").tocsr()
    a12 = mmread("@CMAKE_SOURCE_DIR@/tests/data/bips98_606/A12.mtx").tocsr()
    a21 = mmread("@CMAKE_SOURCE_DIR@/tests/data/bips98_606/A21.mtx").tocsr()
    a22 = mmread("@CMAKE_SOURCE_DIR@/tests/data/bips98_606/A22.mtx").tocsr()
    b = mmread("@CMAKE_SOURCE_DIR@/tests/data/bips98_606/B.mtx").todense()
    b = b / norm(b)

    #create opt instance
    opt = Options()
    opt.adi.output = 0
    opt.nm.output = 0
    opt.adi.res2_tol = 1e-10
    print(opt)

    # create equation
    eqn = EquationGLyapDAE1(opt, e11, a11, a12, a21, a22, b)

    # solve equation
    z, status = lradi(eqn, opt)

    # get residual
    res2 = status.res2_norm
    res2_0 = status.res2_0
    it = status.it
    print("Size of Low Rank Solution Factor Z: %d x %d \n"%(z.shape))
    print("it = %d \t rel_res2 = %e\t res2 = %e \n" % (it, res2 / res2_0, res2))
    print(status)

if __name__ == "__main__":
    main()
