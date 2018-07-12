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
import numpy as np
from scipy.io import mmread
from pymess import EquationGLyapSO2, lradi, Options, MESS_OP_NONE, MESS_LRCFADI_PARA_MINMAX


def main():
    """Solve Lyapunov Equation SO2 System."""

    # read data
    m = mmread("@CMAKE_SOURCE_DIR@/tests/data/TripleChain/M_301.mtx").tocsr()
    d = mmread("@CMAKE_SOURCE_DIR@/tests/data/TripleChain/D_301.mtx").tocsr()
    k = mmread("@CMAKE_SOURCE_DIR@/tests/data/TripleChain/K_301.mtx").tocsr()

    # generate rhs
    n = m.shape[0]
    b = np.random.rand(2 * n, 1)

    # bounds for shift parameters
    lowerbound = 1e-8
    upperbound = 1e+8

    # create options instance
    opt = Options()
    opt.adi.output = 1
    opt.adi.res2_tol = 1e-7
    opt.type = MESS_OP_NONE
    opt.adi.shifts.paratype = MESS_LRCFADI_PARA_MINMAX

    # create equation
    eqn = EquationGLyapSO2(opt, m, d, k, b, lowerbound, upperbound)

    # solve equation
    z, status = lradi(eqn, opt)

    # get residual
    res2 = status.res2_norm
    res2_0 = status.res2_0
    it = status.it
    print("it = %d \t rel_res2 = %e\t res2 = %e \n" % (it, res2 / res2_0, res2))
    print("Size of Low Rank Solution Factor Z: %d x %d \n"%(z.shape))
    print(status)

if __name__ == "__main__":
    main()
