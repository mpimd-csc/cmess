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
from pymess import EquationGRiccatiDAE2, lrnm, Options, MESS_LRCFADI_PARA_ADAPTIVE_V


def main():
    """Solve Riccati Equation DAE2 System."""

    # read data
    m = mmread("@CMAKE_SOURCE_DIR@/tests/data/NSE/NSE_RE_100_lvl1_M.mtx").tocsr()
    a = mmread("@CMAKE_SOURCE_DIR@/tests/data/NSE/NSE_RE_100_lvl1_A.mtx").tocsr()
    g = mmread("@CMAKE_SOURCE_DIR@/tests/data/NSE/NSE_RE_100_lvl1_G.mtx").tocsr()
    b = mmread("@CMAKE_SOURCE_DIR@/tests/data/NSE/NSE_RE_100_lvl1_B.mtx")
    c = mmread("@CMAKE_SOURCE_DIR@/tests/data/NSE/NSE_RE_100_lvl1_C.mtx")
    delta = -0.02

    #create opt instance
    opt = Options()
    opt.adi.output = 0
    opt.nm.output = 0
    opt.nm.res2_tol = 1e-2
    opt.adi.paratype = MESS_LRCFADI_PARA_ADAPTIVE_V

    # create equation
    eqn = EquationGRiccatiDAE2(opt, m, a, g, b, c, delta)

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
