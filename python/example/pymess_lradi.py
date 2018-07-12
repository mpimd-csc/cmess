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

""" to document"""
from __future__ import print_function, division
from scipy.io import mmread
from pymess import EquationGLyap, Options, lradi


def main():
    """Solve standard/generalized Lyapunov equation."""

    # read data
    e = mmread("@CMAKE_SOURCE_DIR@/tests/data/filter2D/filter2D.E").tocsr()
    a = mmread("@CMAKE_SOURCE_DIR@/tests/data/filter2D/filter2D.A").tocsr()
    b = mmread("@CMAKE_SOURCE_DIR@/tests/data/filter2D/filter2D.B")

    # create opt instances
    opt1 = Options()
    opt1.adi.output = 0
    opt1.adi.res2_tol = 1e-6
    print(opt1)
    opt2 = Options()
    opt2.adi.output = 0
    opt2.adi.res2_tol = 1e-6
    print(opt2)

    # create standard and generalized equation
    eqn1 = EquationGLyap(opt1, a, None, b)
    eqn2 = EquationGLyap(opt2, a, e, b)

    # solve both equations
    z1, stat1 = lradi(eqn1, opt1)
    z2, stat2 = lradi(eqn2, opt2)

    print("Standard Lyapunov Equation \n")
    print("Size of Low Rank Solution Factor Z1: %d x %d \n"%(z1.shape))
    print(stat1)
    print("Generalized Lyapunov Equation \n")
    print("Size of Low Rank Solution Factor Z2: %d x %d \n"%(z2.shape))
    print(stat2)

if __name__ == "__main__":
    main()
