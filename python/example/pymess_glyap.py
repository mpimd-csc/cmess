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

"""Example for glyap and res2_glyap."""
from __future__ import print_function, division
from numpy.random import rand
from pymess import glyap, res2_glyap, MESS_OP_NONE

def main():
    """Solve a generalized Lyapunov Equation using glyap."""

    # generate matrices
    n = 100
    a = rand(n, n)
    e = rand(n, n)
    y1 = rand(n, n)
    y1 = y1 + y1.T
    y2 = rand(n, n)
    y2 = y2 + y2.T

    # solve standard and generalized Lyapunov Equation
    x1, ahat, qa = glyap(a, None, y1, MESS_OP_NONE)
    x2, ahat2, ehat2, qa2, qe2 = glyap(a, e, y1, MESS_OP_NONE)

    # solve again using schur decomposition
    x3 = glyap(a, None, y2, MESS_OP_NONE, ahat, None, qa, None)
    x4 = glyap(a, e, y2, MESS_OP_NONE, ahat2, ehat2, qa2, qe2)

    # compute residual
    res1, _, relres1 = res2_glyap(a, None, y1, x1, MESS_OP_NONE)
    res2, _, relres2 = res2_glyap(a, e, y1, x2, MESS_OP_NONE)
    res3, _, relres3 = res2_glyap(a, None, y2, x3, MESS_OP_NONE)
    res4, _, relres4 = res2_glyap(a, e, y2, x4, MESS_OP_NONE)

    print("n = {0:d}".format(n))
    print("(1): residual = {0:e}\t rel. residual = {1:e}".format(res1, relres1))
    print("(2): residual = {0:e}\t rel. residual = {1:e}".format(res2, relres2))
    print("(3): residual = {0:e}\t rel. residual = {1:e}".format(res3, relres3))
    print("(4): residual = {0:e}\t rel. residual = {1:e}".format(res4, relres4))
    print("\n")


if __name__ == "__main__":
    main()
