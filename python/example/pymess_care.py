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

"""Example for care."""
from __future__ import print_function, division
from scipy.io import mmread
from pymess import res2_ric, MESS_OP_TRANSPOSE, care

def main():
    """Solve standard/generalized Riccati equation."""

    # read data
    a = mmread("@CMAKE_SOURCE_DIR@/tests/data/filter2D/filter2D.A")
    e = mmread("@CMAKE_SOURCE_DIR@/tests/data/filter2D/filter2D.E")
    b = mmread("@CMAKE_SOURCE_DIR@/tests/data/filter2D/filter2D.B").todense()
    c = mmread("@CMAKE_SOURCE_DIR@/tests/data/filter2D/filter2D.C").todense()

    # standard
    z = care(a, None, b, c)

    # compute residual
    _, _, rel = res2_ric(a, None, b, c, z, MESS_OP_TRANSPOSE)
    print("standard rel. resisudal=", rel)

    # generalized
    z = care(a, e, b, c)

    # compute residual
    _, _, rel = res2_ric(a, e, b, c, z, MESS_OP_TRANSPOSE)
    print("generalized rel. resisdual=", rel)

if __name__ == "__main__":
    main()
