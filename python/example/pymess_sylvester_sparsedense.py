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

"""Example for sylvester_sparsedense."""
from __future__ import print_function, division
from scipy.io import mmread
from pymess import sylvester_sparsedense, res2_sylv_sd

def main():
    """Solve standard/generalized sparse dense Sylvester equation."""

    # read data
    a = mmread("@CMAKE_SOURCE_DIR@/tests/data/filter2D/filter2D.A")
    f = mmread("@CMAKE_SOURCE_DIR@/tests/data/rand10x10_1.mtx")
    e = mmread("@CMAKE_SOURCE_DIR@/tests/data/filter2D/filter2D.E")
    h = mmread("@CMAKE_SOURCE_DIR@/tests/data/rand10x10_2.mtx")
    m = mmread("@CMAKE_SOURCE_DIR@/tests/data/rand10x1668_1.mtx")


    x = sylvester_sparsedense(a, None, None, h, m)
    _, _, rel = res2_sylv_sd(a, None, None, h, m, x)
    print("standard rel_res         =", rel)

    x = sylvester_sparsedense(a, None, e, h, m)
    _, _, rel = res2_sylv_sd(a, None, e, h, m, x)
    print("semi-generalized rel_res =", rel)

    x = sylvester_sparsedense(a, f, e, h, m)
    _, _, rel = res2_sylv_sd(a, f, e, h, m, x)
    print("general rel_res          =", rel)

if __name__ == "__main__":
    main()
