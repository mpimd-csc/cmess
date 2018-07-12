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

"""Example for dense_nm_gmpare and res2_gmpare."""
from __future__ import print_function, division
from scipy.io import mmread
from pymess import dense_nm_gmpare, res2_gmpare, MESS_OP_TRANSPOSE

def main():
    """Solve a generalized dense positvie negative Algebraic Riccati equation using dense_nm_gmpare."""

    # read data
    a = mmread("@CMAKE_SOURCE_DIR@/tests/data/iss/A.mtx").todense()
    e = mmread("@CMAKE_SOURCE_DIR@/tests/data/iss/E.mtx").todense()
    b = mmread("@CMAKE_SOURCE_DIR@/tests/data/iss/B.mtx").todense()
    c = mmread("@CMAKE_SOURCE_DIR@/tests/data/iss/C.mtx").todense()

    g = c.T*c
    q = b*b.T
    n = a.shape[0]

    # solve generalized riccati equation with minus
    xm, m_absres_c, m_relres_c = dense_nm_gmpare(None, a, e, q, g, plus=False, trans=MESS_OP_TRANSPOSE)

    # solve generalized riccati equation with plus
    xp, p_absres_c, p_relres_c = dense_nm_gmpare(None, a, e, q, g, plus=True, trans=MESS_OP_TRANSPOSE)

    # compute residuals
    m_absres, _, m_relres = res2_gmpare(a, e, q, g, xm, plus=False, trans=MESS_OP_TRANSPOSE)
    p_absres, _, p_relres = res2_gmpare(a, e, q, g, xp, plus=True, trans=MESS_OP_TRANSPOSE)

    print("n = {0:d}".format(n))
    print("Negative Riccati Equation")
    print("dense_nm_gmpare: abs./rel. 2-Norm residual = {0:e} / {1:e}".format(m_absres_c, m_relres_c))
    print("res2_gmpare:     abs./rel. 2-Norm residual = {0:e} / {1:e}".format(m_absres, m_relres))
    print("\n")
    print("Positive Riccati Equation")
    print("dense_nm_gmpare: abs./rel. 2-Norm residual = {0:e} / {1:e}".format(p_absres_c, p_relres_c))
    print("res2_gmpare:     abs./rel. 2-Norm residual = {0:e} / {1:e}".format(p_absres, p_relres))
    print("\n")

if __name__ == "__main__":
    main()
