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

"""Test glyap and gstein function."""
from __future__ import print_function, division
import unittest
from numpy import sqrt
from numpy.random import rand
from pymess import eps, glyap, res2_glyap, gstein, res2_gstein, MESS_OP_NONE, MESS_OP_TRANSPOSE

class TestGlyap3(unittest.TestCase):
    """Test glyap and gstein function for random instances."""
    ns = [1, 2, 3, 10, 50, 100]
    ops = [MESS_OP_NONE, MESS_OP_TRANSPOSE]
    hase = [False, True]
    funcs = [(glyap, res2_glyap), (gstein, res2_gstein)]
    tol = sqrt(eps())
    verbose = True

    def setUp(self):
        """is called before every test function"""
        print("\n")

    def test_glyap3(self):
        """test"""
        print("\n")
        for fun, res in self.funcs:
            for n in self.ns:
                for op in self.ops:
                    for hase in self.hase:

                        # create random matrices
                        a = rand(n, n)
                        y1 = rand(n, n)
                        y2 = rand(n, n)
                        y1 = y1 + y1.T
                        y2 = y2 + y2.T
                        e = rand(n, n) if hase else None

                        # solve equation and get schur decomposition
                        if hase:
                            x1, ahat, ehat, qa, qe = fun(a, e, y1, op)
                        else:
                            x1, ahat, qa = fun(a, e, y1, op)
                            ehat = None
                            qe = None

                        # solve again using schur decomposition
                        x2 = fun(a, e, y2, op, ahat, ehat, qa, qe)

                        # compute residual and print
                        nrmr1, nrmy1, relnrm1 = res(a, e, y1, x1, op)
                        nrmr2, nrmy2, relnrm2 = res(a, e, y2, x2, op)
                        if self.verbose:
                            print("Funtion: {0}".format(fun))
                            print("hasE = {0}".format(hase))
                            print("op   = {0:d}".format(op))
                            print("n    = {0:d}".format(n))
                            print("res2 = {0:e}\t rel = {1:e}\t res2/rel = {2:e}".format(nrmr1, nrmy1, relnrm1))
                            print("res2 = {0:e}\t rel = {1:e}\t res2/rel = {2:e}".format(nrmr2, nrmy2, relnrm2))
                            print("\n")

                        self.assertLessEqual(relnrm1, self.tol)
                        self.assertLessEqual(relnrm2, self.tol)


    def run(self, result=None):
        """ Stop after first error """
        if not result.errors:
            super(TestGlyap3, self).run(result)
