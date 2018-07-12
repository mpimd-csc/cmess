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

"""Test for care function."""
from __future__ import print_function, division
import unittest
from scipy.io import mmread
from scipy.sparse  import csc_matrix, csr_matrix, coo_matrix
import pymess
from pymess import res2_ric, MESS_OP_TRANSPOSE


class Testcare(unittest.TestCase):
    """Test for care function."""
    amtx_rail = '@CMAKE_SOURCE_DIR@/tests/data/Rail/rail_371_c60.A'
    bmtx_rail = '@CMAKE_SOURCE_DIR@/tests/data/Rail/rail_371_c60.B'
    cmtx_rail = '@CMAKE_SOURCE_DIR@/tests/data/Rail/rail_371_c60.C'
    emtx_rail = '@CMAKE_SOURCE_DIR@/tests/data/Rail/rail_371_c60.E'

    amtx_iss = '@CMAKE_SOURCE_DIR@/tests/data/iss/A.mtx'
    bmtx_iss = '@CMAKE_SOURCE_DIR@/tests/data/iss/B.mtx'
    cmtx_iss = '@CMAKE_SOURCE_DIR@/tests/data/iss/C.mtx'
    emtx_iss = '@CMAKE_SOURCE_DIR@/tests/data/iss/E.mtx'

    matformats = [csc_matrix, csr_matrix, coo_matrix]


    def setUp(self):
        """is called before every test function"""
        self.tol = 1e-7

    def template_care(self, hase, mf1, mf2):
        """template to reduce code in the test cases"""
        a = mf1(mmread(self.amtx_rail).tocsr())
        b = mmread(self.bmtx_rail).todense()
        e = mf2(mmread(self.emtx_rail).tocsr()) if hase else None
        c = mmread(self.cmtx_rail).todense()
        #if e is not None:
        #    print("a={0:s} e={1:s} b={2:s} c={3:s}".format(type(a), type(e), type(b), type(c)))
        #else:
        #    print("a={0:s} b={1:s} c={2:s}".format(type(a), type(b), type(c)))

        z = pymess.care(a, e, b, c)
        nrmr, nrmrhs, relnrm = res2_ric(a, e, b, c, z, MESS_OP_TRANSPOSE)
        print("res2 = {0:e}\t rel = {1:e}\t res2/rel = {2:e}".format(nrmr, nrmrhs, relnrm))
        self.assertLessEqual(relnrm, self.tol)

    def template_iss(self):
        """template to reduce code in the test cases"""
        a = mmread(self.amtx_iss).todense()
        b = mmread(self.bmtx_iss).todense()
        e = mmread(self.emtx_iss).todense()
        c = mmread(self.cmtx_iss).todense()
        z = pymess.care(a, e, b, c)
        nrmr, nrmrhs, relnrm = res2_ric(a, e, b, c, z, MESS_OP_TRANSPOSE)
        print("res2 = {0:e}\t rel = {1:e}\t res2/rel = {2:e}".format(nrmr, nrmrhs, relnrm))
        self.assertLessEqual(relnrm, self.tol)


    def test_care_rail(self):
        """test"""
        print("\n")
        for mf1 in self.matformats:
            for mf2 in self.matformats:
                self.template_care(False, mf1, mf2)

    def test_care_generalized_rail(self):
        """test"""
        print("\n")
        for mf1 in self.matformats:
            for mf2 in self.matformats:
                self.template_care(True, mf1, mf2)

    def test_iss(self):
        """test"""
        print("\n")
        self.template_iss()


    def tearDown(self):
        """called after every test function"""
        pass

    def run(self, result=None):
        """ Stop after first error """
        if not result.errors:
            super(Testcare, self).run(result)
