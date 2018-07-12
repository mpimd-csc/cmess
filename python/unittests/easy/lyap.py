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

"""Test for lyap function."""
from __future__ import print_function, division
import unittest
from scipy.io import mmread
from scipy.sparse  import csc_matrix, csr_matrix, coo_matrix
import numpy as np
from pymess import res2_lyap, MESS_OP_NONE
import pymess

class Testlyap(unittest.TestCase):
    """Test for lyap function."""
    amtx = '@CMAKE_SOURCE_DIR@/tests/data/Rail/A.mtx'
    bmtx = '@CMAKE_SOURCE_DIR@/tests/data/Rail/B.mtx'
    emtx = '@CMAKE_SOURCE_DIR@/tests/data/Rail/E.mtx'
    matformats = [csc_matrix, csr_matrix, coo_matrix]

    def setUp(self):
        """is called before every test function"""
        self.a = mmread(self.amtx).tocsr()
        self.b = mmread(self.bmtx)
        self.e = mmread(self.emtx).tocsr()
        self.tol = 1e-10

    def template_lyap(self, hase):
        """template to reduce code in the test cases"""
        #solve lyapunov equation
        a = self.a
        e = self.e if hase else None
        b = self.b
        z = pymess.lyap(a, e, b)
        nrmr, nrmrhs, relnrm = res2_lyap(a, e, b, z, MESS_OP_NONE)
        print("res2 = {0:e}\t rel = {1:e}\t res2/rel = {2:e}".format(nrmr, nrmrhs, relnrm))
        self.assertLessEqual(relnrm, self.tol)

    def test_lyap_filter(self):
        """test"""
        for mf1 in self.matformats:
            for mf2 in self.matformats:
                self.a = mf1(self.a)
                self.e = mf2(self.e)
                self.template_lyap(False)

    def test_lyap_generalized_filter(self):
        """test"""
        for mf1 in self.matformats:
            for mf2 in self.matformats:
                self.a = mf1(self.a)
                self.e = mf2(self.e)
                self.template_lyap(True)

    def template_lyap_dense(self, n, hase):
        """template to reduce code in test cases"""

        # generate stable matrix A
        d = np.diag(np.ravel(-10 * np.random.rand(n, 1) - 10))
        q, _ = np.linalg.qr(np.random.rand(n, n))
        a = q.dot(d).dot(q.T)

        # generate B
        b = np.random.rand(n, np.max([1, np.int(n / 10)]))

        if hase:
            # generate matrix E
            e = np.diag(np.ravel(np.random.rand(n, 1)))

            # solve lyapunov
            z = pymess.lyap(a, e, b)

            # residual
            nrmr, nrmrhs, relnrm = res2_lyap(a, e, b, z, MESS_OP_NONE)

        else:
            # solve lyapunov
            z = pymess.lyap(a, None, b)

            # residual
            nrmr, nrmrhs, relnrm = res2_lyap(a, None, b, z, MESS_OP_NONE)

        print("res2 = {0:e}\t rel = {1:e}\t res2/rel = {2:e}".format(nrmr, nrmrhs, relnrm))
        self.assertLessEqual(relnrm, self.tol)

    def test_lyap_dense_1(self):
        """test"""
        print("\n")
        self.template_lyap_dense(1, 0)

    def test_lyap_dense_2(self):
        """test"""
        print("\n")
        self.template_lyap_dense(2, 0)

    def test_lyap_dense_5(self):
        """test"""
        print("\n")
        self.template_lyap_dense(5, 0)

    def test_lyap_dense_10(self):
        """test"""
        print("\n")
        self.template_lyap_dense(10, 0)

    def test_lyap_dense_100(self):
        """test"""
        print("\n")
        self.template_lyap_dense(100, 0)

    def test_lyap_dense_250(self):
        """test"""
        print("\n")
        self.template_lyap_dense(250, 0)

    def test_lyap_dense_333(self):
        """test"""
        print("\n")
        self.template_lyap_dense(333, 0)

    def test_glyap_dense_1(self):
        """test"""
        print("\n")
        self.template_lyap_dense(1, 1)

    def test_glyap_dense_2(self):
        """test"""
        print("\n")
        self.template_lyap_dense(2, 1)

    def test_glyap_dense_5(self):
        """test"""
        print("\n")
        self.template_lyap_dense(5, 1)

    def test_glyap_dense_10(self):
        """test"""
        print("\n")
        self.template_lyap_dense(10, 1)

    def test_glyap_dense_100(self):
        """test"""
        print("\n")
        self.template_lyap_dense(100, 1)

    def test_glyap_dense_250(self):
        """test"""
        print("\n")
        self.template_lyap_dense(250, 1)

    def test_glyap_dense_333(self):
        """test"""
        print("\n")
        self.template_lyap_dense(333, 1)

    def tearDown(self):
        """called after every test function"""
        pass

    def run(self, result=None):
        """ Stop after first error """
        if not result.errors:
            super(Testlyap, self).run(result)
