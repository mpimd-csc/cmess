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

"""Test for sylvester_sparsedense function."""
from __future__ import print_function, division
import unittest
import numpy as np
from scipy.io import mmread
from scipy.sparse  import csc_matrix, csr_matrix, coo_matrix
import pymess
from pymess import res2_sylv_sd

class Testsylvestersparsedense(unittest.TestCase):
    """Test for sylvester_sparsedense function."""
    amtx = '@CMAKE_SOURCE_DIR@/tests/data/Rail/A.mtx'
    fmtx = '@CMAKE_SOURCE_DIR@/tests/data/rand10x10_1.mtx'
    emtx = '@CMAKE_SOURCE_DIR@/tests/data/Rail/E.mtx'
    hmtx = '@CMAKE_SOURCE_DIR@/tests/data/rand10x10_2.mtx'
    mmtx = '@CMAKE_SOURCE_DIR@/tests/data/rand1357x10.mtx'
    matformats = [csc_matrix, csr_matrix, coo_matrix]

    def setUp(self):
        """is called before every test function"""
        self.a = mmread(self.amtx).tocsc()
        self.f = mmread(self.fmtx)
        self.e = mmread(self.emtx).tocsc()
        self.h = mmread(self.hmtx)
        self.m = mmread(self.mmtx)
        self.tol = 1e-10

    def template_sylvester_sparsedense(self, n, hase, hasf):
        """template to reduce code in the test cases"""
        #solve sylvester equation
        a = self.a
        if hasf and not hase:
            print("ignoring input hasf=True because there is no e (hase=False)")
        f = np.random.rand(n, n) if hase and hasf else None
        e = self.e if hase else None
        h = np.random.rand(n, n)
        m = np.random.rand(1357, n)

        x = pymess.sylvester_sparsedense(a, f, e, h, m)
        nrmr, nrmrhs, relnrm = res2_sylv_sd(a, f, e, h, m, x)
        print("res2 = {0:e}\t nrmrhs = {1:e}\t rel = res2/nrmrhs = {2:e}".format(nrmr, nrmrhs, relnrm))
        self.assertLessEqual(relnrm, self.tol)

    def template_sylvester_sparsedense_fixed(self, hase, hasf):
        """template to reduce code in the test cases"""
        #solve sylvester equation
        a = self.a
        if hasf and not hase:
            print("ignoring input hasf=True because there is no e (hase=False)")
        f = self.f if hase and hasf else None
        e = self.e if hase else None
        h = self.h
        m = self.m

        x = pymess.sylvester_sparsedense(a, f, e, h, m)
        nrmr, nrmrhs, relnrm = res2_sylv_sd(a, f, e, h, m, x)
        print("res2 = {0:e}\t nrmrhs = {1:e}\t rel = res2/nrmrhs = {2:e}".format(nrmr, nrmrhs, relnrm))
        self.assertLessEqual(relnrm, self.tol)

    def test_sylvester_sparsedense_nongeneralized_filter(self):
        """test"""
        for mf1 in self.matformats:
            for mf2 in self.matformats:
                self.a = mf1(self.a)
                self.e = mf2(self.e)
                self.template_sylvester_sparsedense_fixed(False, False)

    def test_sylvester_sparsedense_semigeneralized_filter(self):
        """test"""
        for mf1 in self.matformats:
            for mf2 in self.matformats:
                self.a = mf1(self.a)
                self.e = mf2(self.e)
                self.template_sylvester_sparsedense_fixed(True, False)

    def test_sylvester_sparsedense_generalized_filter(self):
        """test"""
        for mf1 in self.matformats:
            for mf2 in self.matformats:
                self.a = mf1(self.a)
                self.e = mf2(self.e)
                self.template_sylvester_sparsedense_fixed(True, True)


    def test_ngsylvester_sparsedense_1(self):
        """test"""
        print("\n")
        self.template_sylvester_sparsedense(1, 0, 0)

    def test_ngsylvester_sparsedense_2(self):
        """test"""
        print("\n")
        self.template_sylvester_sparsedense(2, 0, 0)

    def test_ngsylvester_sparsedense_5(self):
        """test"""
        print("\n")
        self.template_sylvester_sparsedense(5, 0, 0)

    def test_ngsylvester_sparsedense_10(self):
        """test"""
        print("\n")
        self.template_sylvester_sparsedense(10, 0, 0)

    def test_ngsylvester_sparsedense_25(self):
        """test"""
        print("\n")
        self.template_sylvester_sparsedense(25, 0, 0)

    def test_ngsylvester_sparsedense_33(self):
        """test"""
        print("\n")
        self.template_sylvester_sparsedense(33, 0, 0)

    def test_ngsylvester_sparsedense_100(self):
        """test"""
        print("\n")
        self.template_sylvester_sparsedense(100, 0, 0)

    def test_sgsylvester_sparsedense_1(self):
        """test"""
        print("\n")
        self.template_sylvester_sparsedense(1, 1, 0)

    def test_sgsylvester_sparsedense_2(self):
        """test"""
        print("\n")
        self.template_sylvester_sparsedense(2, 1, 0)

    def test_sgsylvester_sparsedense_5(self):
        """test"""
        print("\n")
        self.template_sylvester_sparsedense(5, 1, 0)

    def test_sgsylvester_sparsedense_10(self):
        """test"""
        print("\n")
        self.template_sylvester_sparsedense(10, 1, 0)

    def test_sgsylvester_sparsedense_25(self):
        """test"""
        print("\n")
        self.template_sylvester_sparsedense(25, 1, 0)

    def test_sgsylvester_sparsedense_33(self):
        """test"""
        print("\n")
        self.template_sylvester_sparsedense(33, 1, 0)

    def test_sgsylvester_sparsedense_100(self):
        """test"""
        print("\n")
        self.template_sylvester_sparsedense(100, 1, 0)

    def test_gsylvester_sparsedense_1(self):
        """test"""
        print("\n")
        self.template_sylvester_sparsedense(1, 1, 1)

    def test_gsylvester_sparsedense_2(self):
        """test"""
        print("\n")
        self.template_sylvester_sparsedense(2, 1, 1)

    def test_gsylvester_sparsedense_5(self):
        """test"""
        print("\n")
        self.template_sylvester_sparsedense(5, 1, 1)

    def test_gsylvester_sparsedense_10(self):
        """test"""
        print("\n")
        self.template_sylvester_sparsedense(10, 1, 1)

    def test_gsylvester_sparsedense_25(self):
        """test"""
        print("\n")
        self.template_sylvester_sparsedense(25, 1, 1)

    def test_gsylvester_sparsedense_33(self):
        """test"""
        print("\n")
        self.template_sylvester_sparsedense(33, 1, 1)

    def test_gsylvester_sparsedense_100(self):
        """test"""
        print("\n")
        self.template_sylvester_sparsedense(100, 1, 1)

    def tearDown(self):
        """called after every test function"""
        pass

    def run(self, result=None):
        """ Stop after first error """
        if not result.errors:
            super(Testsylvestersparsedense, self).run(result)
