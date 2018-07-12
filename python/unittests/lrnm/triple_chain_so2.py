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

"""Test lrnm function for triple chain instance, so2 function handles."""
from __future__ import print_function, division
import unittest
from scipy.io import mmread
import numpy as np
from pymess import EquationGRiccatiSO2, Options, lrnm, MESS_OP_NONE, MESS_OP_TRANSPOSE


class TripleChainSO2(unittest.TestCase):
    """Test lrnm function for triple chain instance, so2 function handles."""
    mmtx = '@CMAKE_SOURCE_DIR@/tests/data/TripleChain/M_301.mtx'
    dmtx = '@CMAKE_SOURCE_DIR@/tests/data/TripleChain/D_301.mtx'
    kmtx = '@CMAKE_SOURCE_DIR@/tests/data/TripleChain/K_301.mtx'

    def setUp(self):
        """is called before every test function"""
        self.opt = Options()
        self.opt.adi.res2_tol = 1e-7
        self.opt.nm.res2_tol = 1e-2
        self.opt.adi.output = 0
        self.opt.nm.output = 0
        self.lowerbound = 1e-8
        self.upperbound = 1e+8
        self.m = mmread(self.mmtx).tocsr()
        self.d = mmread(self.dmtx).tocsr()
        self.k = mmread(self.kmtx).tocsr()

    def template_so2(self, mytype):
        """template to reduce code in the test cases"""
        #set type
        self.opt.type = mytype

        #get matrices
        n = self.m.shape[0]
        if self.opt.type == MESS_OP_NONE:
            b = np.random.rand(2 * n, 1)
            c = np.random.rand(1, n)
        else:
            b = np.random.rand(n, 1)
            c = np.random.rand(1, 2 * n)

        #build equation
        eqn = EquationGRiccatiSO2(self.opt, self.m, self.d, self.k, b, c, self.lowerbound, self.upperbound)

        #get res2 from lrnm
        _, status = lrnm(eqn, self.opt)
        res2 = status.res2_norm
        res2_0 = status.res2_0
        it = status.it
        print("it=%d\t rel_res2=%e\t res2=%e\t tol=%e\n" %(it, res2 / res2_0, res2, self.opt.nm.res2_tol))
        self.assertLessEqual(res2 / res2_0, self.opt.nm.res2_tol)

    def test_lrnm_so2_n(self):
        """test"""
        self.template_so2(MESS_OP_NONE)

    def test_lrnm_so2_t(self):
        """test"""
        self.template_so2(MESS_OP_TRANSPOSE)

    def tearDown(self):
        """called after every test function"""
        pass

    def run(self, result=None):
        """ Stop after first error """
        if not result.errors:
            super(TripleChainSO2, self).run(result)
