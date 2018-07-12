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

"""Test lrnm function for bips instance, dae1 function handles."""
from __future__ import print_function, division
import unittest
from scipy.io import mmread
import numpy as np
from pymess import EquationGRiccatiDAE1, Options, lrnm, MESS_OP_NONE, MESS_OP_TRANSPOSE


class WWDAE1(unittest.TestCase):
    """Test lrnm function for bips instance, dae1 function handles."""
    #matmtx    = '@CMAKE_SOURCE_DIR@/tests/data/siso_ww_vref_6405/{0:s}.mtx'
    matmtx = '@CMAKE_SOURCE_DIR@/tests/data/bips98_606/{0:s}.mtx'

    def setUp(self):
        """is called before every test function"""
        self.opt = Options()
        self.opt.adi.res2_tol = 1e-7
        self.opt.adi.maxit = 1500
        self.opt.nm.res2_tol = 1e-2
        self.opt.nm.maxit = 50
        self.opt.adi.output = 0
        self.opt.nm.output = 0

    def template_dae1(self, mytype):
        """template to reduce code in the test cases"""
        #set type
        self.opt.type = mytype

        #get matrices
        e11 = mmread(self.matmtx.format('E11'))
        a11 = mmread(self.matmtx.format('A11'))
        a12 = mmread(self.matmtx.format('A12'))
        a21 = mmread(self.matmtx.format('A21'))
        a22 = mmread(self.matmtx.format('A22'))

        np.random.seed(1111)
        if mytype == MESS_OP_TRANSPOSE:
            b = np.random.rand(e11.shape[0], 2)
            c = np.random.rand(2, e11.shape[1] + a22.shape[1])
        else:
            b = np.random.rand(e11.shape[0] + a22.shape[0], 2)
            c = np.random.rand(2, e11.shape[1])

        b = b / np.linalg.norm(b, 'fro')
        c = c / np.linalg.norm(c, 'fro')

        #build equation
        eqn = EquationGRiccatiDAE1(self.opt, e11, a11, a12, a21, a22, b, c)

        #get res2 from lrnm
        _, status = lrnm(eqn, self.opt)
        res2 = status.res2_norm
        res2_0 = status.res2_0
        it = status.it
        print("it=%d\t rel_res2=%e\t res2=%e\t tol=%e\n" % (it, res2 / res2_0, res2, self.opt.nm.res2_tol))
        self.assertLessEqual(res2 / res2_0, self.opt.nm.res2_tol)

    #@unittest.skip("Skip is intended")
    def test_lrnm_dae1_n(self):
        """test"""
        self.template_dae1(MESS_OP_NONE)

    #@unittest.skip("Skip is intended")
    def test_lrnm_dae1_t(self):
        """test"""
        self.template_dae1(MESS_OP_TRANSPOSE)

    def tearDown(self):
        """called after every test function"""
        pass

    def run(self, result=None):
        """ Stop after first error """
        if not result.errors:
            super(WWDAE1, self).run(result)
