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

"""Test lradi function for bips instance, dae1 function handles."""
from __future__ import print_function, division
import unittest
from scipy.io import mmread
from pymess import EquationGLyapDAE1, Options, lradi, MESS_OP_NONE, MESS_OP_TRANSPOSE
from pymess import direct_select, MESS_DIRECT_UMFPACK_LU, MESS_DIRECT_SPARSE_LU


class WWDAE1(unittest.TestCase):
    """Test lradi function for bips instance, dae1 function handles."""
    #MATmtx    = '@CMAKE_SOURCE_DIR@/tests/data/siso_ww_vref_6405/{0:s}.mtx'
    matmtx = '@CMAKE_SOURCE_DIR@/tests/data/bips98_606/{0:s}.mtx'

    def setUp(self):
        """is called before every test function"""
        self.opt = Options()
        self.opt.adi.res2_tol = 1e-10
        self.opt.adi.maxit = 1000
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
        b = mmread(self.matmtx.format('B')).todense()
        c = mmread(self.matmtx.format('C')).todense()

        # build equation
        eqn = EquationGLyapDAE1(self.opt, e11, a11, a12, a21, a22, b if (mytype == MESS_OP_NONE) else c)

        #get res2 from lradi
        _, status = lradi(eqn, self.opt)
        res2 = status.res2_norm
        res2_0 = status.res2_0
        it = status.it
        print("it=%d\t rel_res2=%e\t res2=%e\t tol=%e\n" %(it, res2 / res2_0, res2, self.opt.adi.res2_tol))
        self.assertLessEqual(res2 / res2_0, self.opt.adi.res2_tol)

    def test_lradi_dae1_n_umfpack(self):
        """test"""
        direct_select(MESS_DIRECT_UMFPACK_LU)
        self.template_dae1(MESS_OP_NONE)

    def test_lradi_dae1_n_messlu(self):
        """test"""
        direct_select(MESS_DIRECT_SPARSE_LU)
        self.template_dae1(MESS_OP_NONE)

    def test_lradi_dae1_t_umfpack(self):
        """test"""
        direct_select(MESS_DIRECT_UMFPACK_LU)
        self.template_dae1(MESS_OP_TRANSPOSE)

    def test_lradi_dae1_t_messlu(self):
        """test"""
        direct_select(MESS_DIRECT_SPARSE_LU)
        self.template_dae1(MESS_OP_TRANSPOSE)

    def tearDown(self):
        """called after every test function"""
        pass

    def run(self, result=None):
        """ Stop after first error """
        if not result.errors:
            super(WWDAE1, self).run(result)
