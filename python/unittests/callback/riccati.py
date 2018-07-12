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

"""Test callback for Riccat equation."""
from __future__ import print_function, division
import unittest
from scipy.io import mmread
from numpy import log10
from pymess import MESS_LRCFADI_PARA_MINMAX, MESS_LRCFADI_PARA_MINMAX_REAL
from pymess import MESS_LRCFADI_PARA_WACHSPRESS
from pymess import MESS_LRCFADI_PARA_ADAPTIVE_V, MESS_LRCFADI_PARA_ADAPTIVE_Z
from pymess import MESS_OP_NONE, MESS_OP_TRANSPOSE
from pymess import res2_ric, lrnm, Options
from callback.myequation import MyEquation, MyEquation2


class Riccati(unittest.TestCase):
    """Test callback for Riccat equation."""
    amtx = '@CMAKE_SOURCE_DIR@/tests/data/Rail/A.mtx'
    bmtx = '@CMAKE_SOURCE_DIR@/tests/data/Rail/B.mtx'
    cmtx = '@CMAKE_SOURCE_DIR@/tests/data/Rail/C.mtx'
    emtx = '@CMAKE_SOURCE_DIR@/tests/data/Rail/E.mtx'

    def setUp(self):
        """is called before every test function"""
        self.a = mmread(self.amtx).tocsr()
        self.b = mmread(self.bmtx)
        self.c = mmread(self.cmtx)
        self.e = mmread(self.emtx).tocsr()
        self.opt = Options()
        self.opt.type = MESS_OP_TRANSPOSE

    def template_riccati(self, rhs, para):
        """template to reduce code in the test cases"""
        #set some options
        self.opt.adi.output = 0
        self.opt.adi.res2_tol = 1e-10
        self.opt.nm.output = 0
        self.opt.nm.res2_tol = 1e-3

        # setup equation
        if para < 0:
            eqn = MyEquation2(self.opt, self.a, self.e, rhs, self.b, self.c)
            z, status = lrnm(eqn, self.opt)
        else:
            eqn2 = MyEquation(self.opt, self.a, self.e, rhs, self.b, self.c)
            self.opt.adi.shifts.paratype = para
            z, status = lrnm(eqn2, self.opt)

        #get res2 from lrnm
        res2 = status.res2_norm
        res2_0 = status.res2_0
        it = status.it
        print("it=%d\t rel_res2=%e\t res2=%e\t tol=%e\n" %(it, res2 / res2_0, res2, self.opt.nm.res2_tol))
        self.assertLessEqual(res2 / res2_0, self.opt.nm.res2_tol)

        #check norm
        nrmr, nrmrhs, relnrm = res2_ric(self.a, self.e, self.b, self.c, z, self.opt.type)
        print("check:\t rel_res2=%e\t res2=%e\t nrmrhs=%e\n" % (relnrm, nrmr, nrmrhs))
        # compare residual, should be roughly the same magnitude
        self.assertLessEqual(abs(log10(relnrm)-log10(res2/res2_0)), 1)

    def test_n_own(self):
        """test"""
        self.opt.type = MESS_OP_NONE
        self.template_riccati(self.b, -1)

    def test_t_own(self):
        """test"""
        self.opt.type = MESS_OP_TRANSPOSE
        self.template_riccati(self.c.T, -1)

    def test_n_minmax(self):
        """test"""
        self.opt.type = MESS_OP_NONE
        self.template_riccati(self.b, MESS_LRCFADI_PARA_MINMAX)

    def test_t_minmax(self):
        """test"""
        self.opt.type = MESS_OP_TRANSPOSE
        self.template_riccati(self.c.T, MESS_LRCFADI_PARA_MINMAX)

    def test_n_minmaxreal(self):
        """test"""
        self.opt.type = MESS_OP_NONE
        self.template_riccati(self.b, MESS_LRCFADI_PARA_MINMAX_REAL)

    def test_t_minmaxreal(self):
        """test"""
        self.opt.type = MESS_OP_TRANSPOSE
        self.template_riccati(self.c.T, MESS_LRCFADI_PARA_MINMAX_REAL)

    def test_n_wachspress(self):
        """test"""
        self.opt.type = MESS_OP_NONE
        self.template_riccati(self.b, MESS_LRCFADI_PARA_WACHSPRESS)

    def test_t_wachspress(self):
        """test"""
        self.opt.type = MESS_OP_TRANSPOSE
        self.template_riccati(self.c.T, MESS_LRCFADI_PARA_WACHSPRESS)

    def test_n_adaptive_v(self):
        """test"""
        self.opt.type = MESS_OP_NONE
        self.template_riccati(self.b, MESS_LRCFADI_PARA_ADAPTIVE_V)

    def test_t_adaptive_v(self):
        """test"""
        self.opt.type = MESS_OP_TRANSPOSE
        self.template_riccati(self.c.T, MESS_LRCFADI_PARA_ADAPTIVE_V)

    def test_n_adaptive_z(self):
        """test"""
        self.opt.type = MESS_OP_NONE
        self.template_riccati(self.b, MESS_LRCFADI_PARA_ADAPTIVE_Z)

    def test_t_adaptive_z(self):
        """test"""
        self.opt.type = MESS_OP_TRANSPOSE
        self.template_riccati(self.c.T, MESS_LRCFADI_PARA_ADAPTIVE_Z)

    def tearDown(self):
        """called after every test function"""
        pass

    def run(self, result=None):
        """ Stop after first error """
        if not result.errors:
            super(Riccati, self).run(result)
