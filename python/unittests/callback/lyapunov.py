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

r"""Test callback for Lyapunov equation."""
from __future__ import print_function, division
import unittest
from numpy import log10
from scipy.io import mmread
from pymess import MESS_LRCFADI_PARA_MINMAX, MESS_LRCFADI_PARA_MINMAX_REAL
from pymess import MESS_LRCFADI_PARA_WACHSPRESS
from pymess import MESS_LRCFADI_PARA_ADAPTIVE_V, MESS_LRCFADI_PARA_ADAPTIVE_Z
from pymess import res2_lyap, lradi, Options, MESS_OP_TRANSPOSE, MESS_OP_NONE
from callback.myequation import MyEquation, MyEquation2


class Lyapunov(unittest.TestCase):
    """Test callback for Lyapunov equation."""
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

    def template_lyapunov(self, rhs, para):
        """template to reduce code in the test cases"""
        #set some options
        self.opt.adi.output = 0
        self.opt.adi.res2_tol = 1e-3

        # setup equation
        if para < 0:
            eqn = MyEquation2(self.opt, self.a, self.e, rhs, None, None)
            z, status = lradi(eqn, self.opt)
        else:
            eqn2 = MyEquation(self.opt, self.a, self.e, rhs, None, None)
            self.opt.adi.shifts.paratype = para
            z, status = lradi(eqn2, self.opt)

        #get res2 from lradi
        res2 = status.res2_norm
        res2_0 = status.res2_0
        it = status.it
        print("it=%d\t rel_res2=%e\t res2=%e\t tol=%e\n" %(it, res2 / res2_0, res2, self.opt.adi.res2_tol))
        self.assertLessEqual(res2 / res2_0, self.opt.adi.res2_tol)

        if self.opt.type == MESS_OP_NONE:
            nrmr, nrmrhs, relnrm = res2_lyap(self.a, self.e, rhs, z, MESS_OP_NONE)
        else:
            nrmr, nrmrhs, relnrm = res2_lyap(self.a, self.e, rhs.T, z, MESS_OP_TRANSPOSE)

        print("check:\t rel_res2=%e\t res2=%e\t nrmrhs=%e\n" % (relnrm, nrmr, nrmrhs))
        # compare residual, should be roughly the same magnitude
        self.assertLessEqual(abs(log10(relnrm)-log10(res2/res2_0)), 1)

    #@unittest.skip("Skip is intended")
    def test_n_own(self):
        """test"""
        self.opt.type = MESS_OP_NONE
        self.template_lyapunov(self.b, -1)

    #@unittest.skip("Skip is intended")
    def test_t_own(self):
        """test"""
        self.opt.type = MESS_OP_TRANSPOSE
        self.template_lyapunov(self.c.T, -1)

    #@unittest.skip("Skip is intended")
    def test_n_minmax(self):
        """test"""
        self.opt.type = MESS_OP_NONE
        self.template_lyapunov(self.b, MESS_LRCFADI_PARA_MINMAX)

    #@unittest.skip("Skip is intended")
    def test_t_minmax(self):
        """test"""
        self.opt.type = MESS_OP_TRANSPOSE
        self.template_lyapunov(self.c.T, MESS_LRCFADI_PARA_MINMAX)

    #@unittest.skip("Skip is intended")
    def test_n_minmaxreal(self):
        """test"""
        self.opt.type = MESS_OP_NONE
        self.template_lyapunov(self.b, MESS_LRCFADI_PARA_MINMAX_REAL)

    #@unittest.skip("Skip is intended")
    def test_t_minmaxreal(self):
        """test"""
        self.opt.type = MESS_OP_TRANSPOSE
        self.template_lyapunov(self.c.T, MESS_LRCFADI_PARA_MINMAX_REAL)

    #@unittest.skip("Skip is intended")
    def test_n_wachspress(self):
        """test"""
        self.opt.type = MESS_OP_NONE
        self.template_lyapunov(self.b, MESS_LRCFADI_PARA_WACHSPRESS)

    #@unittest.skip("Skip is intended")
    def test_t_wachspress(self):
        """test"""
        self.opt.type = MESS_OP_TRANSPOSE
        self.template_lyapunov(self.c.T, MESS_LRCFADI_PARA_WACHSPRESS)

    #@unittest.skip("Skip is intended")
    def test_n_adaptive_v(self):
        """test"""
        self.opt.type = MESS_OP_NONE
        self.template_lyapunov(self.b, MESS_LRCFADI_PARA_ADAPTIVE_V)

    #@unittest.skip("Skip is intended")
    def test_t_adaptive_v(self):
        """test"""
        self.opt.type = MESS_OP_TRANSPOSE
        self.template_lyapunov(self.c.T, MESS_LRCFADI_PARA_ADAPTIVE_V)

    #@unittest.skip("Skip is intended")
    def test_n_adaptive_z(self):
        """test"""
        self.opt.type = MESS_OP_NONE
        self.template_lyapunov(self.b, MESS_LRCFADI_PARA_ADAPTIVE_Z)

    #@unittest.skip("Skip is intended")
    def test_t_adaptive_z(self):
        """test"""
        self.opt.type = MESS_OP_TRANSPOSE
        self.template_lyapunov(self.c.T, MESS_LRCFADI_PARA_ADAPTIVE_Z)

    def tearDown(self):
        """called after every test function"""
        pass

    def run(self, result=None):
        """ Stop after first error """
        if not result.errors:
            super(Lyapunov, self).run(result)
