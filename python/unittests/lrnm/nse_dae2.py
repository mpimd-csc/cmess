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

"""Test lrnm function for navier stokes instance."""
from __future__ import print_function, division
import unittest
from scipy.io import mmread
from pymess import EquationGRiccatiDAE2, Options, lrnm, MESS_OP_NONE, MESS_OP_TRANSPOSE, MESS_MEMORY_LOW


class NSEDAE2(unittest.TestCase):
    """Test lrnm function for navier stokes instance."""
    matmtx = '@CMAKE_SOURCE_DIR@/tests/data/NSE/NSE_RE_{0:d}_lvl1_{1:s}.mtx'

    def setUp(self):
        """is called before every test function"""
        self.opt = Options()
        self.opt.adi.res2_tol = 1e-6
        self.opt.nm.res2_tol = 5e-1
        self.opt.nm.maxit = 30
        self.opt.adi.output = 0
        self.opt.nm.output = 0
        self.delta = -0.02

    def template_dae2(self, mytype, re, mem_usage):
        """template to reduce code in the test cases"""
        #set type
        self.opt.type = mytype
        self.opt.adi.memory_usage = mem_usage

        #get matrices
        m = mmread(self.matmtx.format(re, 'M'))
        a = mmread(self.matmtx.format(re, 'A'))
        g = mmread(self.matmtx.format(re, 'G'))
        b = mmread(self.matmtx.format(re, 'B'))
        c = mmread(self.matmtx.format(re, 'C'))

        if self.opt.type == MESS_OP_NONE:
            if re >= 300:
                self.opt.nm.k0 = mmread(self.matmtx.format(re, 'Feed1'))
        else:
            if re >= 300:
                self.opt.nm.k0 = mmread(self.matmtx.format(re, 'Feed0'))

        #build equation
        eqn = EquationGRiccatiDAE2(self.opt, m, a, g, b, c, self.delta)

        #get res2 from lrnm
        _, status = lrnm(eqn, self.opt)
        res2 = status.res2_norm
        res2_0 = status.res2_0
        it = status.it
        print("it=%d\t rel_res2=%e\t res2=%e\t tol=%e\n" %(it, res2 / res2_0, res2, self.opt.nm.res2_tol))
        self.assertLessEqual(res2 / res2_0, self.opt.nm.res2_tol)

    #@unittest.skip("Skip is intended")
    def test_n_re100_mem0(self):
        """test"""
        self.template_dae2(MESS_OP_NONE, 100, MESS_MEMORY_LOW)

    #@unittest.skip("Skip is intended")
    def test_t_re100_mem0(self):
        """test"""
        self.template_dae2(MESS_OP_TRANSPOSE, 100, MESS_MEMORY_LOW)

    #@unittest.skip("Skip is intended")
    def test_n_feed_re300_mem0(self):
        """test"""
        self.template_dae2(MESS_OP_NONE, 300, MESS_MEMORY_LOW)

    #@unittest.skip("Skip is intended")
    def test_t_feed_re300_mem0(self):
        """test"""
        self.template_dae2(MESS_OP_TRANSPOSE, 300, MESS_MEMORY_LOW)

    def tearDown(self):
        """called after every test function"""
        pass


    def run(self, result=None):
        """ Stop after first error """
        if not result.errors:
            super(NSEDAE2, self).run(result)
