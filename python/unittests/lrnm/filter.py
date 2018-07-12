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

"""Test lrnm function for filter instance."""
from __future__ import print_function, division
import unittest
from scipy.io import mmread
from pymess import EquationGRiccati, Options, lrnm, MESS_OP_NONE, MESS_OP_TRANSPOSE



class Filter(unittest.TestCase):
    """Test lrnm function for filter instance."""
    amtx = '@CMAKE_SOURCE_DIR@/tests/data/filter2D/filter2D.A'
    bmtx = '@CMAKE_SOURCE_DIR@/tests/data/filter2D/filter2D.B'
    cmtx = '@CMAKE_SOURCE_DIR@/tests/data/filter2D/filter2D.C{0:d}'
    emtx = '@CMAKE_SOURCE_DIR@/tests/data/filter2D/filter2D.E'

    def setUp(self):
        """is called before every test function"""
        self.opt = Options()
        self.opt.adi.output = 0
        self.opt.nm.output = 0
        self.opt.adi.res2_tol = 1e-6
        self.opt.nm.res2_tol = 1e-2
        self.a = mmread(self.amtx).tocsr()
        self.b = mmread(self.bmtx).todense()

    def template_filter(self, mytype, cnum, e):
        """template to reduce code in the test cases"""
        self.opt.type = mytype
        c = mmread(self.cmtx.format(cnum)).todense()
        e = mmread(self.emtx).tocsr() if e else None

        # build equation
        eqn = EquationGRiccati(self.opt, self.a, e, self.b, c)

        #get res2 from lrnm
        _, status = lrnm(eqn, self.opt)
        res2 = status.res2_norm
        res2_0 = status.res2_0
        it = status.it
        print("it=%d\t rel_res2=%e\t res2=%e\t tol=%e\n" %(it, res2 / res2_0, res2, self.opt.nm.res2_tol))
        self.assertLessEqual(res2 / res2_0, self.opt.nm.res2_tol)

    #@unittest.skip("Skip is intended")
    def test_lrnm_riccati_filter_n_c1(self):
        """test"""
        self.template_filter(MESS_OP_NONE, 1, False)

    #@unittest.skip("Skip is intended")
    def test_lrnm_riccati_filter_t_c1(self):
        """test"""
        self.template_filter(MESS_OP_TRANSPOSE, 1, False)

    #@unittest.skip("Skip is intended")
    def test_lrnm_griccati_filter_n_c1(self):
        """test"""
        self.template_filter(MESS_OP_NONE, 1, True)

    #@unittest.skip("Skip is intended")
    def test_lrnm_griccati_filter_t_c1(self):
        """test"""
        self.template_filter(MESS_OP_TRANSPOSE, 1, True)

    def tearDown(self):
        """called after every test function"""
        pass

    def run(self, result=None):
        """ Stop after first error """
        if not result.errors:
            super(Filter, self).run(result)
