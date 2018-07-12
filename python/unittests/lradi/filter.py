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

"""Test lradi function for filter instance."""
from __future__ import print_function, division
import unittest
from scipy.io import mmread
from pymess import EquationGLyap, lradi, Options, MESS_OP_NONE, MESS_OP_TRANSPOSE
from pymess import multidirect_select, MESS_MULTIDIRECT_UMFPACK_LU, MESS_MULTIDIRECT_SPARSE_LU


class Filter(unittest.TestCase):
    """Test lradi function for filter instance."""
    amtx = '@CMAKE_SOURCE_DIR@/tests/data/filter2D/filter2D.A'
    bmtx = '@CMAKE_SOURCE_DIR@/tests/data/filter2D/filter2D.B'
    cmtx = '@CMAKE_SOURCE_DIR@/tests/data/filter2D/filter2D.C{0:s}'
    emtx = '@CMAKE_SOURCE_DIR@/tests/data/filter2D/filter2D.E'

    def setUp(self):
        """is called before every test function"""
        self.opt = Options()
        self.opt.adi.output = 0
        self.opt.nm.output = 0
        self.opt.adi.res2_tol = 1e-10
        self.a = mmread(self.amtx).tocsr()
        self.b = mmread(self.bmtx).todense()
        self.e = mmread(self.emtx).tocsr()

    def template_filter(self, mytype, cnum, e):
        """template to reduce code in the test cases"""
        self.opt.type = mytype
        if cnum is not None:
            c = mmread(self.cmtx.format(cnum)).todense()

        #build equation
        eqn = EquationGLyap(self.opt, self.a, self.e if e else None, self.b if cnum is None else c)

        #get res2 from lradi
        _, status = lradi(eqn, self.opt)
        res2 = status.res2_norm
        res2_0 = status.res2_0
        it = status.it
        print("it=%d\t rel_res2=%e\t res2=%e\t tol=%e\n" % (it, res2 / res2_0, res2, self.opt.adi.res2_tol))
        self.assertLessEqual(res2 / res2_0, self.opt.adi.res2_tol)

    def test_n_umfpack_multi(self):
        """test"""
        multidirect_select(MESS_MULTIDIRECT_UMFPACK_LU)
        self.template_filter(MESS_OP_NONE, None, False)

    def test_n_messlu_multi(self):
        """test"""
        multidirect_select(MESS_MULTIDIRECT_SPARSE_LU)
        self.template_filter(MESS_OP_NONE, None, False)

    def test_t1_umfpack_multi(self):
        """test"""
        multidirect_select(MESS_MULTIDIRECT_UMFPACK_LU)
        self.template_filter(MESS_OP_TRANSPOSE, '1', False)

    def test_t1_messlu_multi(self):
        """test"""
        multidirect_select(MESS_MULTIDIRECT_SPARSE_LU)
        self.template_filter(MESS_OP_TRANSPOSE, '1', False)

    #@unittest.skip("Skip is intended")
    def test_t2(self):
        """test"""
        self.template_filter(MESS_OP_TRANSPOSE, '2', False)

    #@unittest.skip("Skip is intended")
    def test_t3(self):
        """test"""
        self.template_filter(MESS_OP_TRANSPOSE, '3', False)

    #@unittest.skip("Skip is intended")
    def test_t4(self):
        """test"""
        self.template_filter(MESS_OP_TRANSPOSE, '4', False)

    #@unittest.skip("Skip is intended")
    def test_t5(self):
        """test"""
        self.template_filter(MESS_OP_TRANSPOSE, '5', False)

    def test_gn(self):
        """test"""
        self.template_filter(MESS_OP_NONE, None, True)

    def test_gt1(self):
        """test"""
        self.template_filter(MESS_OP_TRANSPOSE, '1', True)

    #@unittest.skip("Skip is intended")
    def test_gt2(self):
        """test"""
        self.template_filter(MESS_OP_TRANSPOSE, '2', True)

    #@unittest.skip("Skip is intended")
    def test_gt3(self):
        """test"""
        self.template_filter(MESS_OP_TRANSPOSE, '3', True)

    #@unittest.skip("Skip is intended")
    def test_gt4(self):
        """test"""
        self.template_filter(MESS_OP_TRANSPOSE, '4', True)

    #@unittest.skip("Skip is intended")
    def test_gt5(self):
        """test"""
        self.template_filter(MESS_OP_TRANSPOSE, '5', True)

    def tearDown(self):
        """called after every test function"""
        multidirect_select(MESS_MULTIDIRECT_UMFPACK_LU)

    def run(self, result=None):
        """ Stop after first error """
        if not result.errors:
            super(Filter, self).run(result)
