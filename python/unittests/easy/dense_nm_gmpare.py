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

"""Test for dense_nm_gmpare function."""
from __future__ import print_function, division
import unittest
from scipy.io import mmread
from pymess import res2_gmpare, dense_nm_gmpare, MESS_OP_NONE, MESS_OP_TRANSPOSE, MESSErrorArguments


class Testdensenmgmpare(unittest.TestCase):
    """Test for dense_nm_gmpare function."""
    amtx_rail = '@CMAKE_SOURCE_DIR@/tests/data/Rail/rail_371_c60.A'
    bmtx_rail = '@CMAKE_SOURCE_DIR@/tests/data/Rail/rail_371_c60.B'
    cmtx_rail = '@CMAKE_SOURCE_DIR@/tests/data/Rail/rail_371_c60.C'
    emtx_rail = '@CMAKE_SOURCE_DIR@/tests/data/Rail/rail_371_c60.E'

    amtx_iss = '@CMAKE_SOURCE_DIR@/tests/data/iss/A.mtx'
    bmtx_iss = '@CMAKE_SOURCE_DIR@/tests/data/iss/B.mtx'
    cmtx_iss = '@CMAKE_SOURCE_DIR@/tests/data/iss/C.mtx'
    emtx_iss = '@CMAKE_SOURCE_DIR@/tests/data/iss/E.mtx'

    amtx_heat_cont = '@CMAKE_SOURCE_DIR@/tests/data/heat-cont/A.mtx'
    bmtx_heat_cont = '@CMAKE_SOURCE_DIR@/tests/data/heat-cont/B.mtx'
    cmtx_heat_cont = '@CMAKE_SOURCE_DIR@/tests/data/heat-cont/C.mtx'
    emtx_heat_cont = '@CMAKE_SOURCE_DIR@/tests/data/heat-cont/E.mtx'

    amtx_build = '@CMAKE_SOURCE_DIR@/tests/data/build/A.mtx'
    bmtx_build = '@CMAKE_SOURCE_DIR@/tests/data/build/B.mtx'
    cmtx_build = '@CMAKE_SOURCE_DIR@/tests/data/build/C.mtx'
    emtx_build = '@CMAKE_SOURCE_DIR@/tests/data/build/E.mtx'


    def setUp(self):
        """is called before every test function"""

        # setup rail matrices
        self.a_rail = mmread(self.amtx_rail).todense()
        self.e_rail = mmread(self.emtx_rail).todense()
        b_rail = mmread(self.bmtx_rail).todense()
        c_rail = mmread(self.cmtx_rail).todense()
        self.q_rail = c_rail.T.dot(c_rail)
        self.g_rail = b_rail.dot(b_rail.T)

        # setup iss matrices
        self.a_iss = mmread(self.amtx_iss).todense()
        self.e_iss = mmread(self.emtx_iss).todense()
        b_iss = mmread(self.bmtx_iss).todense()
        c_iss = mmread(self.cmtx_iss).todense()
        self.q_iss = c_iss.T.dot(c_iss)
        self.g_iss = b_iss.dot(b_iss.T)

        # setup heat matrices
        self.a_heat_cont = mmread(self.amtx_heat_cont).todense()
        self.e_heat_cont = mmread(self.emtx_heat_cont).todense()
        b_heat_cont = mmread(self.bmtx_heat_cont)
        c_heat_cont = mmread(self.cmtx_heat_cont)
        self.q_heat_cont = c_heat_cont.T.dot(c_heat_cont)
        self.g_heat_cont = b_heat_cont.dot(b_heat_cont.T)

        # setup build matrices
        self.a_build = mmread(self.amtx_build).todense()
        self.e_build = mmread(self.emtx_build).todense()
        b_build = mmread(self.bmtx_build)
        c_build = mmread(self.cmtx_build)
        self.q_build = c_build.T.dot(c_build)
        self.g_build = b_build.dot(b_build.T)


        self.tol = 1e-7

    def template_nm_gmpare(self, instance, plus, linesearch, trans):
        """template to reduce code in the test cases"""
        if instance == "rail":
            a = self.a_rail
            e = self.e_rail
            q = self.q_rail
            g = self.g_rail
        elif instance == "iss":
            a = self.a_iss
            e = self.e_iss
            q = self.q_iss
            g = self.g_iss
        elif instance == "heat_cont":
            a = self.a_heat_cont
            e = self.e_heat_cont
            q = self.q_heat_cont
            g = self.g_heat_cont
        elif instance == "build":
            a = self.a_build
            e = self.e_build
            q = self.q_build
            g = self.g_build
        else:
            raise MESSErrorArguments("Instance not implemented.")

        x, _, _ = dense_nm_gmpare(None, a, e, q, g, plus=plus, linesearch=linesearch, trans=trans)
        nrmr, nrmrhs, relnrm = res2_gmpare(a, e, q, g, x, plus=plus, trans=trans)
        print("trans = {0:d}, plus={1:d}, res2 = {2:e}, nrmrhs = {3:e}, rel = res2/nrmrhs = {4:e}".format(trans,\
        plus, nrmr, nrmrhs, relnrm))
        self.assertLessEqual(relnrm, self.tol)


    def test_nm_gmpare_rail(self):
        """test"""
        print("\n")
        for linesearch in [False, True]:
            for trans in [MESS_OP_NONE, MESS_OP_TRANSPOSE]:
                for plus in [False, True]:
                    self.template_nm_gmpare("rail", plus, linesearch, trans)

    def test_nm_gmpare_iss(self):
        """test"""
        print("\n")
        for linesearch in [False, True]:
            for trans in [MESS_OP_NONE, MESS_OP_TRANSPOSE]:
                for plus in [False, True]:
                    self.template_nm_gmpare("iss", plus, linesearch, trans)

    def test_nm_gmpare_heat_cont(self):
        """test"""
        print("\n")
        for linesearch in [False, True]:
            for trans in [MESS_OP_NONE, MESS_OP_TRANSPOSE]:
                for plus in [False, True]:
                    self.template_nm_gmpare("heat_cont", plus, linesearch, trans)

    def test_nm_gmpare_build(self):
        """test"""
        print("\n")
        for linesearch in [False, True]:
            for trans in [MESS_OP_NONE, MESS_OP_TRANSPOSE]:
                for plus in [False, True]:
                    self.template_nm_gmpare("build", plus, linesearch, trans)


    def tearDown(self):
        """called after every test function"""
        pass

    def run(self, result=None):
        """ Stop after first error """
        if not result.errors:
            super(Testdensenmgmpare, self).run(result)
