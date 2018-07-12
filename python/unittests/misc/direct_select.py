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

"""Test of mess_have and direct_select like functions."""
from __future__ import print_function, division
import unittest
import pymess

class DirectSelect(unittest.TestCase):
    """Class for test mess_have and direct_select like functions."""

    def setUp(self):
        """is called before every test function"""
        print("\n")

    #@unittest.skip("Skip is intended")
    def test_mess_have(self):
        """test"""
        self.assertTrue(isinstance(pymess.mess_have_amd(), bool))
        self.assertTrue(isinstance(pymess.mess_have_arpack(), bool))
        self.assertTrue(isinstance(pymess.mess_have_bzip2(), bool))
        self.assertTrue(isinstance(pymess.mess_have_cholmod(), bool))
        self.assertTrue(isinstance(pymess.mess_have_colamd(), bool))
        self.assertTrue(isinstance(pymess.mess_have_csparse(), bool))
        self.assertTrue(isinstance(pymess.mess_have_matio(), bool))
        self.assertTrue(isinstance(pymess.mess_have_mess64(), bool))
        self.assertTrue(isinstance(pymess.mess_have_mklpardiso(), bool))
        self.assertTrue(isinstance(pymess.mess_have_openmp(), bool))
        self.assertTrue(isinstance(pymess.mess_have_superlu(), bool))
        self.assertTrue(isinstance(pymess.mess_have_umfpack(), bool))
        self.assertTrue(isinstance(pymess.mess_have_zlib(), bool))
        self.assertTrue(isinstance(pymess.mess_is_debug(), bool))

    def test_direct_select(self):
        """test"""

        pymess.direct_select(pymess.MESS_DIRECT_DEFAULT_LU)
        pymess.direct_select(pymess.MESS_DIRECT_SPARSE_LU)
        pymess.direct_select(pymess.MESS_DIRECT_LAPACK_LU)
        pymess.direct_select(pymess.MESS_DIRECT_BANDED_LU)

        if pymess.mess_have_umfpack():
            pymess.direct_select(pymess.MESS_DIRECT_UMFPACK_LU)
        else:
            self.assertRaises(Exception, pymess.direct_select, pymess.MESS_DIRECT_UMFPACK_LU)

        if pymess.mess_have_superlu():
            pymess.direct_select(pymess.MESS_DIRECT_SUPERLU_LU)
        else:
            self.assertRaises(Exception, pymess.direct_select, pymess.MESS_DIRECT_SUPERLU_LU)

        if pymess.mess_have_csparse():
            pymess.direct_select(pymess.MESS_DIRECT_CSPARSE_LU)
        else:
            self.assertRaises(Exception, pymess.direct_select, pymess.MESS_DIRECT_CSPARSE_LU)

        if pymess.mess_have_mklpardiso():
            pymess.direct_select(pymess.MESS_DIRECT_MKLPARDISO_LU)
        else:
            self.assertRaises(Exception, pymess.direct_select, pymess.MESS_DIRECT_MKLPARDISO_LU)

        # select the default solver after test
        pymess.direct_select(pymess.MESS_DIRECT_DEFAULT_LU)

    def test_direct_chol_select(self):
        """test"""

        pymess.direct_chol_select(pymess.MESS_DIRECT_DEFAULT_CHOLESKY)
        pymess.direct_chol_select(pymess.MESS_DIRECT_LAPACK_CHOLESKY)

        if pymess.mess_have_csparse():
            pymess.direct_chol_select(pymess.MESS_DIRECT_CSPARSE_CHOLESKY)
        else:
            self.assertRaises(Exception, pymess.direct_chol_select, pymess.MESS_DIRECT_CSPARSE_CHOLESKY)

        if pymess.mess_have_cholmod():
            pymess.direct_chol_select(pymess.MESS_DIRECT_CHOLMOD_CHOLESKY)
        else:
            self.assertRaises(Exception, pymess.direct_chol_select, pymess.MESS_DIRECT_CHOLMOD_CHOLESKY)

        # select the default solver after test
        pymess.direct_chol_select(pymess.MESS_DIRECT_DEFAULT_CHOLESKY)

    def test_multidirect_select(self):
        """test"""

        pymess.multidirect_select(pymess.MESS_MULTIDIRECT_SPARSE_LU)

        if pymess.mess_have_umfpack():
            pymess.multidirect_select(pymess.MESS_MULTIDIRECT_UMFPACK_LU)
        else:
            self.assertRaises(Exception, pymess.multidirect_select, pymess.MESS_MULTIDIRECT_UMFPACK_LU)

        # select the default solver after test
        pymess.multidirect_select(pymess.MESS_MULTIDIRECT_SPARSE_LU)

    def tearDown(self):
        """called after every test function"""
