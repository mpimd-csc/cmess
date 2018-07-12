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

""" Test conversion of matrix types."""
from __future__ import print_function, division
import unittest
from scipy import int8, int16, int32, int64, float32, float64, float128, complex64, complex128, complex256, matrix, \
array
import scipy.sparse as sparse
import pymess


class TestMatrix(unittest.TestCase):
    """Class for test conversion of matrix types."""
    input = matrix([[1.1, 2 + 3j, 2.1 + 5.12321312322j, 4, 5, 65], [10, 1, 2, 2, 3, 5]])
    result = input.H
    dtypes = [int8, int16, int32, int64, float32, float64, complex64, complex128]
    exdtypes = [float128, complex256]

    def setUp(self):
        """is called before every test function"""
        print("\n")

    def template_matrix(self, temp):
        """template to reduce code in the test cases"""
        #print ("Conversion from %s \t %s"%(type(temp),temp.dtype))
        output = pymess.test_matrix(temp)
        tmp = (self.result.astype(temp.dtype) - output) == 0
        self.assertTrue(tmp.all())

    #@unittest.skip("Skip is intended")
    def test_matrix_dense(self):
        """test"""
        for dt in self.dtypes:
            temp = self.input.astype(dt)
            self.template_matrix(temp)
        for dt in self.exdtypes:
            self.assertRaises(Exception, pymess.test_matrix, self.input.astype(dt))

    #@unittest.skip("Skip is intended")
    def test_matrix_ndarray_fortran(self):
        """test"""
        for dt in self.dtypes:
            temp = array(self.input, order='F', dtype=dt)
            self.template_matrix(temp)
        for dt in self.exdtypes:
            self.assertRaises(Exception, pymess.test_matrix, array(self.input, order='F', dtype=dt))

    #@unittest.skip("Skip is intended")
    def test_matrix_ndarray_c(self):
        """test"""
        for dt in self.dtypes:
            temp = array(self.input, order='C', dtype=dt)
            self.template_matrix(temp)
        for dt in self.exdtypes:
            self.assertRaises(Exception, pymess.test_matrix, array(self.input, order='C', dtype=dt))

    #@unittest.skip("Skip is intended")
    def test_matrix_csr(self):
        """test"""
        for dt in self.dtypes:
            temp = sparse.csr_matrix(self.input, dtype=dt)
            self.template_matrix(temp)
        for dt in self.exdtypes:
            self.assertRaises(Exception, pymess.test_matrix, sparse.csr_matrix(self.input, dtype=dt))

    #@unittest.skip("Skip is intended")
    def test_matrix_csc(self):
        """test"""
        for dt in self.dtypes:
            temp = sparse.csc_matrix(self.input, dtype=dt)
            self.template_matrix(temp)
        for dt in self.exdtypes:
            self.assertRaises(Exception, pymess.test_matrix, sparse.csc_matrix(self.input, dtype=dt))

    #@unittest.skip("Skip is intended")
    def test_matrix_coo(self):
        """test"""
        for dt in self.dtypes:
            temp = sparse.coo_matrix(self.input, dtype=dt)
            self.template_matrix(temp)
        for dt in self.exdtypes:
            self.assertRaises(Exception, pymess.test_matrix, sparse.coo_matrix(self.input, dtype=dt))

    def tearDown(self):
        """called after every test function"""
