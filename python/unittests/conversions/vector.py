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

"""Test conversion of vector types."""
from __future__ import print_function, division
import unittest
from scipy import int8, int16, int32, int64, float32, float64, float128, complex64, complex128, complex256, array
from scipy import all as spall
import pymess


class TestVector(unittest.TestCase):
    """Class for test conversion of matrix types."""
    input = array([1.1, 2 + 3j, 2.1 + 5.12321312322j, 4, 5, 65])
    dtypes = [int8, int16, int32, int64, float32, float64, complex64, complex128]
    exdtypes = [float128, complex256]

    def setUp(self):
        """is called before every test function"""
        print("\n")

    def template_vector(self, temp):
        """template to reduce code in the test cases"""
        output = pymess.test_vector(temp)
        self.assertTrue(spall((output - temp) == 0))

    #@unittest.skip("Skip is intended")
    def test_vector(self):
        """test"""
        for dt in self.dtypes:
            temp = self.input.astype(dt)
            self.template_vector(temp)
        for dt in self.exdtypes:
            self.assertRaises(Exception, pymess.test_vector, self.input.astype(dt))

    def tearDown(self):
        """called after every test function"""
