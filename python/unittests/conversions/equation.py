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

""" Tests for :class:`Equation`."""
from __future__ import print_function, division
import unittest
import pymess



class TestEquation(unittest.TestCase):
    """Test if instance of subclass of pymess is correct found"""

    def setUp(self):
        """called before every test function"""
        print("\n")

    def test_class(self):
        """test some expected functionality of :class:`~pymess.equation.Equation`."""
        # test 1
        cl1 = pymess.Equation("test")
        ret = pymess.test_equation(cl1)
        self.assertEqual("test - 1 - 2", ret)

        # test 2
        class Eqn(object):
            """Only for test."""
            def __init__(self, name):
                self.name = name

        cl2 = Eqn("TEST")
        ret = pymess.test_equation(cl2)
        self.assertEqual(None, ret)

        # test 3
        class Eqn2(pymess.Equation):
            """Only for test."""
            def __init__(self, name):
                self.name = name

        cl3 = Eqn2("TEsts2")
        ret = pymess.test_equation(cl3)
        self.assertEqual("TEsts2 - 1 - 2", ret)

        # test 4
        class Eqn3(pymess.Equation):
            """Only for test."""
            def __init__(self, name):
                self.name = name

        ret = pymess.test_equation(Eqn3)
        self.assertEqual(None, ret)

    def tearDown(self):
        """called after every test function"""
