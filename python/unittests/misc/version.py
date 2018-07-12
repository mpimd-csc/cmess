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

"""Test of mess_version like functions."""
from __future__ import print_function, division
import unittest
import sys
import pymess

class Version(unittest.TestCase):
    """Class for test mess_version functions."""

    def setUp(self):
        """is called before every test function"""
        print("\n")

    #@unittest.skip("Skip is intended")
    def test_version(self):
        """test"""
        self.assertIsNone(pymess.mess_version())
        self.assertIsNone(pymess.mess_version_verbose())
        major = pymess.mess_version_major()
        minor = pymess.mess_version_minor()
        patch = pymess.mess_version_patch()
        gitid = pymess.mess_git_id()
        branch = pymess.mess_git_branch()
        self.assertTrue(major >= 0 and isinstance(major, int))
        self.assertTrue(minor >= 0 and isinstance(minor, int))
        self.assertTrue(patch >= 0 and isinstance(patch, int))
        if sys.version_info > (3, 0):
            self.assertTrue(isinstance(gitid, str))
            self.assertTrue(isinstance(branch, str))
        else:
            self.assertTrue(isinstance(gitid, unicode))
            self.assertTrue(isinstance(branch, unicode))


    def tearDown(self):
        """called after every test function"""
