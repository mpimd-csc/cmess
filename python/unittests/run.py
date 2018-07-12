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

"""Run all tests in testsuite."""
from __future__ import print_function, division

import unittest

import sys

import conversions
import misc
import easy
import callback
import lradi
import lrnm
import glyap3


def main():
    """Run all tests in testsuite."""
    #creating test suites and add test cases
    conversionssuite = unittest.TestSuite()
    conversionssuite.addTests(unittest.makeSuite(conversions.TestMatrix))
    conversionssuite.addTests(unittest.makeSuite(conversions.TestMatrix2))
    conversionssuite.addTests(unittest.makeSuite(conversions.TestVector))
    conversionssuite.addTests(unittest.makeSuite(conversions.TestEquation))
    #unittest.TextTestRunner(verbosity=2).run(conversionssuite)

    miscsuite = unittest.TestSuite()
    miscsuite.addTests(unittest.makeSuite(misc.Version))
    miscsuite.addTests(unittest.makeSuite(misc.DirectSelect))
    #unittest.TextTestRunner(verbosity=2).run(miscsuite)

    easysuite = unittest.TestSuite()
    easysuite.addTests(unittest.makeSuite(easy.Testlyap))
    easysuite.addTests(unittest.makeSuite(easy.Testsylvestersparsedense))
    easysuite.addTests(unittest.makeSuite(easy.Testdensenmgmpare))
    easysuite.addTests(unittest.makeSuite(easy.Testcare))
    #unittest.TextTestRunner(verbosity=2).run(easysuite)

    callbacksuite = unittest.TestSuite()
    callbacksuite.addTests(unittest.makeSuite(callback.Lyapunov))
    callbacksuite.addTests(unittest.makeSuite(callback.Riccati))
    #unittest.TextTestRunner(verbosity=2).run(callbacksuite)

    lradisuite = unittest.TestSuite()
    lradisuite.addTests(unittest.makeSuite(lradi.Filter))
    lradisuite.addTests(unittest.makeSuite(lradi.TripleChainSO1))
    lradisuite.addTests(unittest.makeSuite(lradi.TripleChainSO2))
    lradisuite.addTests(unittest.makeSuite(lradi.WWDAE1))
    lradisuite.addTests(unittest.makeSuite(lradi.NSEDAE2))
    #unittest.TextTestRunner(verbosity=2).run(lradisuite)

    lrnmsuite = unittest.TestSuite()
    lrnmsuite.addTests(unittest.makeSuite(lrnm.Filter))
    lrnmsuite.addTests(unittest.makeSuite(lrnm.TripleChainSO1))
    lrnmsuite.addTests(unittest.makeSuite(lrnm.TripleChainSO2))
    lrnmsuite.addTests(unittest.makeSuite(lrnm.WWDAE1))
    lrnmsuite.addTests(unittest.makeSuite(lrnm.NSEDAE2))
    #unittest.TextTestRunner(verbosity=2).run(lrnmsuite)

    glyap3suite = unittest.TestSuite()
    glyap3suite.addTests(unittest.makeSuite(glyap3.TestGlyap3))
    #unittest.TextTestRunner(verbosity=2).run(glyap3suite)

    #set unittest framework
    if len(sys.argv) == 1:
        print("RUN ALL TESTS")
        allsuites = unittest.TestSuite([conversionssuite, miscsuite, easysuite, \
                                        callbacksuite, lradisuite, lrnmsuite,   \
                                        glyap3suite])
    elif len(sys.argv) == 2:
        if sys.argv[1].upper() == "EASY":
            allsuites = unittest.TestSuite([easysuite])
        elif sys.argv[1].upper() == "CONVERSION":
            allsuites = unittest.TestSuite([conversionssuite])
        elif sys.argv[1].upper() == "MISC":
            allsuites = unittest.TestSuite([miscsuite])
        elif sys.argv[1].upper() == "LRADI":
            allsuites = unittest.TestSuite([lradisuite])
        elif sys.argv[1].upper() == "CALLBACK":
            allsuites = unittest.TestSuite([callbacksuite])
        elif sys.argv[1].upper() == "LRNM":
            allsuites = unittest.TestSuite([lrnmsuite])
        elif sys.argv[1].upper() == "GLYAP3":
            allsuites = unittest.TestSuite([glyap3suite])
        else:
            raise ValueError("Wrong TestSuite. Allowed names are: EASY, CONVERSION, LRADI, CALLBACK, LRNM, GLYAP3.")

    else:
        raise ValueError("Wrong number of input arguments.")

    #run unittest
    ret = not unittest.TextTestRunner(verbosity=2).run(allsuites).wasSuccessful()
    sys.exit(ret)

if __name__ == "__main__":
    main()
