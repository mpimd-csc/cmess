%RUN_TESTS  Runs all unit tests for MEX-M.E.S.S.
%
%   RUN_TESTS() runs all unit tests for MEX-M.E.S.S.
%
%   Author: Maximilian Behr (MPI Magdeburg, CSC)

%
%  This program is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation; either version 2 of the License, or
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, see <http://www.gnu.org/licenses/>.
%
%  Copyright (C) Peter Benner, Martin Koehler, Jens Saak and others
%                2009-2018
%

function ret = run_tests()
    %% import test package
    import matlab.unittest.TestSuite;

    %% create suites
    conversion  = TestSuite.fromFolder('@CMAKE_BINARY_DIR@/matlab/mexmess/tests/conversions/');
    misc        = TestSuite.fromFolder('@CMAKE_BINARY_DIR@/matlab/mexmess/tests/misc/');
    easy        = TestSuite.fromFolder('@CMAKE_BINARY_DIR@/matlab/mexmess/tests/easy/');
    lradi       = TestSuite.fromFolder('@CMAKE_BINARY_DIR@/matlab/mexmess/tests/lradi/');
    lrnm        = TestSuite.fromFolder('@CMAKE_BINARY_DIR@/matlab/mexmess/tests/lrnm/');
    callback    = TestSuite.fromFolder('@CMAKE_BINARY_DIR@/matlab/mexmess/tests/callback/');
    sylvester   = TestSuite.fromFolder('@CMAKE_BINARY_DIR@/matlab/mexmess/tests/sylvester/');
    glyap3      = TestSuite.fromFolder('@CMAKE_BINARY_DIR@/matlab/mexmess/tests/glyap3/');

    %% suite conversion
    suite = conversion;
    results = run(suite);
    disp(results);
    ret = any([results.Failed]);
    if ret, warning('conversion testsuite failed.'); return; end

    %% suite misc
    suite = misc;
    results = run(suite);
    disp(results);
    ret = any([results.Failed]);
    if ret, warning('misc testsuite failed.'); return; end

    %% suite easy
    suite = easy;
    results = run(suite);
    disp(results);
    ret = any([results.Failed]);
    if ret, warning('easy testsuite failed.'); return; end

    %% suite lradi
    suite = lradi;
    results = run(suite);
    disp(results);
    ret = any([results.Failed]);
    if ret, warning('lradi testsuite failed.'); return; end

    %% suite lrnm
    suite = lrnm;
    results = run(suite);
    disp(results);
    ret = any([results.Failed]);
    if ret, warning('lrnm testsuite failed.'); return; end

    %% suite callback
    suite = callback;
    results = run(suite);
    disp(results);
    ret = any([results.Failed]);
    if ret, warning('callback testsuite failed.'); return; end

    %% suite sylvester
    suite = sylvester;
    results = run(suite);
    disp(results);
    ret = any([results.Failed]);
    if ret, warning('sylvester testsuite failed.'); return; end

    %% suite glyap3
    suite = glyap3;
    results = run(suite);
    disp(results);
    ret = any([results.Failed]);
    if ret, warning('glyap3 testsuite failed.'); return; end

end
