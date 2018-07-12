% TEST_MESS_HAVE test various small functions of MEX-M.E.S.S.
%
%   TEST_MESS_HAVE tests the functions
%
%       <a href="matlab:help mess_version">mess_version</a>
%       <a href="matlab:help mess_version_verbose">mess_version_verbose</a>
%       <a href="matlab:help mess_version_major">mess_version_major</a>
%       <a href="matlab:help mess_version_minor">mess_version_minor</a>
%       <a href="matlab:help mess_version_patch">mess_version_patch</a>
%       <a href="matlab:help mess_git_id">mess_git_id</a>
%       <a href="matlab:help mess_git_branch">mess_git_branch</a>
%       <a href="matlab:help mess_call">mess_call</a> with different arguments
%
%   Author: Maximilian Behr (MPI Magdeburg, CSC)

%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, see <http://www.gnu.org/licenses/>.
%
% Copyright (C) Peter Benner, Martin Koehler, Jens Saak and others
%               2009-2018
%

classdef test_mess_have < matlab.unittest.TestCase

    methods(Test)
        function test_mess_version(test)
            mess_version();
            mess_version_verbose();
            major = mess_version_major();
            minor = mess_version_minor();
            patch = mess_version_patch();
            branch = mess_git_branch();
            gitid = mess_git_id();

            test.verifyClass(major,'double');
            test.verifyClass(minor,'double');
            test.verifyClass(patch,'double');
            test.verifyClass(branch,'char');
            test.verifyClass(gitid,'char');
            test.verifyTrue(0<=major);
            test.verifyTrue(0<=minor);
            test.verifyTrue(0<=patch);
        end

        function test_mess_have_funcs(test)
            test.verifyClass(mess_call('mess_is_debug'),        'logical');
            test.verifyClass(mess_call('mess_have_zlib'),       'logical');
            test.verifyClass(mess_call('mess_have_bzip2'),      'logical');
            test.verifyClass(mess_call('mess_have_umfpack'),    'logical');
            test.verifyClass(mess_call('mess_have_amd'),        'logical');
            test.verifyClass(mess_call('mess_have_colamd'),     'logical');
            test.verifyClass(mess_call('mess_have_cholmod'),    'logical');
            test.verifyClass(mess_call('mess_have_csparse'),    'logical');
            test.verifyClass(mess_call('mess_have_superlu'),    'logical');
            test.verifyClass(mess_call('mess_have_mklpardiso'), 'logical');
            test.verifyClass(mess_call('mess_have_arpack'),     'logical');
            test.verifyClass(mess_call('mess_have_matio'),      'logical');
            test.verifyClass(mess_call('mess_have_openmp'),     'logical');
            test.verifyClass(mess_call('mess_have_mess64'),     'logical');
        end

    end

end

