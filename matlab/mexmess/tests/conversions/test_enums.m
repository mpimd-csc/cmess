% TEST_ENUMS test the enumeration type conversion between MEX-M.E.S.S. and C-M.E.S.S.
%
%   TEST_ENUMS calls for each enumeration type T a internal routine
%   mex_test_T which converts the enum from MEX-M.E.S.S. to C-M.E.S.S. and
%   back. The converted enum is returned and should be the same as the input.
%
%   TEST_ENUMS properties:
%       output          - logical, turn output on/off.
%
%   TEST_ENUMS methods:
%       test_mess_enum_t                - conversion of various enumeration types
%       mess_equation_t                 - conversion of <a href=matlab:help mess_equation_t>mess_equation_t</a>
%       mess_memusage_t                 - conversion of <a href=matlab:help mess_memusage_t>mess_memusage_t</a>
%       mess_operation_t                - conversion of <a href=matlab:help mess_operation_t>mess_operation_t</a>
%       mess_parameter_t                - conversion of <a href=matlab:help mess_parameter_t>mess_parameter_t</a>
%       mess_residual_t                 - conversion of <a href=matlab:help mess_residual_t>mess_residual_t</a>
%       mess_direct_cholpackage_t       - conversion of <a href=matlab:help mess_direct_cholpackage_t>mess_direct_cholpackage_t</a>
%       mess_direct_lupackage_t         - conversion of <a href=matlab:help mess_direct_lupackage_t>mess_direct_lupackage_t</a>
%       mess_multidirect_t              - conversion of <a href=matlab:help mess_multidirect_t>mess_multidirect_t</a>
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

classdef test_enums < matlab.unittest.TestCase

    properties
       output   % logical, turn output on/off.
    end

    methods
        function test_mess_enum_t(test,test_func,enum1)
            enum2 = mess_call(test_func,enum1);
            if test.output
                fprintf('C:%s - MATLAB:%s\n',char(enum2),char(enum1));
            end
            verifyTrue(test,enum1==enum2);
        end
    end

    methods(TestMethodSetup)
        function prepare_test(test)
            test.output = false;
        end
    end


    methods(Test)

        function mess_equation_t(test)
            for enum = enumeration('mess_equation_t')'
                test.test_mess_enum_t('mex_test_mess_equation_t', enum);
            end
        end

        function mess_memusage_t(test)
            for enum = enumeration('mess_memusage_t')'
                test.test_mess_enum_t('mex_test_mess_memusage_t', enum);
            end
        end

        function mess_operation_t(test)
            for enum = enumeration('mess_operation_t')'
                test.test_mess_enum_t('mex_test_mess_operation_t', enum);
            end
        end

        function mess_parameter_t(test)
            for enum = enumeration('mess_parameter_t')'
                test.test_mess_enum_t('mex_test_mess_parameter_t', enum);
            end
        end

        function mess_residual_t(test)
            for enum = enumeration('mess_residual_t')'
                test.test_mess_enum_t('mex_test_mess_residual_t', enum);
            end
        end

        function mess_direct_cholpackage_t(test)
            for enum = enumeration('mess_direct_cholpackage_t')'
                test.test_mess_enum_t('mex_test_mess_direct_cholpackage_t', enum);
            end
        end

        function mess_direct_lupackage_t(test)
            for enum = enumeration('mess_direct_lupackage_t')'
                test.test_mess_enum_t('mex_test_mess_direct_lupackage_t', enum);
            end
        end

        function mess_multidirect_t(test)
            for enum = enumeration('mess_multidirect_t')'
                test.test_mess_enum_t('mex_test_mess_multidirect_t', enum);
            end
        end

        function mess_norm_t(test)
            for enum = enumeration('mess_norm_t')'
                test.test_mess_enum_t('mex_test_mess_norm_t', enum);
            end
        end

    end

end



