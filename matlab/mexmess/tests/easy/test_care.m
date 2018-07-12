% TEST_CARE test for <a href="matlab:help mess_care">mess_care</a>,
%
%   TEST_CARE tests the  <a href="matlab:help mess_care">mess_care</a> solver.
%   TEST_CARE solves  standard and generalized Riccati eqauation
%   and computes the residual using <a href="matlab:help res2_ric">res2_ric</a>.
%
%   TEST_CARE properties:
%       A           - sparse, real n-by-n matrix of Riccati equation.
%       B           - dense, real n-by-p matrix of Riccati equation.
%       C           - dense, real n-by-q matrix of Riccati equation.
%       E           - sparse, real, regular n-by-n matrix of Riccati equation.
%
%       Please note that the eigenvalues of A / (A,E) must lie in the left open
%       complex halfplane.
%
%   TEST_CARE methods:
%       template_care               - solves Riccati equation and computes the residual.
%       generalized_riccati         - solves generalized Riccati equation.
%       standard_riccati            - solve standard Riccati equation.
%       generalized_dual_riccati    - solves generalized dual Riccati equation.
%       standard_dual_riccati       - solves standard dual Riccati equation.
%
%   See also MESS_CARE.
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

classdef test_care < matlab.unittest.TestCase

    properties
        A = mmread('@CMAKE_SOURCE_DIR@/tests/data/Rail/A.mtx');
        B = full(mmread('@CMAKE_SOURCE_DIR@/tests/data/Rail/B.mtx'));
        C = full(mmread('@CMAKE_SOURCE_DIR@/tests/data/Rail/C.mtx'));
        E = mmread('@CMAKE_SOURCE_DIR@/tests/data/Rail/E.mtx');
    end

    methods
        function template_care(test,op,generalized)

            if op == mess_operation_t.MESS_OP_NONE
                if generalized
                    Z = mess_care(test.A',test.E',test.C',test.B');
                    res = res2_ric(test.A,test.E,test.B,test.C,Z,op);
                else
                    Z = mess_care(test.A',[],test.C',test.B');
                    res = res2_ric(test.A,[],test.B,test.C,Z,op);
                end
                rhs = norm(full(test.B),2)^2;
            else
                if generalized
                    Z = mess_care(test.A,test.E,test.B,test.C);
                    res = res2_ric(test.A,test.E,test.B,test.C,Z,op);
                else
                    Z = mess_care(test.A,[],test.B,test.C);
                    res = res2_ric(test.A,[],test.B,test.C,Z,op);
                end
                rhs = norm(full(test.C),2)^2;
            end

            verifyLessThanOrEqual(test, res/rhs, 1e-10);
        end
    end

    methods(Test)

        function generalized_riccati(test)
            test.template_care(mess_operation_t.MESS_OP_NONE,true);
        end

        function standard_riccati(test)
            test.template_care(mess_operation_t.MESS_OP_NONE,false);
        end

        function generalized_dual_riccati(test)
            test.template_care(mess_operation_t.MESS_OP_TRANSPOSE,true);
        end

        function standard_dual_riccati(test)
            test.template_care(mess_operation_t.MESS_OP_TRANSPOSE,false);
        end

    end

end
