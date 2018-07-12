% TEST_LYAP test for <a href="matlab:help mess_lyap">mess_lyap</a>,
%
%   TEST_LYAP tests the  <a href="matlab:help mess_lyap">mess_lyap</a>
%   solver. TEST_LYAP solves  standard and generalized Lyapunov eqauation
%   and computes the residual using <a href="matlab:help res2_lyap">res2_lyap</a>.
%
%   TEST_LYAP properties:
%       A           - sparse, real n-by-n matrix of Lyapunov equation.
%       B           - dense, real n-by-p matrix of Lyapunov equation.
%       C           - dense, real n-by-q matrix of Lyapunov equation.
%       E           - sparse, real, regular n-by-n matrix of Lyapunov equation.
%
%       Please note that the eigenvalues of A / (A,E) must lie in the left open
%       complex halfplane.
%
%   TEST_LYAP methods:
%       template_lyap               - solves Lyapunov equation and computes the residual.
%       generalized_lyapunov        - solves generalized Lyapunov equation.
%       standard_lyapunov           - solve standard Lyapunov equation.
%       generalized_dual_lyapunov   - solves generalized dual Lyapunov equation.
%       standard_dual_lyapunov      - solves standard dual Lyapunov equation.
%
%   See also MESS_LYAP.
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


classdef test_lyap < matlab.unittest.TestCase

    properties
        A = mmread('@CMAKE_SOURCE_DIR@/tests/data/Rail/A.mtx');
        B = full(mmread('@CMAKE_SOURCE_DIR@/tests/data/Rail/B.mtx'));
        C = full(mmread('@CMAKE_SOURCE_DIR@/tests/data/Rail/C.mtx'));
        E = mmread('@CMAKE_SOURCE_DIR@/tests/data/Rail/E.mtx');
    end

    methods
        function template_lyap(test,op,generalized)
            if op == mess_operation_t.MESS_OP_NONE
                if generalized
                    Z = mess_lyap(test.A,test.E,test.B);
                    res = res2_lyap(test.A,test.E,test.B,Z,op);
                else
                    Z = mess_lyap(test.A,[],test.B);
                    res = res2_lyap(test.A,[],test.B,Z,op);
                end
                rhs = norm(full(test.B),2)^2;
            else
                if generalized
                    Z = mess_lyap(test.A',test.E',test.C');
                    res = res2_lyap(test.A,test.E,test.C,Z,op);
                else
                    Z = mess_lyap(test.A',[],test.C');
                    res = res2_lyap(test.A,[],test.C,Z,op);
                end
                rhs = norm(full(test.C),2)^2;
            end

            verifyLessThanOrEqual(test, res/rhs, 1e-10);
        end
    end

    methods(Test)

        function generalized_lyapunov(test)
            test.template_lyap(mess_operation_t.MESS_OP_NONE,true);
        end

        function standard_lyapunov(test)
            test.template_lyap(mess_operation_t.MESS_OP_NONE,false);
        end

        function generalized_dual_lyapunov(test)
            test.template_lyap(mess_operation_t.MESS_OP_TRANSPOSE,true);
        end

        function standard_dual_lyapunov(test)
            test.template_lyap(mess_operation_t.MESS_OP_TRANSPOSE,false);
        end

    end

end
