% TEST_LRNM_FILTER test for <a href="matlab:help mess_lrnm">mess_lrnm</a> using <a href="matlab:help mess_equation_griccati">mess_equation_griccati</a>.
%
%   TEST_LRNM_FILTER tests the <a href="matlab:help mess_lrnm">mess_lrnm</a> solver.
%   TEST_LRNM_FILTER creates a standard/generalized Riccati equation using <a href="matlab:help mess_equation_griccati">mess_equation_griccati</a>
%   solves the equation and computes the residual using  <a href="matlab:help res2_ric">res2_ric</a>.
%   The residual is compared against the residual obtained from <a href="matlab:help mess_lrnm">mess_lrnm</a> contained in the <a href="matlab:help mess_status">mess_status</a> instance.
%
%   TEST_LRNM_FILTER properties:
%       A           - sparse, real n-by-n matrix of Riccati equation.
%       B           - dense, real n-by-p matrix of Riccati equation.
%       C           - dense, real n-by-q matrix of Riccati equation.
%       E           - sparse, real, regular n-by-n matrix of Riccati equation.
%       opt         - <a href="matlab:help mess_options">mess_options</a> instance.
%
%       Please not that the eigenvalues of A / (A,E) must lie in the left open
%       complex halfplane.
%
%   TEST_LRNM_FILTER methods:
%       template_lrnm                       - solves Riccati equation and computes the residual.
%       filter_generalized_riccati          - solves generalized Riccati equation.
%       filter_standard_riccati             - solves standard Riccati equation.
%       filter_generalized_dual_riccati     - solves generalized dual Riccati equation.
%       filter_standard_dual_riccati        - solves standard dual Riccati equation.
%
%   See also MESS_LRNM, MESS_EQUATION_GRICCATI.
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

classdef test_lrnm_filter < matlab.unittest.TestCase

    properties
        A   = mmread('@CMAKE_SOURCE_DIR@/tests/data/filter2D/filter2D.A');
        B   = full(mmread('@CMAKE_SOURCE_DIR@/tests/data/filter2D/filter2D.B'));
        C   = full(mmread('@CMAKE_SOURCE_DIR@/tests/data/filter2D/filter2D.C'));
        E   =      mmread('@CMAKE_SOURCE_DIR@/tests/data/filter2D/filter2D.E');
        opt = mess_options();
    end

    methods
        function template_lrnm(test, generalized)
            test.opt.nm.output = true;
            test.opt.adi.output = true;

            % create equation
            if test.opt.type == mess_operation_t.MESS_OP_NONE
                if generalized
                    eqn = mess_equation_griccati(test.opt,test.A,test.E,test.B,test.C);
                else
                    eqn = mess_equation_griccati(test.opt,test.A,[],test.B,test.C);
                end
            else
                if generalized
                    eqn = mess_equation_griccati(test.opt,test.A,test.E,test.B,test.C);
                else
                    eqn = mess_equation_griccati(test.opt,test.A,[],test.B,test.C);
                end
            end

            % solve riccait equation
            [Z,status] = mess_lrnm(eqn,test.opt);

            % compute residual
            if test.opt.type == mess_operation_t.MESS_OP_NONE
                if generalized
                    res2 = res2_ric(test.A,test.E,test.B,test.C,Z,test.opt.type);
                else
                    res2 = res2_ric(test.A,[],test.B,test.C,Z,test.opt.type);
                end
                rhs = norm(full(test.B),2)^2;
            else
                if generalized
                    res2 = res2_ric(test.A,test.E,test.B,test.C,Z,test.opt.type);
                else
                    res2 = res2_ric(test.A,[],test.B,test.C,Z,test.opt.type);
                end
                rhs = norm(full(test.C),2)^2;
            end

            % output
            if test.opt.adi.output || test.opt.nm.output
                fprintf('absolute residual: MEX-M.E.S.S.=%e\t C-M.E.S.S.=%e\n',res2, status.res2_norm);
                fprintf('relative residual: MEX-M.E.S.S.=%e\t C-M.E.S.S.=%e\n',res2/rhs, status.res2_norm/status.res2_0);
            end

            % test residual
            verifyLessThanOrEqual(test,  status.res2_norm/status.res2_0, test.opt.nm.res2_tol);

            % compare residual from res2_ric
            verifyLessThanOrEqual(test,  abs(status.res2_norm/status.res2_0 - res2/rhs), sqrt(eps));
            verifyLessThanOrEqual(test,  abs(status.res2_norm/status.res2_0 - res2/rhs)/(res2/rhs), 1e-2);
            verifyLessThanOrEqual(test,  abs(status.res2_0 - rhs), sqrt(eps));
            verifyLessThanOrEqual(test,  abs(status.res2_0 - rhs)/rhs, 1e-2);

        end
    end


    methods(Test)

        function filter_generalized_riccati(test)
           test.opt.type = mess_operation_t.MESS_OP_TRANSPOSE;
           test.template_lrnm(true);
        end

        function filter_standard_riccati(test)
           test.opt.type = mess_operation_t.MESS_OP_TRANSPOSE;
           test.template_lrnm(false);
        end

        function filter_generalized_dual_riccati(test)
           test.opt.type = mess_operation_t.MESS_OP_NONE;
           test.template_lrnm(true);
        end

        function filter_standard_dual_riccati(test)
           test.opt.type = mess_operation_t.MESS_OP_NONE;
           test.template_lrnm(false);
        end

    end

end



