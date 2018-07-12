% TEST_LRADI_FILTER test for <a href="matlab:help mess_lradi">mess_lradi</a> using <a href="matlab:help mess_equation_glyap">mess_equation_glyap</a>.
%
%   TEST_LRADI_FILTER tests the <a href="matlab:help mess_lradi">mess_lradi</a> solver.
%   TEST_LRADI_FILTER creates a standard/generalized Lyapunov equation using <a href="matlab:help mess_equation_glyap">mess_equation_glyap</a>
%   solves the equation and computes the residual using  <a href="matlab:help res2_lyap">res2_lyap</a>.
%   The residual is compared against the residual obtained from
%   <a href="matlab:help mess_lradi">mess_lradi</a> contained in the <a href="matlab:help mess_status">mess_status</a> instance.
%
%   TEST_LRADI_FILTER properties:
%       A           - sparse, real n-by-n matrix of Lyapunov equation.
%       B           - dense, real n-by-p matrix of Lyapunov equation.
%       C           - dense, real n-by-q matrix of Lyapunov equation.
%       E           - sparse, real, regular n-by-n matrix of Lyapunov equation.
%       opt         - <a href="matlab:help mess_options">mess_options</a> instance.
%
%       Please not that the eigenvalues of A / (A,E) must lie in the left open
%       complex halfplane.
%
%   TEST_LRADI_FILTER methods:
%       template_lradi              - solves Lyapunov equation and computes the residual.
%       generalized_lyapunov        - solves generalized Lyapunov equation.
%       standard_lyapunov           - solve standard Lyapunov equation.
%       generalized_dual_lyapunov   - solves generalized dual Lyapunov equation.
%       standard_dual_lyapunov      - solves standard dual Lyapunov equation.
%
%   See also MESS_LRADI, MESS_EQUATION_GLYAP.
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

classdef test_lradi_filter < matlab.unittest.TestCase

    properties
        A =      mmread('@CMAKE_SOURCE_DIR@/tests/data/filter2D/filter2D.A');
        B = full(mmread('@CMAKE_SOURCE_DIR@/tests/data/filter2D/filter2D.B'));
        C = full(mmread('@CMAKE_SOURCE_DIR@/tests/data/filter2D/filter2D.C'));
        E =      mmread('@CMAKE_SOURCE_DIR@/tests/data/filter2D/filter2D.E');
        opt = mess_options();
    end

    methods
        function template_lradi(test,generalized)
            % create lyapunov equation
            if generalized
                if test.opt.type == mess_operation_t.MESS_OP_NONE
                    eqn = mess_equation_glyap(test.opt,test.A,test.E,test.B);
                else
                    eqn = mess_equation_glyap(test.opt,test.A,test.E,test.C);
                end
            else
                if test.opt.type == mess_operation_t.MESS_OP_NONE
                    eqn = mess_equation_glyap(test.opt,test.A,[],test.B);
                else
                    eqn = mess_equation_glyap(test.opt,test.A,[],test.C);
                end
            end

            % solve lyapunov equation
            [Z,status] = mess_lradi(eqn,test.opt);

            % compute residual
            if generalized
                if test.opt.type == mess_operation_t.MESS_OP_NONE
                    res2 = res2_lyap(test.A,test.E,test.B,Z,test.opt.type);
                    rhs = norm(test.B,2)^2;
                else
                    res2 = res2_lyap(test.A,test.E,test.C,Z,test.opt.type);
                    rhs = norm(test.C,2)^2;
                end
            else
                if test.opt.type == mess_operation_t.MESS_OP_NONE
                    res2 = res2_lyap(test.A,[],test.B,Z,test.opt.type);
                    rhs = norm(test.B,2)^2;
                else
                    res2 = res2_lyap(test.A,[],test.C,Z,test.opt.type);
                    rhs = norm(test.C,2)^2;
                end
            end

            % print
            if test.opt.adi.output
                fprintf('absolute residual: MEX-M.E.S.S.=%e\t C-M.E.S.S.=%e\n',res2, status.res2_norm);
                fprintf('relative residual: MEX-M.E.S.S.=%e\t C-M.E.S.S.=%e\n',res2/rhs, status.res2_norm/status.res2_0);
            end

            % test residual
            verifyLessThanOrEqual(test,  status.res2_norm/status.res2_0, test.opt.adi.res2_tol);

            % compare residual from res2_ric
            verifyLessThanOrEqual(test,  abs(status.res2_norm/status.res2_0 - res2/rhs), sqrt(eps));
            verifyLessThanOrEqual(test,  abs(status.res2_norm/status.res2_0 - res2/rhs)/(res2/rhs), 1e-2);
            verifyLessThanOrEqual(test,  abs(status.res2_0 - rhs), sqrt(eps));
            verifyLessThanOrEqual(test,  abs(status.res2_0 - rhs)/rhs, 1e-2);

            if test.opt.adi.output
                disp(status);
            end
        end

    end
    methods(Test)

        function generalized_lyapunov(test)
            test.opt.type = mess_operation_t.MESS_OP_NONE;
            test.opt.adi.output = true;
            test.template_lradi(true);
        end

        function generalized_dual_lyapunov(test)
            test.opt.type = mess_operation_t.MESS_OP_TRANSPOSE;
            test.opt.adi.output = true;
            test.template_lradi(true);
        end

        function standard_lyapunov(test)
            test.opt.type = mess_operation_t.MESS_OP_NONE;
            test.opt.adi.output = true;
            test.template_lradi(false);
        end

        function standard_dual_lyapunov(test)
            test.opt.type = mess_operation_t.MESS_OP_TRANSPOSE;
            test.opt.adi.output = true;
            test.template_lradi(false);
        end

    end

end

