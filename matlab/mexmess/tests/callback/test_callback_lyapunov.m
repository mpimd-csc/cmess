% TEST_CALLBACK_LYAPUNOV test the user defined callback equation from <a href="matlab:help myequation">myequation</a>.
%
%   TEST_CALLBACK_LYAPUNOV tests the user defined callback equation and solves
%   standard/generalized Lyapunov equations  <a href="matlab:help mess_lradi">mess_lradi</a>.
%
%   See also TEST_CALLBACK_RICCATI, MYEQUATION.
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

classdef test_callback_lyapunov < matlab.unittest.TestCase

    properties
        A = mmread('@CMAKE_SOURCE_DIR@/tests/data/Rail/A.mtx');
        B = full(mmread('@CMAKE_SOURCE_DIR@/tests/data/Rail/B.mtx'));
        C = full(mmread('@CMAKE_SOURCE_DIR@/tests/data/Rail/C.mtx'));
        E = mmread('@CMAKE_SOURCE_DIR@/tests/data/Rail/E.mtx');
        opt = mess_options();
    end

    methods
        function template_lradi(test,generalized)

            % create lyapunov equation
            if generalized
                if test.opt.type == mess_operation_t.MESS_OP_NONE
                    eqn = myequation(mess_equation_t.MESS_EQN_GLYAP,test.opt,test.A,test.E,test.B,[]);
                else
                    eqn = myequation(mess_equation_t.MESS_EQN_GLYAP,test.opt,test.A,test.E,test.C,[]);
                end
            else
                if test.opt.type == mess_operation_t.MESS_OP_NONE
                    eqn = myequation(mess_equation_t.MESS_EQN_LYAP, test.opt,test.A,[],test.B,[]);
                else
                    eqn = myequation(mess_equation_t.MESS_EQN_LYAP,test.opt,test.A,[],test.C,[]);
                end
            end

            % solve lyapunov equation
            [Z,status] = mess_lradi(eqn,test.opt);

            % compute residual
            if generalized
                if test.opt.type == mess_operation_t.MESS_OP_NONE
                    res2 = res2_lyap(test.A,test.E,test.B,Z,mess_operation_t.MESS_OP_NONE);
                    rhs = norm(test.B,2)^2;
                else
                    res2 = res2_lyap(test.A,test.E,test.C,Z,mess_operation_t.MESS_OP_TRANSPOSE);
                    rhs = norm(test.C,2)^2;
                end
            else
                if test.opt.type == mess_operation_t.MESS_OP_NONE
                    res2 = res2_lyap(test.A,[],test.B,Z,mess_operation_t.MESS_OP_NONE);
                    rhs = norm(test.B,2)^2;
                else
                    res2 = res2_lyap(test.A,[],test.C,Z,mess_operation_t.MESS_OP_TRANSPOSE);
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
            test.opt.adi.output = false;
            test.template_lradi(true);
        end

        function generalized_dual_lyapunov(test)
            test.opt.type = mess_operation_t.MESS_OP_TRANSPOSE;
            test.opt.adi.output = false;
            test.template_lradi(true);
        end

        function standard_lyapunov(test)
            test.opt.type = mess_operation_t.MESS_OP_NONE;
            test.opt.adi.output = false;
            test.template_lradi(false);
        end

        function standard_dual_lyapunov(test)
            test.opt.type = mess_operation_t.MESS_OP_TRANSPOSE;
            test.opt.adi.output = false;
            test.template_lradi(false);
        end

    end

end


