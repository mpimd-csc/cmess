% TEST_CALLBACK_RICCATI test the user defined callback equation from <a href="matlab:help myequation">myequation</a>.
%
%   TEST_CALLBACK_RICCATI tests the user defined callback equation and solves
%   standard/generalized Riccati equations  <a href="matlab:help mess_lrnm">mess_lrnm</a>.
%
%   See also TEST_CALLBACK_LYAPUNOV, MYEQUATION.
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

classdef test_callback_riccati < matlab.unittest.TestCase

    properties
        A = mmread('@CMAKE_SOURCE_DIR@/tests/data/Rail/A.mtx');
        B = full(mmread('@CMAKE_SOURCE_DIR@/tests/data/Rail/B.mtx'));
        C = full(mmread('@CMAKE_SOURCE_DIR@/tests/data/Rail/C.mtx'));
        E = mmread('@CMAKE_SOURCE_DIR@/tests/data/Rail/E.mtx');
        opt = mess_options();
    end

    methods
        function template_lrnm(test,generalized)

            % create lyapunov equation
            if generalized
                eqn = myequation(mess_equation_t.MESS_EQN_GRICCATI,test.opt,test.A,test.E,test.B,test.C);
            else
                eqn = myequation(mess_equation_t.MESS_EQN_RICCATI, test.opt,test.A,[],test.B,test.C);
            end

            % solve lyapunov equation
            [Z,status] = mess_lrnm(eqn,test.opt);

            % compute residual
            if generalized
                if test.opt.type == mess_operation_t.MESS_OP_NONE
                    res2 = res2_ric(test.A,test.E,test.B,test.C,Z,test.opt.type);
                    rhs = norm(full(test.B),2)^2;
                else
                    res2 = res2_ric(test.A,test.E,test.B,test.C,Z,test.opt.type);
                    rhs = norm(full(test.C),2)^2;
                end
            else
                if test.opt.type == mess_operation_t.MESS_OP_NONE
                    res2 = res2_ric(test.A,[],test.B,test.C,Z,test.opt.type);
                    rhs = norm(full(test.B),2)^2;
                else
                    res2 = res2_ric(test.A,[],test.B,test.C,Z,test.opt.type);
                    rhs = norm(full(test.C),2)^2;
                end
            end

            % print
            if test.opt.adi.output
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


            if test.opt.nm.output
                disp(status);
            end
        end

    end
    methods(Test)

        function generalized_dual_riccati(test)
            test.opt.type = mess_operation_t.MESS_OP_NONE;
            test.opt.nm.output = false;
            test.opt.adi.output = false;
            test.template_lrnm(true);
        end

        function generalized_riccati(test)
            test.opt.type = mess_operation_t.MESS_OP_TRANSPOSE;
            test.opt.nm.output = false;
            test.opt.adi.output = false;
            test.template_lrnm(true);
        end

        function standard_dual_riccati(test)
            test.opt.type = mess_operation_t.MESS_OP_NONE;
            test.opt.nm.output = false;
            test.opt.adi.output = false;
            test.template_lrnm(false);
        end

        function standard_riccati(test)
            test.opt.type = mess_operation_t.MESS_OP_TRANSPOSE;
            test.opt.nm.output = false;
            test.opt.adi.output = false;
            test.template_lrnm(false);
        end

    end

end


