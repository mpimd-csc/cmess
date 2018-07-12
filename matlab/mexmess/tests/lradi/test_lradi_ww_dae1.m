% TEST_LRADI_WW_DAE1 test for <a href="matlab:help mess_lradi">mess_lradi</a> using <a href="matlab:help mess_equation_glyap_dae1">mess_equation_glyap_dae1</a>.
%
%   TEST_LRADI_WW_DAE1 tests the <a href="matlab:help mess_lradi">mess_lradi</a> solver.
%   TEST_LRADI_WW_DAE1 creates a Lyapunov equation using <a href="matlab:help mess_equation_glyap_dae1">mess_equation_glyap_dae1</a> and solves it.
%
%   TEST_LRADI_WW_DAE1 properties:
%       MATmtx      - function handle, to get path of test matrices
%       E11         - matrix
%       A11         - matrix
%       A12         - matrix
%       A21         - matrix
%       A22         - matrix
%       B           - matrix
%       C           - matrix
%       opt         - <a href="matlab:help mess_options">mess_options</a> instance.
%
%   See also MESS_LRADI, MESS_EQUATION_GLYAP_DAE1.
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

classdef test_lradi_ww_dae1 < matlab.unittest.TestCase

    properties
        MATmtx
        E11
        A11
        A12
        A21
        A22
        B
        C
        opt
    end

    methods(TestMethodSetup)
        function prepare_test(test)

            test.MATmtx = @(name) sprintf('@CMAKE_SOURCE_DIR@/tests/data/bips98_606/%s.mtx',name);

            % read matrices
            test.E11 = mmread(test.MATmtx('E11'));
            test.A11 = mmread(test.MATmtx('A11'));
            test.A12 = mmread(test.MATmtx('A12'));
            test.A21 = mmread(test.MATmtx('A21'));
            test.A22 = mmread(test.MATmtx('A22'));
            test.B   = full(mmread(test.MATmtx('B')));
            test.C   = full(mmread(test.MATmtx('C')));

            % set random number generator seed
            rng(1);

            test.opt    = mess_options();
            test.opt.nm.output = false;
            test.opt.adi.output = false;

        end
    end

    methods
        function template_dae1(test)

            % create lyapunov equation
            if test.opt.type == mess_operation_t.MESS_OP_NONE
                eqn = mess_equation_glyap_dae1(test.opt,test.E11,test.A11, ...
                                                       test.A12,test.A21, ...
                                                       test.A22,test.B);
            else
                eqn = mess_equation_glyap_dae1(test.opt,test.E11,test.A11, ...
                                                       test.A12,test.A21, ...
                                                       test.A22,test.C);
            end

            % solve lyapunov equation
            [Z,status] = mess_lradi(eqn,test.opt);

            if test.opt.adi.output || test.opt.nm.output
                size(Z)
                fprintf('absolute/relative residual: %e\t %e\n',status.res2_norm,status.res2_norm/status.res2_0);
                disp(status);
            end
            verifyLessThanOrEqual(test,  status.res2_norm/status.res2_0, test.opt.adi.res2_tol);
        end
    end

    methods(Test)

        function dae1_lyapunov(test)
            test.opt.type = mess_operation_t.MESS_OP_NONE;
            test.template_dae1();
        end

        function dae1_dual_lyapunov(test)
            test.opt.type = mess_operation_t.MESS_OP_TRANSPOSE;
            test.template_dae1();
        end

    end
end
