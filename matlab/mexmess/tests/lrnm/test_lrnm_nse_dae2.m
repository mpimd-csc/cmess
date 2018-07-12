% TEST_LRNM_NSE_DAE2 test for <a href="matlab:help mess_lrnm">mess_lrnm</a> using <a href="matlab:help mess_equation_griccati_dae2">mess_equation_griccati_dae2</a>.
%
%   TEST_LRNM_NSE_DAE2 tests the <a href="matlab:help mess_lrnm">mess_lrnm</a> solver.
%   TEST_LRNM_NSE_DAE2 creates a Riccati equation using <a href="matlab:help mess_equation_griccati_dae2">mess_equation_griccati_dae2</a> and solves it.
%
%   TEST_LRNM_NSE_DAE2 properties:
%       RE          - scalar, real positive Reynolds Number
%       delta       - sclar, negativ delta, standard value: -0.02.
%       MATmtx      - function handle, to get path of test matrices
%       M           - real, sparse, symmetric positive definite nv-by-nv matrix.
%       A           - real, sparse, nv-by-nv matrix.
%       G           - real, sparse, full rank, nv-by-np matrix.
%       B           - real, dense,  nv-by-p matrix.
%       C           - real, dense,  q-by-nv matrix.
%       K0          - stabilizing feedback.
%       K1          - stabilizing feedback.
%       opt         - <a href="matlab:help mess_options">mess_options</a> instance.
%
%
%   See also MESS_LRNM, MESS_EQUATION_GRICCATI_DAE2.
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

classdef test_lrnm_nse_dae2 < matlab.unittest.TestCase

   properties
        RE
        delta
        MATmtx
        M
        A
        G
        B
        C
        K0
        K1
        opt
    end

    methods
        function template_lrnm_dae2(test)

            % build equation
            eqn = mess_equation_griccati_dae2(test.opt,test.M,test.A,test.G,test.B,test.C,test.delta);

            if test.opt.type == mess_operation_t.MESS_OP_NONE
                test.opt.nm.K0 = test.K1;
            else
                test.opt.nm.K0 = test.K0;
            end


            % solve riccati equation
            [Z,status] = mess_lrnm(eqn,test.opt);

            % check resdiual
            verifyLessThanOrEqual(test,  status.res2_norm/status.res2_0, test.opt.nm.res2_tol);

            if test.opt.nm.output || test.opt.adi.output
                size(Z)
                fprintf('absolute/relative residual: %e\t %e\n',status.res2_norm,status.res2_norm/status.res2_0);
                disp(status);
            end

        end
    end

    methods(TestMethodSetup)

        function prepare_test(test)
            % setup reynoldsnumber and delta
            test.RE = 300;
            test.delta = -0.02;
            test.MATmtx = @(RE,name) sprintf('@CMAKE_SOURCE_DIR@/tests/data/NSE/NSE_RE_%d_lvl1_%s.mtx',RE,name);

            % read matrices
            test.M = mmread(test.MATmtx(test.RE,'M'));
            test.A = mmread(test.MATmtx(test.RE,'A'));
            test.G = mmread(test.MATmtx(test.RE,'G'));
            test.B = full(mmread(test.MATmtx(test.RE,'B')));
            test.C = full(mmread(test.MATmtx(test.RE,'C')));

            % get feedback matrix
            if test.RE >= 300
                test.K0 = mmread(test.MATmtx(test.RE,'Feed0'));
                test.K1 = mmread(test.MATmtx(test.RE,'Feed1'));
            else
                test.K0 = [];
                test.K1 = [];
            end

            % options structure
            test.opt = mess_options();
            test.opt.nm.output = true;
            test.opt.adi.output = false;
            test.opt.nm.res2_tol = 1e-3;
        end
    end



    methods(Test)

        function test_riccati_dae2(test)
            test.opt.type = mess_operation_t.MESS_OP_TRANSPOSE;
            test.template_lrnm_dae2();
        end

        function test_dual_riccati_dae2(test)
            test.opt.type = mess_operation_t.MESS_OP_NONE;
            test.template_lrnm_dae2();
        end

    end

end
