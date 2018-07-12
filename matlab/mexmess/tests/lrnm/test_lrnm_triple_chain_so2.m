% TEST_LRNM_TRIPLE_CHAIN_SO2 test for <a href="matlab:help mess_lrnm">mess_lrnm</a> using <a href="matlab:help mess_equation_griccati_so2">mess_equation_griccati_so2</a>.
%
%   TEST_LRNM_TRIPLE_CHAIN_SO2 tests the <a href="matlab:help mess_lrnm">mess_lrnm</a> solver.
%   TEST_LRNM_TRIPLE_CHAIN_SO2 creates a Riccati equation using <a href="matlab:help mess_equation_griccati_so2">mess_equation_griccati_so2</a> and solves it.
%
%   TEST_LRNM_TRIPLE_CHAIN_SO2 properties:
%       M           - real, sparse, n-by-n matrix.
%       D           - real, sparse, n-by-n matrix.
%       K           - real, sparse, n-by-n matrix.
%       lowerbound  - real, positive, scalar, lower bound for shift parameter criterion.
%       upperbound  - real, positive, scalar, upper bound for shift parameter criterion.
%       B           - real, dense, 2n-by-p or p-by-n matrix.
%       C           - real, dense,  q-by-n or 2n-by-q matrix.
%       opt         - <a href="matlab:help mess_options">mess_options</a> instance.
%
%   See also MESS_LRNM, MESS_EQUATION_GRICCATI_SO2.
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

classdef test_lrnm_triple_chain_so2 < matlab.unittest.TestCase

    properties
        M
        D
        K
        lowerbound
        upperbound
        opt
        B
        C
    end

    methods(TestMethodSetup)
        function prepare_test(test)
            test.M = mmread('@CMAKE_SOURCE_DIR@/tests/data/TripleChain/M_301.mtx');
            test.D = mmread('@CMAKE_SOURCE_DIR@/tests/data/TripleChain/D_301.mtx');
            test.K = mmread('@CMAKE_SOURCE_DIR@/tests/data/TripleChain/K_301.mtx');
            test.lowerbound = 1e-8;
            test.upperbound = 1e+8;

            % rest random number generator seed
            rng(1);

            % options
            test.opt = mess_options();
            test.opt.adi.output = true;
            test.opt.nm.output = true;

        end
    end

    methods
        function template_so2(test)
            % build equation
            eqn = mess_equation_griccati_so2(test.opt,test.M,test.D,test.K,test.B,test.C,test.lowerbound, test.upperbound);

            % solve riccati equation
            [Z,status] = mess_lrnm(eqn,test.opt);

            if test.opt.adi.output || test.opt.nm.output
                size(Z)
                fprintf('absolute/relative residual: %e\t %e\n',status.res2_norm,status.res2_norm/status.res2_0);
                disp(status);
            end
            verifyLessThanOrEqual(test,  status.res2_norm/status.res2_0, test.opt.nm.res2_tol);
        end

    end

    methods(Test)

        function test_riccati_so2(test)
            test.opt.type = mess_operation_t.MESS_OP_NONE;
            test.B = rand(2*length(test.M),1);
            test.C = rand(1,length(test.M));
            test.template_so2();
        end

        function test_dual_riccati_so2(test)
            test.opt.type = mess_operation_t.MESS_OP_TRANSPOSE;
            test.B = rand(length(test.M),1);
            test.C = rand(1,2*length(test.M));
            test.template_so2();
        end

    end
end

