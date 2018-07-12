% TEST_LRADI_TRIPLE_CHAIN_SO1 test for <a href="matlab:help mess_lradi">mess_lradi</a> using <a href="matlab:help mess_equation_glyap_so1">mess_equation_glyap_so1</a>.
%
%   TEST_LRADI_TRIPLE_CHAIN_SO1 tests the <a href="matlab:help mess_lradi">mess_lradi</a> solver.
%   TEST_LRADI_TRIPLE_CHAIN_SO1 creates a Lyapunov equation using <a href="matlab:help mess_equation_glyap_so1">mess_equation_glyap_so1</a> and solves it.
%
%   TEST_LRADI_TRIPLE_CHAIN_SO1 properties:
%       M           - real, sparse, n-by-n matrix.
%       D           - real, sparse, n-by-n matrix.
%       K           - real, sparse, n-by-n matrix.
%       lowerbound  - real, positive, scalar, lower bound for shift parameter criterion.
%       upperbound  - real, positive, scalar, upper bound for shift parameter criterion.
%       B           - real, dense, 2n-by-p matrix.
%       C           - real, dense,  p-by-2n matrix.
%       opt         - <a href="matlab:help mess_options">mess_options</a> instance.
%
%   See also MESS_LRADI, MESS_EQUATION_GLYAP_SO1.
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

classdef test_lradi_triple_chain_so1 < matlab.unittest.TestCase

    properties
        M
        D
        K
        lowerbound
        upperbound
        B
        C
        opt
    end

    methods(TestMethodSetup)
        function prepare_test(test)
            test.M =      mmread('@CMAKE_SOURCE_DIR@/tests/data/TripleChain/M_301.mtx');
            test.D =      mmread('@CMAKE_SOURCE_DIR@/tests/data/TripleChain/D_301.mtx');
            test.K =      mmread('@CMAKE_SOURCE_DIR@/tests/data/TripleChain/K_301.mtx');
            test.lowerbound = 1e-8;
            test.upperbound = 1e+8;

            % rest random number generator seed
            rng(1);

            test.B      = rand(2*length(test.M),1);
            test.C      = rand(1,2*length(test.M));
            test.opt    = mess_options();
            test.opt.nm.output = false;
            test.opt.adi.output = false;

        end
    end


    methods
        function template_so1(test)

            % create lyapunov equation
            if test.opt.type == mess_operation_t.MESS_OP_NONE
                % build equation
                eqn = mess_equation_glyap_so1(test.opt,test.M,test.D,test.K,test.B,test.lowerbound, test.upperbound);
            else
                eqn = mess_equation_glyap_so1(test.opt,test.M,test.D,test.K,test.C,test.lowerbound, test.upperbound);
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

        function so1_lyapunov(test)
            test.opt.type = mess_operation_t.MESS_OP_NONE;
            test.template_so1();
        end

        function so1_dual_lyapunov(test)
            test.opt.type = mess_operation_t.MESS_OP_TRANSPOSE;
            test.template_so1();
        end


    end
end


