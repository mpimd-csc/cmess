% TEST_GSTEIN test for <a href="matlab:help mess_gstein">mess_gstein</a> using <a href="matlab:help mess_gstein">mess_gstein</a>.
%
%   TEST_GSTEIN tests the <a href="matlab:help mess_gstein">mess_gstein</a> solver.
%   TEST_GSTEIN solves a standard/generalized Stein equation using <a href="matlab:help mess_gstein">mess_gstein</a>
%   and computes the residual.
%   TEST_GSTEIN checks if the following dense standard / generalized Stein equations solved properly:
%
%   A  X A' -    X      = - Y      (1.1)
%   A' X A  -    X      = - Y      (1.2)
%   A  X A' - E  X E'   = - Y      (2.1)
%   A' X A  - E' X E    = - Y      (2.2)
%
%   TEST_GSTEIN properties:
%       ns          - vector, dimension of small dense matrix
%       general     - vector, with 0 or 1, 1 a generalized Lyapunov equation is solved
%       seed        - random seed
%       tol         - tolerance for absolute and relative residual
%       verbose     - logical
%
%   TEST_GSTEIN methods:
%       test_gstein_solve_only   - solves Stein equation and compute residual.
%       test_gstein_schur_solve  - solves Stein equation and solve again.
%
%   See also MESS_GSTEIN, RES2_STEIN.
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

classdef test_gstein < matlab.unittest.TestCase

    properties
        ns       = [1, 3, 5, 10, 40, 50, 100];
        seed     = 1;
        tol      = sqrt(eps);
        verbose  = false;
    end

    methods(Test)

        function test_gstein_solve_only(test)
            rng(test.seed);

            for n = test.ns
                for o = enumeration('mess_operation_t')'
                    for hasE = 0:1

                        %% create test data
                        A = rand(n,n);
                        if hasE
                            E = rand(n,n);
                        else
                            E = [];
                        end
                        Y1 = rand(n,n); Y1 = Y1 + Y1';
                        Y2 = rand(n,n); Y2 = Y2 + Y2';

                        %% solve only
                        X1 = mess_gstein(A, E, Y1);
                        X2 = mess_gstein(A, E, Y2, o);

                        %% compute residual
                        if hasE
                            res2_1 = res2_stein(A, E, Y1, X1, mess_operation_t.MESS_OP_NONE);
                            res2_2 = res2_stein(A, E, Y2, X2, o);
                        else
                            res2_1 = res2_stein(A, Y1, X1, mess_operation_t.MESS_OP_NONE);
                            res2_2 = res2_stein(A, Y2, X2, o);
                        end

                        rel2_1 = res2_1 / norm(Y1);
                        rel2_2 = res2_2 / norm(Y2);

                        if test.verbose
                            fprintf('n      = %d\n', n);
                            fprintf('hasE   = %d\n', hasE);
                            fprintf('tol    = %e\n', test.tol);
                            fprintf('res2_1 = %e\t rel2_1 = %e\n', res2_1, rel2_1);
                            fprintf('res2_2 = %e\t rel2_2 = %e\n', res2_2, rel2_2);
                        end

                        %% check residual
                        verifyLessThanOrEqual(test, res2_1, test.tol);  verifyLessThanOrEqual(test, rel2_1, test.tol);
                        verifyLessThanOrEqual(test, res2_2, test.tol);  verifyLessThanOrEqual(test, rel2_2, test.tol);
                    end
                end
            end
        end


        function test_gstein_schur_solve(test)
            rng(test.seed);

            for n = test.ns
                for o = enumeration('mess_operation_t')'
                    for hasE = 0:1

                        %% create test data
                        A = rand(n,n);
                        if hasE
                            E = rand(n,n);
                        else
                            E = [];
                        end
                        Y1 = rand(n,n); Y1 = Y1 + Y1';
                        Y2 = rand(n,n); Y2 = Y2 + Y2';
                        Y3 = rand(n,n); Y3 = Y3 + Y3';
                        Y4 = rand(n,n); Y4 = Y4 + Y4';

                        %% solve and get schur decomposition, solve again
                        if hasE
                            [X1, Ahat1, Ehat1, QA1, QE1]  = mess_gstein(A, E, Y1);
                            [X2, Ahat2, Ehat2, QA2, QE2]  = mess_gstein(A, E, Y2, o);
                            X3  = mess_gstein(A, E, Ahat1, Ehat1, QA1, QE1, Y3);
                            X4  = mess_gstein(A, E, Ahat2, Ehat2, QA2, QE2, Y4, o);
                        else
                            [X1, Ahat1, QA1] = mess_gstein(A, E, Y1);
                            [X2, Ahat2, QA2] = mess_gstein(A, E, Y2, o);
                            X3  = mess_gstein(A, Ahat1, QA1, Y3);
                            X4  = mess_gstein(A, Ahat2, QA2, Y4, o);
                        end

                        %% compute residual
                        if hasE
                            res2_1 = res2_stein(A, E, Y1, X1, mess_operation_t.MESS_OP_NONE);
                            res2_2 = res2_stein(A, E, Y2, X2, o);
                            res2_3 = res2_stein(A, E, Y3, X3, mess_operation_t.MESS_OP_NONE);
                            res2_4 = res2_stein(A, E, Y4, X4, o);
                        else
                            res2_1 = res2_stein(A, Y1, X1, mess_operation_t.MESS_OP_NONE);
                            res2_2 = res2_stein(A, Y2, X2, o);
                            res2_3 = res2_stein(A, Y3, X3, mess_operation_t.MESS_OP_NONE);
                            res2_4 = res2_stein(A, Y4, X4, o);
                        end

                        rel2_1 = res2_1 / norm(Y1);
                        rel2_2 = res2_2 / norm(Y2);
                        rel2_3 = res2_3 / norm(Y3);
                        rel2_4 = res2_4 / norm(Y4);

                        if test.verbose
                            fprintf('n      = %d\n', n);
                            fprintf('hasE   = %d\n', hasE);
                            fprintf('tol    = %e\n', test.tol);
                            fprintf('res2_1 = %e\t rel2_1 = %e\n', res2_1, rel2_1);
                            fprintf('res2_2 = %e\t rel2_2 = %e\n', res2_2, rel2_2);
                            fprintf('res2_3 = %e\t rel2_3 = %e\n', res2_3, rel2_3);
                            fprintf('res2_3 = %e\t rel2_4 = %e\n', res2_4, rel2_4);
                        end

                        %% check residual
                        verifyLessThanOrEqual(test, res2_1, test.tol);  verifyLessThanOrEqual(test, rel2_1, test.tol);
                        verifyLessThanOrEqual(test, res2_2, test.tol);  verifyLessThanOrEqual(test, rel2_2, test.tol);
                        verifyLessThanOrEqual(test, res2_3, test.tol);  verifyLessThanOrEqual(test, rel2_3, test.tol);
                        verifyLessThanOrEqual(test, res2_4, test.tol);  verifyLessThanOrEqual(test, rel2_4, test.tol);

                    end
                end
            end
        end

    end
end

