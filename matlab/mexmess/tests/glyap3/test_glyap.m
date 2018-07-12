% TEST_GLYAP test for <a href="matlab:help mess_glyap">mess_glyap</a> using <a href="matlab:help mess_glyap">mess_glyap</a>.
%
%   TEST_GLYAP tests the <a href="matlab:help mess_glyap">mess_glyap</a> solver.
%   TEST_GLYAP solves a standard/generalized Lyapunov equation using <a href="matlab:help mess_glyap">mess_glyap</a>
%   and computes the residual.
%   TEST_GLYAP checks if the following dense standard / generalized Lyapunov equations solved properly:
%
%   A  X    +    X A'  = - Y      (1.1)
%   A' X    +    X A   = - Y      (1.2)
%   A  X E' + E  X A'  = - Y      (2.1)
%   A' X E  + E' X A   = - Y      (2.2)
%
%   TEST_GLYAP properties:
%       ns          - vector, dimension of small dense matrix
%       general     - vector, with 0 or 1, 1 a generalized Lyapunov equation is solved
%       seed        - random seed
%       tol         - tolerance for absolute and relative residual
%       verbose     - logical
%
%   TEST_GLYAP methods:
%       res2_glyap              - computes the residual of a standard / generalized Lyapunov equation.
%       test_glyap_solve_only   - solves Lyapunov equation and compute residual.
%       test_glyap_schur_solve  - solves Lyapunov equation  and solve again.
%
%   See also MESS_GLYAP, RES2_LYAP.
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

classdef test_glyap < matlab.unittest.TestCase

    properties
        ns       = [1, 3, 5, 10, 40, 50, 100];
        seed     = 1;
        tol      = sqrt(eps);
        verbose  = false;
    end

    methods(Test)

        function test_glyap_solve_only(test)
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
                        X1 = mess_glyap(A, E, Y1);
                        X2 = mess_glyap(A, E, Y2, o);

                        %% compute residual
                        res2_1 = res2_glyap(A, E, Y1, X1, mess_operation_t.MESS_OP_NONE);
                        rel2_1 = res2_1 / norm(Y1);

                        res2_2 = res2_glyap(A, E, Y2, X2, o);
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


        function test_glyap_schur_solve(test)
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
                            [X1, Ahat1, Ehat1, QA1, QE1]  = mess_glyap(A, E, Y1);
                            [X2, Ahat2, Ehat2, QA2, QE2]  = mess_glyap(A, E, Y2, o);
                            X3  = mess_glyap(A, E, Ahat1, Ehat1, QA1, QE1, Y3);
                            X4  = mess_glyap(A, E, Ahat2, Ehat2, QA2, QE2, Y4, o);
                        else
                            [X1, Ahat1, QA1] = mess_glyap(A, E, Y1);
                            [X2, Ahat2, QA2] = mess_glyap(A, E, Y2, o);
                            X3  = mess_glyap(A, Ahat1, QA1, Y3);
                            X4  = mess_glyap(A, Ahat2, QA2, Y4, o);
                        end

                        %% compute residual
                        res2_1 = res2_glyap(A, E, Y1, X1, mess_operation_t.MESS_OP_NONE);
                        rel2_1 = res2_1 / norm(Y1);

                        res2_2 = res2_glyap(A, E, Y2, X2, o);
                        rel2_2 = res2_2 / norm(Y2);

                        res2_3 = res2_glyap(A, E, Y3, X3, mess_operation_t.MESS_OP_NONE);
                        rel2_3 = res2_3 / norm(Y3);

                        res2_4 = res2_glyap(A, E, Y4, X4, o);
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

