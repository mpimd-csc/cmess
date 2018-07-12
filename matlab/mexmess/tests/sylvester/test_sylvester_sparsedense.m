% TEST_SYLVESTER_SPARSEDENSE test for <a href="matlab:help mess_sylvester_sparsedense">mess_sylvester_sparsedense</a>,
%
%   TEST_SYLVESTER_SPARSEDENSE tests the <a href="matlab:help mess_sylvester_sparsedense">mess_sylvester_sparsedense</a> solver.
%   TEST_SYLVESTER_SPARSEDENSE solves sparse dense Sylvester eqauation
%   and computes the residual using <a href="matlab:help res2_syl">res2_syl</a>.
%   TEST_SYLVESTER_SPARSEDENSE checks if the following sparse dense
%   Sylvester equations are solved properly:
%
%    A X   +   X H + M = 0         (1)
%    A X   + E X H + M = 0         (2)
%    A X F + E X H + M = 0         (3)
%
%   TEST_SYLVESTER_SPARSEDENSE properties:
%       FDMA        - sparse test matrix
%       BIPSA11     - sparse test matrix
%       BIPSE11     - sparse test matrix
%       ms          - vector, contains dimension of small dense matrix
%       cpxsA       - vector, with 0 or 1, 1 if A should be a complex converted to complex
%       cpxsE       - vector, with 0 or 1, 1 if E should be a complex converted to complex
%       cpxsF       - vector, with 0 or 1, 1 if F should be a complex converted to complex
%       cpxsH       - vector, with 0 or 1, 1 if H should be a complex converted to complex
%       cpxsM       - vector, with 0 or 1, 1 if M should be a complex converted to complex
%       seed        - random seed
%       tol         - tolerance for absolute and relative residual
%       verbose     - logical
%
%   TEST_SYLVESTER_SPARSEDENSE methods:
%       standard_sparsedense            - solves (1) and compute residual.
%       halfgeneralized_sparsedense     - solves (2) and compute residual.
%       generalized_sparsedense         - solves (3) and compute residual.
%
%   See also MESS_SYLVESTER_SPARSEDENSE, RES2_SYL.
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

classdef test_sylvester_sparsedense < matlab.unittest.TestCase

    properties
        FDMA = mmread('@CMAKE_SOURCE_DIR@/tests/data/FDM-CD/A25.mtx');
        BIPSA11 = mmread('@CMAKE_SOURCE_DIR@/tests/data/bips98_606/A11.mtx');
        BIPSE11 = mmread('@CMAKE_SOURCE_DIR@/tests/data/bips98_606/E11.mtx');
        ms = [1,2,10];
        cpxsA = [0,1];
        cpxsE = [0,1];
        cpxsF = [0,1];
        cpxsH = [0,1];
        cpxsM = [0,1];
        seed = 1;
        tol = sqrt(eps);
        verbose = true;
    end

    methods(Test)

        function standard_sparsedense(test)
            rng(test.seed);
            for cpxA = test.cpxsA
                for cpxH = test.cpxsH
                    for cpxM = test.cpxsM
                        for m = test.ms
                            %% convert A to complex if necessary
                            A = test.FDMA + cpxA*1i*test.FDMA;
                            n = size(A,1);

                            %% create random H and M matrix
                            M = rand(n,m)+cpxM*1i*rand(n,m);
                            H = rand(m,m)+cpxH*1i*rand(m,m);

                            %% solve standard sparse dense Sylvester equation
                            X = mess_sylvester_sparsedense(A,H,M);

                            %% compute residual and check
                            res = res2_syl(A,H,M,X);
                            rel = res / norm(M);
                            if(test.verbose || res > test.tol || rel > test.tol)
                                fprintf('Size(A)=%d-by-%d\n',size(A,1),size(A,2));
                                fprintf('Size(H)=%d-by-%d\n',size(H,1),size(H,2));
                                fprintf('abs. 2-Norm residual = %e\n',res);
                                fprintf('rel. 2-Norm residual = %e\n',rel);
                                fprintf('tol                  = %e\n',test.tol);
                            end

                            verifyLessThanOrEqual(test, res, test.tol);
                            verifyLessThanOrEqual(test, rel, test.tol);
                        end
                    end
                end
            end
        end

        function halfgeneralized_sparsedense(test)
            rng(test.seed);
            for cpxA = test.cpxsA
                for cpxE = test.cpxsE
                    for cpxH = test.cpxsH
                        for cpxM = test.cpxsM
                            for m = test.ms
                                %% convert A to convert if necessary
                                A = test.BIPSA11 + cpxA*1i*test.BIPSA11;
                                n = size(A,1);

                                %% convert E to convert if necessary
                                E = test.BIPSE11 + cpxA*1i*test.BIPSE11;

                                %% create random H and M matrix
                                M = rand(n,m)+cpxM*1i*rand(n,m);
                                H = rand(m,m)+cpxH*1i*rand(m,m);

                                %% solve standard sparse dense Sylvester equation
                                X = mess_sylvester_sparsedense(A,E,H,M);

                                %% compute residual and check
                                res = res2_syl(A,E,H,M,X);
                                rel = res / norm(M);
                                if(test.verbose || res > test.tol || rel > test.tol)
                                    fprintf('Size(A)=%d-by-%d\n',size(A,1),size(A,2));
                                    fprintf('Size(A)=%d-by-%d\n',size(E,1),size(E,2));
                                    fprintf('Size(H)=%d-by-%d\n',size(H,1),size(H,2));
                                    fprintf('abs. 2-Norm residual = %e\n',res);
                                    fprintf('rel. 2-Norm residual = %e\n',rel);
                                    fprintf('tol                  = %e\n',test.tol);
                                end

                                verifyLessThanOrEqual(test, res, test.tol);
                                verifyLessThanOrEqual(test, rel, test.tol);
                            end
                        end
                    end
                end
            end
        end

        function generalized_sparsedense(test)
            rng(test.seed);
            for cpxA = test.cpxsA
                for cpxE = test.cpxsE
                    for cpxF = test.cpxsF
                        for cpxH = test.cpxsH
                            for cpxM = test.cpxsM
                                for m = test.ms
                                    %% convert A to convert if necessary
                                    A = test.BIPSA11 + cpxA*1i*test.BIPSA11;
                                    n = size(A,1);

                                    %% convert E to convert if necessary
                                    E = test.BIPSE11 + cpxA*1i*test.BIPSE11;

                                    %% create random H and M matrix
                                    M = rand(n,m)+cpxM*1i*rand(n,m);
                                    H = rand(m,m)+cpxH*1i*rand(m,m);
                                    F = rand(m,m)+cpxF*1i*rand(m,m);

                                    %% solve standard sparse dense Sylvester equation
                                    X = mess_sylvester_sparsedense(A,F,E,H,M);

                                    %% compute residual and check
                                    res = res2_syl(A,F,E,H,M,X);
                                    rel = res / norm(M);
                                    if(test.verbose || res > test.tol || rel > test.tol)
                                        fprintf('Size(A)=%d-by-%d\n',size(A,1),size(A,2));
                                        fprintf('Size(A)=%d-by-%d\n',size(E,1),size(E,2));
                                        fprintf('Size(H)=%d-by-%d\n',size(H,1),size(H,2));
                                        fprintf('abs. 2-Norm residual = %e\n',res);
                                        fprintf('rel. 2-Norm residual = %e\n',rel);
                                        fprintf('tol                  = %e\n',test.tol);
                                    end

                                    verifyLessThanOrEqual(test, res, test.tol);
                                    verifyLessThanOrEqual(test, rel, test.tol);
                                end
                            end
                        end
                    end
                end
            end
        end


    end

end
