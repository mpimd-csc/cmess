% TEST_SYLVESTER_DENSE test for <a href="matlab:help mess_sylvester_dense">mess_sylvester_dense</a>,
%
%   TEST_SYLVESTER_DENSE tests the <a href="matlab:help mess_sylvester_dense">mess_sylvester_dense</a> solver.
%   TEST_SYLVESTER_DENSE solves sparse dense Sylvester eqauation
%   and computes the residual using <a href="matlab:help res2_syl">res2_syl</a>.
%   TEST_SYLVESTER_DENSE checks if the following sparse dense
%   Sylvester equations are solved properly:
%
%    A X   +   X H + M = 0         (1)
%    A X F + E X H + M = 0         (2)
%
%   TEST_SYLVESTER_DENSE properties:
%       ns          - vector, contains dimension of dense matrices A and E
%       ms          - vector, contains dimension of desne matrices F and H
%       cpxsA       - vector, with 0 or 1, 1 if A should be a complex converted to complex
%       cpxsE       - vector, with 0 or 1, 1 if E should be a complex converted to complex
%       cpxsF       - vector, with 0 or 1, 1 if F should be a complex converted to complex
%       cpxsH       - vector, with 0 or 1, 1 if H should be a complex converted to complex
%       cpxsM       - vector, with 0 or 1, 1 if M should be a complex converted to complex
%       seed        - random seed
%       tol         - tolerance for absolute and relative residual
%       verbose     - logical
%
%   TEST_SYLVESTER_DENSE methods:
%       standard_dense            - solves (1) and compute residual.
%       generalized_dense         - solves (2) and compute residual.
%
%   See also MESS_SYLVESTER_DENSE, RES2_SYL.
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

classdef test_sylvester_dense < matlab.unittest.TestCase

    properties
        ns = [1,2,10];
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

        function standard_dense(test)
            rng(test.seed);
            for cpxA = test.cpxsA
                for cpxH = test.cpxsH
                    for cpxM = test.cpxsM
                        for n = test.ns
                            for m = test.ms
                                %% create random A, H and M matrix
                                A = rand(n,n) + cpxA*1i*rand(n,n);
                                H = rand(m,m) + cpxH*1i*rand(m,m);
                                M = rand(n,m) + cpxM*1i*rand(n,m);

                                %% solve standard sparse dense Sylvester equation
                                X = mess_sylvester_dense(A,H,M);

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
        end

        function generalized_dense(test)
            rng(test.seed);
            for cpxA = test.cpxsA
                for cpxE = test.cpxsE
                    for cpxF = test.cpxsF
                        for cpxH = test.cpxsH
                            for cpxM = test.cpxsM
                                for n = test.ns
                                    for m = test.ms
                                        %% create random A, E, H, F and M matrix
                                        A = rand(n,n) + cpxA*1i*rand(n,n);
                                        E = rand(n,n) + cpxE*1i*rand(n,n);
                                        H = rand(m,m) + cpxH*1i*rand(m,m);
                                        F = rand(m,m) + cpxF*1i*rand(m,m);
                                        M = rand(n,m) + cpxM*1i*rand(n,m);

                                        %% solve standard sparse dense Sylvester equation
                                        X = mess_sylvester_dense(A,F,E,H,M);

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
end
