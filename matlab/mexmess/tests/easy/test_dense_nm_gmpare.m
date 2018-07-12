% TEST_DENSE_NM_GMPARE test for <a href="matlab:help mess_dense_nm_gmpare">mess_dense_nm_gmpare</a>,
%
%   TEST_DENSE_NM_GMPARE tests the  <a href="matlab:help mess_dense_nm_gmpare">mess_dense_nm_gmpare</a> solver.
%   TEST_DENSE_NM_GMPARE solves standard / generalized negative / positive Riccati equations.
%
%   TEST_DENSE_NM_GMPARE properties:
%       A           - dense, real n-by-n matrix of Riccati equation.
%       B           - dense, real n-by-p matrix of Riccati equation (used as G = B*B').
%       C           - dense, real n-by-q matrix of Riccati equation (used as Q = C'*C).
%       E           - dense, real, regular n-by-n matrix of Riccati equation.
%
%       Please note that the eigenvalues of A / (A,E) must lie in the left open
%       complex halfplane.
%
%   TEST_DENSE_NM_GMPARE methods:
%       test_gmpare - solves Riccati equation and computes the residual.
%
%   See also RES2_GMPARE.
%
%   Author: Maximilian Behr (MPI Magdeburg, CSC), Robert Dykstra

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

classdef test_dense_nm_gmpare < matlab.unittest.TestCase

    properties
        A = full(mmread('@CMAKE_SOURCE_DIR@/tests/data/iss/A.mtx'));
        B = full(mmread('@CMAKE_SOURCE_DIR@/tests/data/iss/B.mtx'));
        C = full(mmread('@CMAKE_SOURCE_DIR@/tests/data/iss/C.mtx'));
        E = full(mmread('@CMAKE_SOURCE_DIR@/tests/data/iss/E.mtx'));
        test_tol = sqrt(eps);
        verbose = true;
        trans_ts = [mess_operation_t.MESS_OP_NONE, mess_operation_t.MESS_OP_TRANSPOSE];
        generalized_ts = [false, true];
        linesearch_ts = [false, true];
        plus_ts = [false, true];
        norm_ts = [mess_norm_t.MESS_2_NORM, mess_norm_t.MESS_FROBENIUS_NORM];
    end

    methods(Test)
        function test_gmpare(test)

            G = test.B*test.B';
            Q = test.C'*test.C;

            for trans = test.trans_ts
                for generalized = test.generalized_ts
                    for plus = test.plus_ts
                        for nrm = test.norm_ts
                            for linesearch = test.linesearch_ts

                                %% set matrices
                                A_ = test.A;
                                if generalized
                                    E_ = test.E;
                                else
                                    E_ = [];
                                end

                                %% solve riccati
                                [X, absres, relres] = mess_dense_nm_gmpare(A_, E_, Q, G, trans, 'plus', plus, 'linesearch', linesearch, 'norm', nrm);

                                %% print info and check residual
                                if test.verbose
                                    fprintf('Positive            = %d\n', plus);
                                    fprintf('LINESEARCH          = %d\n', linesearch);
                                    fprintf('GENERALIZED         = %d\n', generalized);
                                    fprintf('Size X              = %d-x-%d\n', size(X,1), size(X,1));
                                    fprintf('abs./rel. residual  = %e / %e\n', absres, relres);
                                    fprintf('\n');
                                end

                                %% check stuff
                                verifyLessThanOrEqual(test, absres, test.test_tol);
                                verifyLessThanOrEqual(test, relres, test.test_tol);

                            end
                        end
                    end
                end
            end
        end
    end
end
