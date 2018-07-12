%EXAMPLE_DENSE_NM_GMPARE  Example for the Riccati Equation solver <a href="matlab:help mess_dense_nm_gmpare">mess_dense_nm_gmpare</a>,
%
%   EXAMPLE_DENSE_NM_GMPARE() is an example to demostrate how to solve a positive algebraic Riccati Equation using <a href="matlab:help mess_dense_nm_gmpare">mess_dense_nm_gmpare</a>,
%   This example shows the following steps:
%
%       - read matrices A, E, B, C from a file
%       - compute G and Q (B and C are no square matrices)
%       - call <a href="matlab:help mess_dense_nm_gmpare">mess_dense_nm_gmpare</a> with optional arguments
%       - call <a href="matlab:help mess_dense_nm_gmpare">mess_dense_nm_gmpare</a> to solve A' X E  + E' X A  + E' X G X E + Q = 0
%       - call <a href="matlab:help mess_dense_nm_gmpare">mess_dense_nm_gmpare</a> to solve A  X E' + E  X A' + E  X G X E' + Q = 0
%       - compute and print the residual
%
%   See also EXAMPLE_CARE, MESS_NORM_T, RES2_GMPARE.
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


function ret = example_dense_nm_gmpare()

    ret = 0;

    %% read matrices
    A = full(mmread('@CMAKE_SOURCE_DIR@/tests/data/iss/A.mtx'));
    B = full(mmread('@CMAKE_SOURCE_DIR@/tests/data/iss/B.mtx'));
    C = full(mmread('@CMAKE_SOURCE_DIR@/tests/data/iss/C.mtx'));
    E = full(mmread('@CMAKE_SOURCE_DIR@/tests/data/iss/E.mtx'));

    %% compute Q and G
    Q = C'*C;
    G = B*B';

    for op = [mess_operation_t.MESS_OP_NONE, mess_operation_t.MESS_OP_TRANSPOSE]

        %% solve generalized riccati equation
        [X, absres, relres] = mess_dense_nm_gmpare(A, E, Q, G, op, 'plus', true, 'absres_tol', 1e-10, 'relres_tol', 1e-11);

        %% compute residual again to compare
        [absres2, relres2] = res2_gmpare(A, E, Q, G, X, op, true);

        %% print residual
        fprintf('abs./rel. 2-Norm residual from mess_dense_nm_gmpare = %e \t %e \n', absres, relres);
        fprintf('abs./rel. 2-Norm residual from res2_gmpare          = %e \t %e \n', absres2, relres2);

        ret = ret + relres2 < 1e-9;

    end

end
