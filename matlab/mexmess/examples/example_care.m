%EXAMPLE_CARE  Example for the Riccati Equation solver <a href="matlab:help mess_care">mess_care</a>,
%
%   EXAMPLE_CARE() is an example to demostrate how to solve a Riccati Equation.
%   This example shows the following steps:
%
%       - read matrices A, E, B, C from a file
%
%       - call <a href="matlab:help mess_care">mess_care</a> to solve A' X E + E' X A - E' X B B' X E  = -C' C
%       - call <a href="matlab:help res2_ric">res2_ric</a> to compute the residual  ||A' X E + E' X A - E' X B B' X E + C' C||_2
%       - print the relative and absolute 2-Norm residual for the generalized Riccati equation
%
%       - call <a href="matlab:help mess_care">mess_care</a> to solve A' X + X A - X B B' X  = -C' C
%       - call <a href="matlab:help res2_ric">res2_ric</a> to compute the residual  ||A' X + X A - X B B' X + C' C||_2
%       - print the relative and absolute 2-Norm residual for the standard Riccati equation
%
%   See also EXAMPLE_LYAP, EXAMPLE_LRADI_RAIL, EXAMPLE_LRNM_RAIL.
%
%   Author: Maximilian Behr (MPI Magdeburg, CSC)

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


function ret = example_care()

    %% read matrices
    A = mmread('@CMAKE_SOURCE_DIR@/tests/data/Rail/A.mtx');
    B = full(mmread('@CMAKE_SOURCE_DIR@/tests/data/Rail/B.mtx'));
    C = full(mmread('@CMAKE_SOURCE_DIR@/tests/data/Rail/C.mtx'));
    E = mmread('@CMAKE_SOURCE_DIR@/tests/data/Rail/E.mtx');

    %% solve generalized riccati equation
    Z = mess_care(A,E,B,C);

    %% compute residual
    res = res2_ric(A,E,B,C,Z,mess_operation_t.MESS_OP_TRANSPOSE);
    rhs = norm(C,2)^2;
    fprintf('absolute/relative residual = %e \t %e \n',res,res/rhs);

    %% return value
    ret = ~(res<eps^(1/4) && res/rhs<eps^(1/4));

    %% solve riccati equation
    Z = mess_care(A,[],B,C);

    %% compute residual
    res = res2_ric(A,[],B,C,Z,mess_operation_t.MESS_OP_TRANSPOSE);
    rhs = norm(C,2)^2;
    fprintf('absolute/relative residual = %e \t %e \n',res,res/rhs);

    %% return value
    ret = ret || ~(res<eps^(1/4) && res/rhs<eps^(1/4));

end
