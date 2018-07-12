%EXAMPLE_LYAP  Example for the Lyapunov Equation solver <a href="matlab:help mess_lyap">mess_lyap</a>,
%
%   EXAMPLE_LYAP() is an example to demostrate how to solve a Lyapunov Equation.
%   This example shows the following steps:
%
%       - read matrices A, E, B from a file
%
%       - call <a href="matlab:help mess_lyap">mess_lyap</a> to solve A X E' + E X A' = -B B'
%       - call <a href="matlab:help res2_lyap">res2_lyap</a> to compute the residual  ||A X E' + E X A' + B B'||_2
%       - print the relative and absolute 2-Norm residual for the generalized Lyapunov equation
%
%       - call <a href="matlab:help mess_lyap">mess_lyap</a> to solve A X + X A' = -B B'
%       - call <a href="matlab:help res2_lyap">res2_lyap</a> to compute the residual  ||A X + X A' + B B'||_2
%       - print the relative and absolute 2-Norm residual for the standard Lyapunov equation
%
%   See also EXAMPLE_CARE, EXAMPLE_LRADI_RAIL, EXAMPLE_LRNM_RAIL.
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


function ret = example_lyap()

    %% read matrices
    A = mmread('@CMAKE_SOURCE_DIR@/tests/data/Rail/A.mtx');
    B = full(mmread('@CMAKE_SOURCE_DIR@/tests/data/Rail/B.mtx'));
    E = mmread('@CMAKE_SOURCE_DIR@/tests/data/Rail/E.mtx');

    %% solve generalized lyapunov equation
    Z = mess_lyap(A,E,B);

    %% compute residual
    res = res2_lyap(A,E,B,Z,mess_operation_t.MESS_OP_NONE);
    rhs = norm(B,2)^2;
    fprintf('absolute/relative residual = %e \t %e \n',res,res/rhs);

    %% return value
    ret = ~(res<eps^(1/4) && res/rhs<eps^(1/4));

    %% solve standard lyapunov equation
    Z = mess_lyap(A,[],B);

    %% compute residual
    res = res2_lyap(A,[],B,Z,mess_operation_t.MESS_OP_NONE);
    rhs = norm(B,2)^2;
    fprintf('absolute/relative residual = %e \t %e \n',res,res/rhs);

    %% return value
    ret = ret || ~(res<eps^(1/4) && res/rhs<eps^(1/4));

end
