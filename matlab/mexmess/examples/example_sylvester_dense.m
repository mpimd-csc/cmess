%EXAMPLE_SYLVESTER_DENSE  Example for the dense Sylvester equation solver <a href="matlab:help mess_sylvester_dense">mess_sylvester_dense</a>,
%
%   EXAMPLE_SYLVESTER_DENSE() is an example to demostrate how to solve a sparse dense Sylvester equation.
%   This example shows the following steps:
%
%       - call <a href="matlab:help mess_sylvester_dense">mess_sylvester_dense</a> to solve
%
%           A X   +   X H + M = 0         (1)
%           A X F + E X H + M = 0         (2)
%
%       A, E, F, H are dense quadratic matrices. M is also a dense matrix.
%
%       - call <a href="matlab:help res2_syl">res2_syl</a> to compute the residual
%       - print the relative and absolute 2-Norm residual for the dense Sylvester equation
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


function ret = example_sylvester_dense()

    %% generate random matrices
    rng(1);
    n = 200;
    m = 10;
    A = rand(n,n);
    E = rand(n,n);
    F = rand(m,m);
    H = rand(m,m);
    M = rand(n,m);

    %% solve dense Sylvester equation (1)
    X1 = mess_sylvester_dense(A,H,M);

    %% solve dense Sylvester equation (2)
    X2 = mess_sylvester_dense(A,F,E,H,M);

    %% compute residual
    res_1 = res2_syl(A,H,M,X1);
    res_2 = res2_syl(A,F,E,H,M,X2);
    rhs = norm(M);
    fprintf('absolute/relative residual of (1) = %e \t %e \n',res_1,res_1/rhs);
    fprintf('absolute/relative residual of (2) = %e \t %e \n',res_2,res_2/rhs);

    %% return value
    ret = ~(res_1<eps^(1/4) && res_1/rhs<eps^(1/4));
    ret = ret || ~(res_2<eps^(1/4) && res_2/rhs<eps^(1/4));

end

