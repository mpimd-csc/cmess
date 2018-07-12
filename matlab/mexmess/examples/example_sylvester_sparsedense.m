%EXAMPLE_SYLVESTER_SPARSEDENSE  Example for the sparse dense Sylvester equation solver <a href="matlab:help mess_sylvester_sparsedense">mess_sylvester_sparsedense</a>,
%
%   EXAMPLE_SYLVESTER_SPARSEDENSE() is an example to demostrate how to solve a sparse dense Sylvester equation.
%   This example shows the following steps:
%
%       - call <a href="matlab:help mess_sylvester_sparsedense">mess_sylvester_sparsedense</a> to solve
%
%           A X   +   X H + M = 0         (1)
%           A X   + E X H + M = 0         (2)
%           A X F + E X H + M = 0         (3)
%
%       A and E are sparse quadratic and F and H are small dense quadratic matrices.
%       M is also dense.
%
%       - call <a href="matlab:help res2_syl">res2_syl</a> to compute the residual
%       - print the relative and absolute 2-Norm residual for the sparse dense Sylvester equation
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


function ret = example_sylvester_sparsedense()

    %% generate random matrices
    rng(1);
    n = 200;
    m = 10;
    A = sprand(n,n,0.5);
    E = sprand(n,n,0.5);
    F = rand(m,m);
    H = rand(m,m);
    M = rand(n,m);

    %% solve sparse-dense Sylvester equation (1)
    X1 = mess_sylvester_sparsedense(A,H,M);

    %% solve sparse-dense Sylvester equation (2)
    X2 = mess_sylvester_sparsedense(A,E,H,M);

    %% solve sparse-dense Sylvester equation (3)
    X3 = mess_sylvester_sparsedense(A,F,E,H,M);


    %% compute residual
    res_1 = res2_syl(A,H,M,X1);
    res_2 = res2_syl(A,E,H,M,X2);
    res_3 = res2_syl(A,F,E,H,M,X3);
    rhs = norm(M);
    fprintf('absolute/relative residual of (1) = %e \t %e \n',res_1,res_1/rhs);
    fprintf('absolute/relative residual of (2) = %e \t %e \n',res_2,res_2/rhs);
    fprintf('absolute/relative residual of (3) = %e \t %e \n',res_3,res_3/rhs);

    %% return value
    ret = ~(res_1<eps^(1/4) && res_1/rhs<eps^(1/4));
    ret = ret || ~(res_2<eps^(1/4) && res_2/rhs<eps^(1/4));
    ret = ret || ~(res_3<eps^(1/4) && res_3/rhs<eps^(1/4));

end

