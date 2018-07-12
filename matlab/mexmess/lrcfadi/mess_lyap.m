%MESS_LYAP  Algebraic Lyapunov Equation Solver
%
%   MESS_LYAP solves the standard/generalized Lyapunov Equations
%
%   A X   +   X A    = - B B'
%   A X E + E X A'  = - B B'
%
%   A is real, sparse or dense, n-by-n matrix
%   E is real, sparse or dense, regular, n-by-n matrix or empty
%   B is real, n-by-q dense matrix
%
%   The eigenvalues of A / (A,E) must lie in the left open complex halfplane.
%   The algorithm is efficient if p is much smaller than n.
%
%   Z = MESS_LYAP(A,[],B) solves the standard Lyapunov Equation.
%   The solution is returned as a Cholesky factorization X=ZZ'.
%
%   Z = MESS_LYAP(A,E,B) solves the generalized Lyapunov Equation.
%   The solution is returned as a Cholesky factorization X=ZZ'.
%
%   Examples / Tests:
%       <a href="matlab:help example_lyap">example_lyap</a>, <a href="matlab:edit example_lyap">example_lyap.m</a>
%       <a href="matlab:help test_lyap">test_lyap</a>, <a href="matlab:edit test_lyap">test_lyap.m</a>
%
%   See also MESS_CARE, RES2_LYAP, MESS_EQUATION_GLYAP, MESS_LRADI.
%
%   References:
%   <a href="matlab:disp(mess_cite('Kue16'))">[1]</a>, <a href="matlab:disp(mess_cite('BenKS12'))">[2]</a>, <a href="matlab:disp(mess_cite('BenKS13a'))">[3]</a>,
%   <a href="matlab:disp(mess_cite('BenKS13'))">[4]</a>, <a href="matlab:disp(mess_cite('BenKS14b'))">[5]</a>, <a href="matlab:disp(mess_cite('Kue16'))">[6]</a>,
%   <a href="matlab:disp(mess_cite('Saa09'))">[7]</a>
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

function Z = mess_lyap(A,E,B)
    Z = mess_call('mess_lyap',A,E,B);
end

