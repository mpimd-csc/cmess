%MESS_CARE  Algebraic Riccati Equation Solver
%
%   MESS_CARE solves the standard/generalized Riccati Equations
%
%   A' X   +     X A   -   X B B' X   + C' C = 0
%   A' X E + E' X A   - E' X B B' X E + C' C = 0
%
%   A is real, sparse or dense, n-by-n matrix
%   E is real, sparse or dense, regular, n-by-n matrix or empty
%   B is real, n-by-q dense matrix
%   C is real, p-by-n dense matrix
%
%   The eigenvalues of A / (A,E) must lie in the left open complex halfplane.
%   The algorithm is efficient if p and q is much smaller than n.
%
%   Z = MESS_CARE(A,[],B,C) solves the standard Riccati Equation
%   for the stabilizing solution X. The solution is returned as
%   a Cholesky factorization X=ZZ'.
%
%   Z = MESS_CARE(A,E,B,C) solves the generalized Riccati Equation
%   for the stabilizing solution X. The solution is returned as
%   a Cholesky factorization X=ZZ'.
%
%   Examples and Tests:
%       <a href="matlab:help example_care">example_care</a>, <a href="matlab:edit example_care">example_care.m</a>
%       <a href="matlab:help test_care">test_care</a>, <a href="matlab:edit test_care">test_care.m</a>
%
%   See also MESS_LYAP, RES2_RIC, MESS_EQUATION_GRICCATI, MESS_LRNM.
%
%   References:
%   <a href="matlab:disp(mess_cite('BenS10'))">[1]</a>, <a href="matlab:disp(mess_cite('BenHSetal16'))">[2]</a>, <a href="matlab:disp(mess_cite('Ham82a'))">[3]</a>,
%   <a href="matlab:disp(mess_cite('Saa09'))">[4]</a>
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

function Z = mess_care(A,E,B,C)
    Z = mess_call('mess_care',A,E,B,C);
end

