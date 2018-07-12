% MESS_GSTEIN Algebraic Stein Equation Solver
%
%   MESS_GSTEIN solves the standard/generalized Stein Equations
%
%   A  X A' -    X    = - Y      (1.1)
%   A' X A  -    X    = - Y      (1.2)
%   A  X A' - E  X E' = - Y      (2.1)
%   A' X A  - E' X E  = - Y      (2.2)
%
%   A is real, dense, n-by-n matrix
%   E is real, dense, regular, n-by-n matrix or empty
%   Y is real, dense, symmmetrix, n-by-n matrix
%
%   The eigenvalues of A / (A, E) must fullfill the condition
%   lambda_i *lambda_j != 1.
%
%   Please not that MESS_GSTEIN does not check the symmetry of Y.
%
%   MESS_GSTEIN uses internally a BLAS level 3 block algorithm.
%   MESS_GSTEIN computes the real Schur decomposition of A / (A, E) in order to solve the equation.
%   The factorization can be returned too.
%   You can reuse the real Schur decomposition in subsequent calls of MESS_GSTEIN.
%   This is usefull if you want to solve several standard / generalized Lyapunov
%   equation with different right hand sides.
%
%   The real Schur decomposition of A fullfills: A = QA Ahat QA'.
%
%   The generalized real Schur decomposition of (A,E) fullfills: A = QA Ahat QE' and E = QA Ehat QE'.
%
%   Solve:
%
%       X = MESS_GSTEIN(A, [], Y) solves (1.1).
%
%       X = MESS_GSTEIN(A, [], Y, MESS_OP_TRANSPOSE) solves (1.2).
%
%       X = MESS_GSTEIN(A, E, Y) solves (2.1).
%
%       X = MESS_GSTEIN(A, E, Y, MESS_OP_TRANSPOSE) solves (2.2).
%
%
%   Solve and get Schur decomposition:
%
%       [X, Ahat, QA] = MESS_GSTEIN(A, [], Y) solves (1.1) and returns the real Schur decomposition.
%
%       [X, Ahat, QA] = MESS_GSTEIN(A, [], Y, MESS_OP_TRANSPOSE) solves (1.2)  and returns the real Schur decomposition.
%
%       [X, Ahat, Ehat, QA, QE] = MESS_GSTEIN(A, E, Y) solves (2.1) and returns the generalized real Schur decomposition.
%
%       [X, Ahat, Ehat, QA, QE] = MESS_GSTEIN(A, E, Y, MESS_OP_TRANSPOSE) solves (2.2) and returns the generalized real Schur decomposition.
%
%
%   Reuse Schur decomposition and solve:
%
%       X = MESS_GSTEIN(A, Ahat, QA, Y) solves (1.1) using the real Schur decomposition.
%
%       X = MESS_GSTEIN(A, Ahat, QA, Y, MESS_OP_TRANSPOSE) solves (1.2) using real Schur decomposition.
%
%       X = MESS_GSTEIN(A, E, Ahat, Ehat, QA, QE, Y) solves (2.1) using the generalized real Schur decomposition.
%
%       X = MESS_GSTEIN(A, E, Ahat, Ehat, QA, QE, Y, MESS_OP_TRANSPOSE) solves (2.2) using the generalized real Schur decomposition.
%
%
%   Examples / Tests:
%       <a href="matlab:help example_gstein">example_gstein</a>, <a href="matlab:edit example_gstein">example_gstein.m</a>
%       <a href="matlab:help test_gstein">test_glyap</a>, <a href="matlab:edit test_gstein">test_gstein.m</a>
%
%   See also RES2_LYAP, MESS_LYAP, MESS_LRADI, MESS_GLYAP.
%
%   References:
%   <a href="matlab:disp(mess_cite('KoeS16a'))">[1]</a>
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

function varargout = mess_gstein(varargin)
    [varargout{1:nargout}] =  mess_call('mess_gstein',varargin{:});
end
