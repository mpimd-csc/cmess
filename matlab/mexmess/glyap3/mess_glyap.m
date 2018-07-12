% MESS_GLYAP Algebraic Lyapunov Equation Solver
%
%   MESS_GLYAP solves the standard/generalized Lyapunov Equations
%
%   A  X    +    X A'  = - Y      (1.1)
%   A' X    +    X A   = - Y      (1.2)
%   A  X E' + E  X A'  = - Y      (2.1)
%   A' X E  + E' X A   = - Y      (2.2)
%
%   A is real, dense, n-by-n matrix
%   E is real, dense, regular, n-by-n matrix or empty
%   Y is real, dense, symmmetrix, n-by-n matrix
%
%   The eigenvalues of A / (A, E) must fullfill the condition
%   lambda_i + lambda_j != 0.
%
%   In contrast to MESS_LYAP no stability is requiered
%   and the only condition on the right hand side is symmetry.
%   Please not that MESS_GLYAP does not check the symmetry of Y.
%
%   MESS_GLYAP uses internally a BLAS level 3 block algorithm.
%   MESS_GLYAP computes the real Schur decomposition of A / (A, E) in order to solve the equation.
%   The factorization can be returned too.
%   You can reuse the real Schur decomposition in subsequent calls of MESS_GLYAP.
%   This is usefull if you want to solve several standard / generalized Lyapunov
%   equation with different right hand sides.
%
%   The real Schur decomposition of A fullfills: A = QA Ahat QA'.
%
%   The generalized real Schur decomposition of (A,E) fullfills: A = QA Ahat QE' and E = QA Ehat QE'.
%
%   Solve:
%
%       X = MESS_GLYAP(A, [], Y) solves (1.1).
%
%       X = MESS_GLYAP(A, [], Y, MESS_OP_NONE) solves (1.1).
%
%       X = MESS_GLYAP(A, [], Y, MESS_OP_TRANSPOSE) solves (1.2).
%
%       X = MESS_GLYAP(A, E, Y) solves (2.1).
%
%       X = MESS_GLYAP(A, E, Y, MESS_OP_NONE) solves (2.1).
%
%       X = MESS_GLYAP(A, E, Y, MESS_OP_TRANSPOSE) solves (2.2).
%
%
%   Solve and get Schur decomposition:
%
%       [X, Ahat, QA] = MESS_GLYAP(A, [], Y) solves (1.1) and returns the real Schur decomposition.
%
%       [X, Ahat, QA] = MESS_GLYAP(A, [], Y, MESS_OP_NONE) solves (1.1) and returns the real Schur decomposition.
%
%       [X, Ahat, QA] = MESS_GLYAP(A, [], Y, MESS_OP_TRANSPOSE) solves (1.2)  and returns the real Schur decomposition.
%
%       [X, Ahat, Ehat, QA, QE] = MESS_GLYAP(A, E, Y) solves (2.1) and returns the generalized real Schur decomposition.
%
%       [X, Ahat, Ehat, QA, QE] = MESS_GLYAP(A, E, Y, MESS_OP_NONE) solves (2.1) and returns the generalized real Schur decomposition.
%
%       [X, Ahat, Ehat, QA, QE] = MESS_GLYAP(A, E, Y, MESS_OP_TRANSPOSE) solves (2.2) and returns the generalized real Schur decomposition.
%
%
%   Reuse Schur decomposition and solve:
%
%       X = MESS_GLYAP(A, Ahat, QA, Y) solves (1.1) using the real Schur decomposition.
%
%       X = MESS_GLYAP(A, Ahat, QA, Y, MESS_OP_NONE) solves (1.1) using the real Schur decomposition.
%
%       X = MESS_GLYAP(A, Ahat, QA, Y, MESS_OP_TRANSPOSE) solves (1.2) using real Schur decomposition.
%
%       X = MESS_GLYAP(A, E, Ahat, Ehat, QA, QE, Y) solves (2.1) using the generalized real Schur decomposition.
%
%       X = MESS_GLYAP(A, E, Ahat, Ehat, QA, QE, Y, MESS_OP_NONE) solves (2.1) using the generalized real Schur decomposition.
%
%       X = MESS_GLYAP(A, E, Ahat, Ehat, QA, QE, Y, MESS_OP_TRANSPOSE) solves (2.2) using the generalized real Schur decomposition.
%
%
%   Examples / Tests:
%       <a href="matlab:help example_glyap">example_glyap</a>, <a href="matlab:edit example_glyap">example_glyap.m</a>
%       <a href="matlab:help test_glyap">test_glyap</a>, <a href="matlab:edit test_glyap">test_glyap.m</a>
%
%   See also RES2_LYAP, MESS_LYAP, MESS_GSTEIN.
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

function varargout = mess_glyap(varargin)
    [varargout{1:nargout}] =  mess_call('mess_glyap',varargin{:});
end
