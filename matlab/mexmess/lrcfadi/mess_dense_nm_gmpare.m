%MESS_DENSE_NM_GMPARE dense positive and negative algebraic Riccati Equation solver
%   based on Newtons method.
%
%   MESS_DENSE_NM_GMPARE solves the standard/generalized Riccati Equations
%
%   A' X    +    X A    +/-       X G X    + Q = 0      (1.1)
%   A' X E  + E' X A    +/-    E' X G X E  + Q = 0      (1.2)
%   A  X    +    X A'   +/-       X G X    + Q = 0      (2.1)
%   A  X E' + E  X A'   +/-    E  X G X E' + Q = 0      (2.2)
%
%
%   where A, E, Q, G are dense real square matrices, Q and G symmetric positive definite.
%
%   The eigenvalues of A / (A,E) must lie in the left open complex halfplane.
%
%   A is real, dense n-by-n matrix
%   E is real, dense n-by-n matrix or empty
%   G is real, dense n-by-n matrix symmetric positive semidefinite
%   Q is real, dense n-by-n matrix symmetric positive semidefinite
%   X0 is real, dense n-by-n matrix symmetric positive semidefinite, or empty
%
%
%   [X, abs_res, relres] = MESS_DENSE_NM_GMPARE(A, [], Q, G) solves (1.1) and returns 2-norm residual.
%
%   [X, abs_res, relres] = MESS_DENSE_NM_GMPARE(A, [], Q, G, MESS_OP_NONE) solves (2.1) and returns 2-norm residual.
%
%   [X, abs_res, relres] = MESS_DENSE_NM_GMPARE(A, E, Q, G) solves (1.2) and returns 2-norm residual.
%
%   [X, abs_res, relres] = MESS_DENSE_NM_GMPARE(A, E, Q, G, MESS_OP_NONE) solves (2.1) and returns 2-norm residual.
%
%   [X, abs_res, relres] = MESS_DENSE_NM_GMPARE(A, [], Q, G, MESS_OP_TRANSPOSE) solves (1.1) and returns 2-norm residual.
%
%   [X, abs_res, relres] = MESS_DENSE_NM_GMPARE(A, E, Q, G, MESS_OP_TRANSPOSE) solves (1.2) and returns 2-norm residual.
%
%   In addition you can provide key value pairs as optional arguments to the function:
%
%   'plus'          - logical, if true positive Algebraic Riccati equation is solved otherwise negative, default false
%   'linesearch'    - logical, turn linesearch on / off, default false
%   'maxit'         - positive scalar, maximum number of newton steps, default 50
%   'norm'          - MESS_NORM_T, indicates which norm is used for residual computation, default MESS_2_NORM
%   'absres_tol'    - positive scalar, desired tolerance for absolute residual, default 1e-10
%   'relres_tol'    - positive scalar, desired tolerance for absolute residual, default 1e-11
%   'X0'            - real dense symmetric n-by-n matrix, initial guess for Newton method, default zero matrix
%
%   Example:
%       <a href="matlab:help example_dense_nm_gmpare">example_dense_nm_gmpare</a>
%
%   See also RES2_GMPARE, MESS_CARE, MESS_NORM_T.
%
%   References:
%   <a href="matlab:disp(mess_cite('Ben97b'))">[1]</a>, <a href="matlab:disp(mess_cite('Ben98'))">[2]</a>, <a href="matlab:disp(mess_cite('JonV06'))">[3]</a>
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

function [X, absres, relres] = mess_dense_nm_gmpare(A, E, Q, G, varargin)

    %% check requiered inputs arguments
    validateattributes(A, {'double'}, {'square', 'nonsparse', '2d'}, mfilename, inputname(1));
    validateattributes(Q, {'double'}, {'square', 'nonsparse', '2d'}, mfilename, inputname(3));
    validateattributes(G, {'double'}, {'square', 'nonsparse', '2d'}, mfilename, inputname(4));
    n = size(A,1);
    if ~isempty(E)
        validateattributes(E, {'double'}, {'square','nonsparse','2d'}, mfilename, inputname(2));
        assert(n==size(E,1), 'E has wrong size');
    end
    assert(n==size(Q,1), 'Q has wrong size');
    assert(n==size(G,1), 'G has wrong size');
    if ~verLessThan('matlab','8.4')
        assert(issymmetric(Q), 'Q is not symmetric');
        assert(issymmetric(G), 'G is not symmetric');
    end

    %% check optional and parameter input arguments
    p = inputParser;
    p.FunctionName = mfilename;
    addOptional(p, 'op', mess_operation_t.MESS_OP_NONE, @(x) isa(x,'mess_operation_t'));
    addParameter(p, 'plus', false, @(x) islogical(x));
    addParameter(p, 'linesearch', false, @(x) islogical(x));
    addParameter(p, 'maxit', 50, @(x) isscalar(x) && isreal(x) && x>0);
    addParameter(p, 'norm', mess_norm_t.MESS_2_NORM, @(x) isa(x,'mess_norm_t'));
    addParameter(p, 'absres_tol', 1e-10, @(x) isscalar(x) && isreal(x) && x>0);
    addParameter(p, 'relres_tol', 1e-11, @(x) isscalar(x) && isreal(x) && x>0);
    addParameter(p, 'X0', [], @(x) isscalar(x) && isreal(x) && x>0);
    parse(p,varargin{:});

    %% call mess_dense_gmpare
    [X, absres, relres] = mess_call('mess_dense_nm_gmpare', p.Results.X0, A, E, Q, G, p.Results.plus, p.Results.linesearch, p.Results.op, p.Results.maxit, ...
    p.Results.norm, p.Results.absres_tol, p.Results.relres_tol);

end


