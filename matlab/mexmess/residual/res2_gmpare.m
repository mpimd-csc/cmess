%RES2_GMPARE computes the 2-norm residual of a positive/negative algebraic Riccati equation
%
%   RES2_GMPARE computes the 2-norm residual of a Riccati equation for a
%   given solution X.
%
%   The 2-norm residual of the following Riccati equations can be computed:
%
%
%   A' X    +    X A    +/-       X G X    + Q = 0      (1.1)
%   A' X E  + E' X A    +/-    E' X G X E  + Q = 0      (1.2)
%   A  X    +    X A'   +/-       X G X    + Q = 0      (2.1)
%   A  X E' + E  X A'   +/-    E  X G X E' + Q = 0      (2.2)
%
%
%   Arguments:
%       A       - real, sparse or dense, n-by-n matrix
%       E       - real, sparse or dense, n-by-n matrix or empty
%       Q       - real, n-by-n dense matrix
%       G       - real, n-by-n dense matrix
%       X       - real, n-by-n dense matrixa
%       op      - MESS_OPERATION_T
%       plus    - logical, true if positive, false if minus
%
%
%   Calling Examples:
%
%   [abs_res, relres] = RES2_GMPARE(A, [], Q, G, X, MESS_OP_TRANSPOSE, false) residual of (1.1) with minus.
%
%   [abs_res, relres] = RES2_GMPARE(A, [], Q, G, X, MESS_OP_TRANSPOSE, true)  reisdual of (1.1) with plus.
%
%   [abs_res, relres] = RES2_GMPARE(A,  E, Q, G, X, MESS_OP_NONE, false) residual of (2.2) with minus.
%
%   [abs_res, relres] = RES2_GMPARE(A,  E, Q, G, X, MESS_OP_NONE, true)  reisdual of (2.2) with plus.
%
%
%   Example:
%       <a href="matlab:help example_dense_nm_gmpare">example_dense_nm_gmpare</a>
%
%   See also RES2_GMPARE, MESS_CARE, MESS_DENSE_NM_GMPARE.
%
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

function [absres2, relres2] = res2_gmpare(A, E, Q, G, X, op, plus)

    %% check input arguments
    validateattributes(A, {'double'}, {'square', 'nonsparse', '2d'}, mfilename, inputname(1));
    validateattributes(Q, {'double'}, {'square', 'nonsparse', '2d'}, mfilename, inputname(3));
    validateattributes(G, {'double'}, {'square', 'nonsparse', '2d'}, mfilename, inputname(4));
    validateattributes(X, {'double'}, {'square', 'nonsparse', '2d'}, mfilename, inputname(4));
    n = size(A,1);
    if ~isempty(E)
        validateattributes(E, {'double'}, {'square','nonsparse','2d'}, mfilename, inputname(2));
        assert(n==size(E,1), 'E has wrong size');
    end
    assert(n==size(Q,1), 'Q has wrong size');
    assert(n==size(G,1), 'G has wrong size');
    assert(n==size(X,1), 'X has wrong size');

    if ~isa(op,'mess_operation_t')
        error('op must be a mess_operation_t');
    end
    assert(islogical(plus),'plus is no logical');

    if ~verLessThan('matlab','8.4')
        assert(issymmetric(Q), 'Q is not symmetric');
        assert(issymmetric(G), 'G is not symmetric');
        assert(issymmetric(X), 'X is not symmetric');
    end

    %% compute residual, largest eigenvalue in magnitude
    if plus
        plus_scalar = 1;
    else
        plus_scalar = -1;
    end

    opts.isreal = true;
    opts.issym = true;

    if op == mess_operation_t.MESS_OP_NONE
        dual1 = @(x) x;
        dual2 = @(x) transpose(x);
    else
        dual1 = @(x) transpose(x);
        dual2 = @(x) x;
    end

    if isempty(E)
        absres2 = abs(eigs(@(x) dual1(A)*(X*x) + X*(dual2(A)*x) + plus_scalar*(X*(G*(X*x)))  + Q*x, length(A), 1, 'LM', opts));
    else
        absres2 = abs(eigs(@(x) dual1(A)*(X*(dual2(E)*x)) +  dual1(E)*(X*(dual2(A)*x)) + plus_scalar*(dual1(E)*(X*(G*(X*(dual2(E)*x))))) + Q*x, ...
                  length(A), 1, 'LM', opts));
    end

    rhs = abs(eigs(@(x)  Q*x, length(A), 1, 'LM', opts));

    relres2 = absres2 / rhs;

end






