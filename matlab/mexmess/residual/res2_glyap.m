%RES2_GLYAP computes the 2-norm residual of a dense Lyapunov equation
%
%   RES2_GLYAP computes the 2-norm residual of a dense Lyapunov equation for a
%   given solution X.
%
%   The 2-norm residual of the following Lyapunov equations can be computed:
%
%       A  X    +     X A'   = - Y      (1.1)
%       A  X E' +  E  X A'   = - Y      (1.2)
%       A' X    +     X A'   = - Y      (2.1)
%       A' X E  +  E' X A'   = - Y      (2.2)
%
%   A is real, dense n-by-n matrix
%   E is real, dense n-by-n matrix
%   X is real, dense n-by-n matrix
%   Y is real, dense n-by-n matrix
%
%   res2 = RES2_GLYAP(A, Y, X) computes the 2-norm residual of (1.1).
%   res2 = RES2_GLYAP(A, Y, X, MESS_OP_NONE) computes the 2-norm residual of (1.1).
%   res2 = RES2_GLYAP(A, Y, X, MESS_OP_TRANSPOSE) computes the 2-norm residual of (1.2).
%   res2 = RES2_GLYAP(A, E, Y, X) computes the 2-norm residual of (2.1).
%   res2 = RES2_GLYAP(A, E, Y, X, MESS_OP_NONE) computes the 2-norm residual of (2.1).
%   res2 = RES2_GLYAP(A, E, Y, X, MESS_OP_TRANSPOSE) computes the 2-norm residual of (2.2).
%
%   Examples and Tests:
%       <a href="matlab:help example_glyap">example_glyap</a>, <a href="matlab:edit example_glyap">example_glyap.m</a>
%       <a href="matlab:help test_glyap">test_glyap</a>, <a href="matlab:edit test_glyap">test_glyap.m</a>
%
%   See also MESS_GLYAP.
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

function res2 = res2_glyap(A, E, Y, X, OP)

    %% check input arguments
    switch nargin
        case 3
            A_  = A;
            E_  = [];
            Y_  = E;
            X_  = Y;
            OP_ = mess_operation_t.MESS_OP_NONE;
        case 4
            if isa(X, 'mess_operation_t')
                A_  = A;
                E_  = [];
                Y_  = E;
                X_  = Y;
                OP_ = X;
            else
                A_  = A;
                E_  = E;
                Y_  = Y;
                X_  = X;
                OP_ = mess_operation_t.MESS_OP_NONE;
            end
        case 5
            A_  = A;
            E_  = E;
            Y_  = Y;
            X_  = X;
            OP_ = OP;
        otherwise
            error('invalid calling sequence');
    end

        %% check types of arguments
        assert(ismatrix(A_) && ~issparse(A_) && size(A_,1)==size(A_,2),'A is no square dense matrix');
        assert(ismatrix(Y_) && ~issparse(Y_) && size(Y_,1)==size(Y_,2),'Y is no square dense matrix');
        assert(ismatrix(X_) && ~issparse(X_) && size(X_,1)==size(X_,2),'X is no square dense matrix');
        if ~isempty(E_)
            assert(ismatrix(E_) && ~issparse(E_) && size(E_,1)==size(E_,2),'E is no square dense matrix');
        end

        %% compute the 2-Norm
        if OP_ == mess_operation_t.MESS_OP_NONE
            if isempty(E_)
                res2 = norm(A_*X_ +  X_*A_' + Y_);
            else
                res2 = norm(A_*X_*E_' + E_*X_*A_' + Y_);
            end
        else
            if isempty(E_)
                res2 = norm(A_'*X_ +  X_*A_ + Y_);
            else
                res2 = norm(A_'*X_*E_ + E_'*X_*A_ + Y_);
            end
        end
    end





