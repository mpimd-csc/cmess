%RES2_SYL computes the 2-norm residual of a Sylvester equation
%
%   RES2_SYL computes the 2-norm residual of a sparse-dense or dense Sylvester equation for a
%   given solution X.
%
%   The 2-norm residual of the following Sylvester equations can be computed:
%
%       A X   +   X H + M = 0         (1)
%       A X   + E X H + M = 0         (2)
%       A X F + E X H + M = 0         (3)
%
%   A is real or complex, sparse or dense n-by-n matrix
%   E is real or complex, sparse or dense n-by-n matrix
%   H is real or complex, dense m-by-m matrix
%   F is real or complex, dense m-by-m matrix
%   M is real or complex, dense n-by-m matrix
%   X is real or complex, dense n-by-m matrix
%
%   res2 = RES2_SYL(A,H,M,X) computes the 2-norm residual of (1).
%
%   res2 = RES2_SYL(A,E,H,M,X) computes the 2-norm residual of (2).
%
%   res2 = RES2_SYL(A,F,E,H,M,X) computes the 2-norm residual of (3).
%
%   Examples and Tests:
%       <a href="matlab:help example_sylvester_sparsedense">example_sylvester_sparsedense</a>, <a href="matlab:edit example_sylvester_sparsedense">example_sylvester_sparsedense.m</a>
%       <a href="matlab:help test_sylvester_sparsedense">test_sylvester_sparsedense</a>, <a href="matlab:edit test_sylvester_sparsedense">test_sylvester_sparsedense.m</a>
%
%   See also MESS_SYLVESTER_SPARSEDENSE.
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

function res2 = res2_syl(A,F,E,H,M,X)


    %% check input arguments
    switch nargin
        case 4
            A_ = A;
            F_ = [];
            E_ = [];
            H_ = F;
            M_ = E;
            X_ = H;
        case 5
            A_ = A;
            F_ = [];
            E_ = F;
            H_ = E;
            M_ = H;
            X_ = M;
        case 6
            A_ = A;
            F_ = F;
            E_ = E;
            H_ = H;
            M_ = M;
            X_ = X;
        otherwise
            error('invalid calling sequence');
    end

        assert(ismatrix(A_) && size(A_,1)==size(A_,2),'A is no square matrix');
        assert(ismatrix(H_) && ~issparse(H_) && size(H_,1)==size(H_,2),'H is no dense square matrix');
        assert(ismatrix(X_) && ~issparse(X_) && size(X_,1)==size(A_,1) && size(X_,2)==size(H_,2),'X is no dense matrix of suiteable dimesion');
        assert(ismatrix(M_) && ~issparse(M_) && all(size(M_)==size(X_)),'M is no dense matrix of suiteable dimesion');

        if nargin >= 5
            assert(ismatrix(E_) && size(E_,1)==size(E_,2) && all(size(A_)==size(E_)),'E is no square matrix of suiteable dimension');
        end

        if nargin == 6
            assert(ismatrix(F_) && ~issparse(F_) && all(size(F_)==size(H_)),'F is no dense square matrix of suiteable dimension');
        end

        %% compute the 2-Norm
        if nargin == 4
            X1 = full(A_*X_);
            X2 = full(X_*H_);
        elseif nargin == 5
            X1 = full(A_*X_);
            X2 = full(E_*X_)*H_;
        else
            X1 = full(A_*X_)*F_;
            X2 = full(E_*X_)*H_;
        end

        res2 = norm(X1+X2+M_);

    end





