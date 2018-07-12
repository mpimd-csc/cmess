%RES2_LYAP  computes the 2-norm residual of a Lyapunov equation
%
%   RES2_LYAP computes the 2-norm residual of a Lyapunov equation for a
%   given low-rank solution factor Z.
%
%   The 2-norm residual of the following Lyapunov equations can be computed:
%
%       A ZZ' + ZZ' A' + BB'            (standard Lyapunov)
%
%       A' ZZ' + ZZ' A   + B'B          (dual Lyapunov)
%
%       A ZZ' E' + E ZZ' A'   + BB'    (generalized Lyapunov)
%
%       A' ZZ' E + E' ZZ' A   + C'C    (generalized dual Lyapunov)
%
%   op must be an instance of <a href="matlab:help mess_operation_t">mess_operation_t</a>.
%   If op is <a href="matlab:help mess_operation_t.MESS_OP_NONE">mess_operation_t.MESS_OP_NONE</a>
%   the standard/generalized Lyapunov equation residual is computed.
%
%   If op is <a href="matlab:help mess_operation_t.MESS_OP_TRANSPOSE">mess_operation_t.MESS_OP_TRANSPOSE</a>
%   the dual/generalized dual Lyapunov equation residual is computed.
%
%   Arguments:
%       A   - real, sparse or dense, n-by-n matrix
%       E   - real, sparse or dense, n-by-n matrix or empty
%       B   - real, n-by-q or q-by-n dense matrix and q<n
%       Z   - real, n-by-m dense matrix, low-rank solution factor
%       op  - instance of  <a href="matlab:help mess_operation_t">mess_operation_t</a>
%       v0  - real, n-by-1 dense vector, initial vector for <a href="matlab:help eigs">eigs</a>, optional
%
%   [res2,v] = RES2_LYAP(A,E,B,Z,op,v0) computes the 2-norm residual of a generalized/generalized dual Lyapunov equation
%   using v0 as an initial vector and returns in addition a corresponding eigenvector v.
%
%   [res2,v] = RES2_LYAP(A,E,B,Z,op) computes the 2-norm residual of a generalized/generalized dual Lyapunov equation
%   and returns in addition a corresponding eigenvector v.
%
%   [res2,v] = RES2_LYAP(A,[],B,Z,op,v0) computes the 2-norm residual of a standard/dual Lyapunov equation
%   using v0 as an initial vector and returns in addition a corresponding eigenvector v.
%
%   [res2,v] = RES2_LYAP(A,[],B,Z,op) computes the 2-norm residual of a standard/dual Lyapunov equation
%   and returns in addition a corresponding eigenvector v.
%
%   Examples and Tests:
%       <a href="matlab:help example_lyap">example_lyap</a>, <a href="matlab:edit example_lyap">example_lyap.m</a>
%       <a href="matlab:help test_lyap">test_lyap</a>, <a href="matlab:edit test_lyap">test_lyap.m</a>
%
%   See also MESS_LYAP, RES2_RIC.
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

function [res2,v] = res2_lyap(A,E,B,Z,op,v0)

    if ~isa(op,'mess_operation_t')
        error('op must be a mess_operation_t');
    end

    if op ~= mess_operation_t.MESS_OP_NONE  && op ~= mess_operation_t.MESS_OP_TRANSPOSE
        error('op must be MESS_OP_NONE or MESS_OP_TRANSPOSE');
    end

    assert(ismatrix(A) && isreal(A),'A is no real matrix');
    assert(ismatrix(B) && isreal(B),'B is no real matrix');
    assert(ismatrix(Z) && isreal(Z),'Z is no real matrix');

    if op == mess_operation_t.MESS_OP_NONE
        assert(size(B,1)>size(B,2),'B is oversized');
    else
        assert(size(B,1)<size(B,2),'B is oversized');
    end

    %% setup opts
    opts.isreal = true;
    opts.issym = true;

    %% set v0 if given
    if nargin == 6
        assert(~isempty(v0),'v0 is not empty');
        assert(isvector(v0),'v0 is no vector');
        assert(isreal(v0),'v0 is not real');
        opts.v0 = v0;
    end

    if isempty(E)
        if op == mess_operation_t.MESS_OP_NONE
            [v,d] = eigs(@(x) A*(Z*(Z'*x)) + Z*(Z'*(A'*x)) + B*(B'*x),length(A),1,'LM',opts);
        else
            [v,d] = eigs(@(x) A'*(Z*(Z'*x)) + Z*(Z'*(A*x)) + B'*(B*x),length(A),1,'LM',opts);
        end
    else
        assert(ismatrix(E) && isreal(E),'E is no real matrix');
        if op == mess_operation_t.MESS_OP_NONE
            [v,d] = eigs(@(x) A*(Z*(Z'*(E'*x))) + E*(Z*(Z'*(A'*x))) + B*(B'*x),length(A),1,'LM',opts);
        else
            [v,d] = eigs(@(x) A'*(Z*(Z'*(E*x))) + E'*(Z*(Z'*(A*x))) + B'*(B*x),length(A),1,'LM',opts);
        end
    end

    res2 = abs(d);













