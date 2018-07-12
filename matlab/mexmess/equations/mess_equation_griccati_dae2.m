% MESS_EQUATION_GRICCATI_DAE2 Class represents a Index 2 system Riccati Equation.
%
%   MESS_EQUATION_GRICCATI_DAE2 represents a Index 2 system:
%
%   | M  , 0 | | v' |   | A  , G | | v |   | B |
%   |        | |    | = |        | |   | + |   | u,
%   | 0  , 0 | | 0  |   | G' , 0 | | p |   | 0 |
%
%                y    =  C v.
%
%   The DAE can be reformulated as an ODE and we solve the corresponding
%   Lyapunov equations.
%
%   MESS_EQUATION_GRICCATI_DAE2 properties:
%       M       - real, sparse, symmetric positive definite nv-by-nv matrix.
%       A       - real, sparse, nv-by-nv matrix.
%       G       - real, sparse, full rank, nv-by-np matrix.
%       B       - real, dense, nv-by-q or q-by-nv matrix.
%       C       - real, dense, p-by-nv or nv-by-p matrix.
%       delta   - real and negativ scalar, standard value: -0.02.
%       dim     - state space dimension of the underlying ODE.
%
%   See also MESS_EQUATION_GLYAP_DAE2, MESS_LRNM, MESS_OPTIONS.
%
%   References:
%   <a href="matlab:disp(mess_cite('Wei16'))">[1]</a>, <a href="matlab:disp(mess_cite('BenSSetal13'))">[2]</a>, <a href="matlab:disp(mess_cite('BaeBSetal15'))">[3]</a>
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

classdef mess_equation_griccati_dae2 < base_equation

    properties
        M       % real, sparse, symmetric positive definite nv-by-nv matrix.
        A       % real, sparse, nv-by-nv matrix.
        G       % real, sparse, full rank, nv-by-np matrix.
        B       % real, dense, nv-by-q or q-by-nv matrix.
        C       % real, dense, p-by-nv or nv-by-p matrix.
        delta   % real and negativ scalar, standard value: -0.02.
        dim     % state space dimension of the underlying ODE.
    end

    methods

        function obj = mess_equation_griccati_dae2(opt, M, A, G, B, C, delta)
            % Create instance from given arguments.
            %   eqn = MESS_EQUATION_GRICCATI_DAE2(opt, M, A, G, B, C, delta) returns
            %   instance from given arguments.

            % check if opt is mess_options instance
            assert(isa(opt,'mess_options'),'opt is assumed to be a mess_options_instance');

            % check matrix M
            assert(~isempty(M),'M is empty');
            assert(ismatrix(M),'M is no matrix');
            assert(isreal(M),'N is not real');
            assert(isnumeric(M),'M is not numeric');
            assert(issparse(M),'M is not sparse');
            assert(numel(unique(size(M)))==1,'M is not square');

            % check matrix A
            assert(~isempty(A),'A is empty');
            assert(ismatrix(A),'A is no matrix');
            assert(isreal(A),'A is not real');
            assert(isnumeric(A),'A is not numeric');
            assert(issparse(A),'A is not sparse');
            assert(numel(unique(size(A)))==1,'A is not square');

            % check matrix G
            assert(~isempty(G),'G is empty');
            assert(ismatrix(G),'G is no matrix');
            assert(isreal(G),'G is not real');
            assert(isnumeric(G),'G is not numeric');
            assert(issparse(G),'G is not sparse');

            %% assert same size
            assert(all(size(M)==size(A)),'M and A have different size');

            % lowerbound and upperbound checks
            assert(isreal(delta), 'delta is not real');
            assert(isscalar(delta), 'delta is not scalar');

            % check matrix B
            assert(~isempty(B),'B is not empty');
            assert(ismatrix(B),'B is no matrix');
            assert(isreal(B),'B is no real');
            assert(isnumeric(B),'B is not numeric');
            assert(~issparse(B),'B is not dense');

            % check matrix C
            assert(~isempty(C),'C is not empty');
            assert(ismatrix(C),'C is no matrix');
            assert(isreal(C),'C is no real');
            assert(isnumeric(C),'C is not numeric');
            assert(~issparse(C),'C is not dense');

            % set dimension
            obj.dim = length(M);

            % set fields
            obj.M = M;
            obj.A = A;
            obj.G = G;
            obj.B = B;
            obj.C = C;
            obj.delta = delta;

        end

    end
end


