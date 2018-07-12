% MESS_EQUATION_GLYAP Class represents a standard / generalized Lyapunov Equation.
%
%   MESS_EQUATION_GLYAP represents standard / generalized Lyapunov Equations:
%
%   A' X    +    X A  = -C' C
%   A  X    +    X A' = -B  B'
%   A' X E  + E' X A  = -C' C
%   A  X E' + E  X A' = -B  B'
%
%   arising from a standard / generalized ODE System:
%
%       x' = A x + B u,
%       y  = C x.
%
%     E x' = A x + B u,
%       y  = C x.
%
%   MESS_EQUATION_GLYAP properties:
%       E       - real, sparse, regular n-by-n matrix or empty.
%       A       - real, sparse, stable, n-by-n matrix.
%       RHS     - real, dense,  n-by-p or p-by-n matrix.
%       dim     - state space dimension.
%
%   The spectrum of A or (A,E) must be contained in the left open
%   complex halfplane and E must be regular.
%
%   See also MESS_EQUATION_GRICCATI, MESS_LRADI, MESS_OPTIONS.
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

classdef mess_equation_glyap < base_equation

    properties
        E           % real, sparse, regular n-by-n matrix or empty.
        A           % real, sparse, stable, n-by-n matrix.
        RHS         % real, dense,  n-by-p or p-by-n matrix.
        dim         % state space dimension of the underlying ODE.
    end

    methods

        function obj = mess_equation_glyap(opt, A, E, B)
            % Create instance from given arguments.
            %   eqn = MESS_EQUATION_GLYAP(opt, A, [], B) creates a standard Lyapunov equation.
            %   If opt.type is MESS_OP_TRANSPOSE, then  A' X + X A  = -B' B  is created.
            %   If opt.type is MESS_OP_NONE, then       A  X + X A' = -C  C' is created.
            %
            %   eqn = MESS_EQUATION_GLYAP(opt, A, E, B) creates a generalized Lyapunov equation.
            %   If opt.type is MESS_OP_TRANSPOSE, then  A' X E  + E' X A  = -B' B  is created.
            %   If opt.type is MESS_OP_NONE, then       A  X E' + E  X A' = -B  B' is created.
            %

            % check if opt is mess_options instance
            assert(isa(opt,'mess_options'),'opt is assumed to be a mess_options_instance');

            % check matrix A
            assert(~isempty(A),'A is empty');
            assert(ismatrix(A),'A is no matrix');
            assert(isreal(A),'A is not real');
            assert(isnumeric(A),'A is not numeric');
            assert(issparse(A),'A is not sparse');
            assert(numel(unique(size(A)))==1,'A is not square');

            % check matrix E
            if ~isempty(E)
                assert(ismatrix(E),'E is no matrix');
                assert(isreal(E),'E is not real');
                assert(isnumeric(E),'E is not numeric');
                assert(issparse(E),'E is not sparse');
                assert(all(size(A)==size(E)),'A and E have different size');
            end

            % check matrix B
            assert(~isempty(B),'B is not empty');
            assert(ismatrix(B),'B is no matrix');
            assert(isreal(B),'B is no real');
            assert(isnumeric(B),'B is not numeric');
            assert(~issparse(B),'B is not dense');

            % set dimension
            obj.dim = length(A);

            % set fields
            obj.A = A;
            obj.E = E;
            obj.RHS = B;

        end

    end
end
