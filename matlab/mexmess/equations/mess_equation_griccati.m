% MESS_EQUATION_GRICCATI Class represents a standard / generalized Riccati Equation.
%
%   MESS_EQUATION_GRICCATI represents standard / generalized Riccati Equations:
%
%   A' X    +    X A   -    X B  B' X    = -C' C
%   A  X    +    X A'  -    X C' C  X    = -B  B'
%   A' X E  + E' X A   - E' X B  B' X E  = -C' C
%   A  X E' + E  X A'  - E  X C' C  X E' = -B  B'
%
%   arising from a standard / generalized ODE System:
%
%       x' = A x + B u,
%       y  = C x.
%
%     E x' = A x + B u,
%       y  = C x.
%
%   MESS_EQUATION_GRICCATI properties:
%       E       - real, sparse, regular n-by-n matrix or empty.
%       A       - real, sparse, stable, n-by-n matrix.
%       B       - real, dense,  n-by-p or p-by-n matrix.
%       C       - real, dense,  q-by-n or n-by-q matrix.
%       dim     - state space dimension.
%
%   The spectrum of A or (A,E) must be contained in the left open
%   complex halfplane and E must be regular.
%
%   See also MESS_EQUATION_GLYAP, MESS_LRNM, MESS_OPTIONS.
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

classdef mess_equation_griccati < base_equation

    properties
        E       % real, sparse, regular n-by-n matrix or empty.
        A       % real, sparse, stable, n-by-n matrix.
        B       % real, dense,  n-by-p or p-by-n matrix.
        C       % real, dense,  q-by-n or n-by-q matrix.
        dim     % state space dimension.
    end

    methods

        function obj = mess_equation_griccati(opt, A, E, B, C)
            % Create instance from given arguments.
            %   eqn = MESS_EQUATION_GRICCATI(opt, A, [], B) creates a standard Riccati equation.
            %   If opt.type is MESS_OP_TRANSPOSE, then   A' X + X A  -  X B  B' X  = -C' C  is created.
            %   If opt.type is MESS_OP_NONE, then        A  X + X A' -  X C' C  X  = -B  B' is created.
            %
            %   eqn = MESS_EQUATION_GRICCATI(opt, A, E, B) creates a generlized Riccati equation.
            %   If opt.type is MESS_OP_TRANSPOSE, then  A' X E  + E' X A  - E' X B  B' X E  = -C' C  is created.
            %   If opt.type is MESS_OP_NONE, then       A  X E' + E  X A' - E  X C' C  X E' = -B  B' is created.
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

            % check matrix B
            assert(~isempty(C),'C is not empty');
            assert(ismatrix(C),'C is no matrix');
            assert(isreal(C),'C is no real');
            assert(isnumeric(C),'C is not numeric');
            assert(~issparse(C),'C is not dense');

            % set dimension
            obj.dim = length(A);

            % set fields
            obj.A = A;
            obj.E = E;
            obj.B = B;
            obj.C = C;

        end

    end

end
