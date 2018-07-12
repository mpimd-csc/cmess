% MESS_EQUATION_GLYAP_DAE1 Class represents a Index 1 system Lyapunov Equation.
%
%   MESS_EQUATION_GLYAP_DAE1 represents a Index 1 system:
%
%   | E11, 0 | | x1' |   | A11, A12 | | x1 |   | B |
%   |        | |     | = |          | |    | + |   | u
%   | 0  , 0 | | x2' |   | A21, A22 | | x2 |   | 0 |
%
%   The DAE can be reformulated as an ODE and we solve the corresponding
%   Lyapunov equations.
%
%   MESS_EQUATION_GLYAP_DAE1 properties:
%       E11     - real, sparse, n1-by-n1 matrix.
%       A11     - real, sparse, n1-by-n1 matrix.
%       A12     - real, sparse, n1-by-n2 matrix.
%       A21     - real, sparse, n2-by-n1 matrix.
%       A22     - real, sparse, n2-by-n2 matrix.
%       B       - real, dense,  n1-by-p or p-by-n1 matrix.
%       dim     - state space dimension of the underlying ODE.
%
%   See also MESS_EQUATION_GRICCATI_DAE1, MESS_LRADI, MESS_OPTIONS.
%
%   References:
%   <a href="matlab:disp(mess_cite('morGugSW13'))">[1]</a>
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

classdef mess_equation_glyap_dae1 < base_equation

    properties
        E11     % mess_options instance.
        A11     % real, sparse, n1-by-n1 matrix.
        A12     % real, sparse, n1-by-n2 matrix.
        A21     % real, sparse, n2-by-n1 matrix.
        A22     % real, sparse, regular, n2-by-n2 matrix.
        B       % real, dense,  n1-by-p or p-by-n1 matrix.
        dim     % state space dimension of the underlying ODE.
    end

    methods

        function obj = mess_equation_glyap_dae1(opt, E11, A11, A12, A21, A22, B)
            % Create instance from given arguments.
            %   eqn = MESS_EQUATION_GLYAP_DAE1(opt, E11, A11, A12, A21, A22) returns
            %   instance from given arguments.

            % check if opt is mess_options instance
            assert(isa(opt,'mess_options'),'opt is assumed to be a mess_options_instance');

            % check matrix E11
            assert(~isempty(E11),'E11 is empty');
            assert(ismatrix(E11),'E11 is no matrix');
            assert(isreal(E11),'E11 is not real');
            assert(isnumeric(E11),'E11 is not numeric');
            assert(issparse(E11),'E11 is not sparse');
            assert(numel(unique(size(E11)))==1,'E11 is not square');

            % check matrix A11
            assert(~isempty(A11),'A11 is empty');
            assert(ismatrix(A11),'A11 is no matrix');
            assert(isreal(A11),'A11 is not real');
            assert(isnumeric(A11),'A11 is not numeric');
            assert(issparse(A11),'A11 is not sparse');
            assert(numel(unique(size(A11)))==1,'A11 is not square');

            % check matrix A12
            assert(~isempty(A12),'A12 is empty');
            assert(ismatrix(A12),'A12 is no matrix');
            assert(isreal(A12),'A12 is not real');
            assert(isnumeric(A12),'A12 is not numeric');
            assert(issparse(A12),'A12 is not sparse');

            % check matrix A21
            assert(~isempty(A21),'A21 is empty');
            assert(ismatrix(A21),'A21 is no matrix');
            assert(isreal(A21),'A21 is not real');
            assert(isnumeric(A21),'A21 is not numeric');
            assert(issparse(A21),'A21 is not sparse');

            % check matrix A22
            assert(~isempty(A22),'A22 is empty');
            assert(ismatrix(A22),'A22 is no matrix');
            assert(isreal(A22),'A22 is not real');
            assert(isnumeric(A22),'A22 is not numeric');
            assert(issparse(A22),'A22 is not sparse');

            % check matrix B
            assert(~isempty(B),'B is not empty');
            assert(ismatrix(B),'B is no matrix');
            assert(isreal(B),'B is no real');
            assert(isnumeric(B),'B is not numeric');
            assert(~issparse(B),'B is not dense');

            % set dimension
            obj.dim = length(E11);

            % set fields
            obj.E11 = E11;
            obj.A11 = A11;
            obj.A12 = A12;
            obj.A21 = A21;
            obj.A22 = A22;
            obj.B = B;

        end

    end
end


