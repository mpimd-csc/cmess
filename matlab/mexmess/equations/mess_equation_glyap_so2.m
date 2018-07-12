% MESS_EQUATION_GLYAP_SO2 Class represents a Second Order system Lyapunov Equation.
%
%   MESS_EQUATION_GLYAP_SO2 represents a Second Order system:
%
%   M x''  +  D x' + K x = Bu
%
%   The Second Order System can be reformulated as a First order system:
%
%   |  I, 0 | | x'  |   |  0,  I | | x  |   | 0 |
%   |       | |     | = |        | |    | + |   | u.
%   |  0, M | | x'' |   | -K, -D | | x' |   | B |
%
%   MESS_EQUATION_GLYAP_SO2 creates an instance which represents
%   the corresponding Lyapunov equation for the above ODE.
%
%   The linearization of the second order system leads to quadratic shifts.
%   In order to avoid numerical problems one can define a lower and upper bounds
%   for the absolute value of the shift parameter p used in the
%   mess_lradi method. Shifts that not fullfill the inequality
%   lowerbound < |p| < upperbound are automatically sorted out.
%   Common choices are lowerbound=1e-5 and upperbound=1e+5.
%
%   MESS_EQUATION_GLYAP_SO2 properties:
%       M           - real, sparse, n-by-n matrix.
%       D           - real, sparse, n-by-n matrix.
%       K           - real, sparse, n-by-n matrix.
%       B           - real, dense, 2n-by-p or p-by-2n matrix.
%       lowerbound  - real, positive, scalar, lower bound for shift parameter criterion.
%       upperbound  - real, positive, scalar, upper bound for shift parameter criterion.
%       dim         - state space dimension of the underlying ODE.
%
%   See also MESS_EQUATION_GLYAP_SO1, MESS_LRADI, MESS_OPTIONS.
%
%   References:
%   <a href="matlab:disp(mess_cite('Saa09'))">[1]</a>
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

classdef mess_equation_glyap_so2 < base_equation

    properties
        M           % real, sparse, n-by-n matrix.
        D           % real, sparse, n-by-n matrix.
        K           % real, sparse, n-by-n matrix.
        B           % real, dense, 2n-by-p or p-by-2n matrix.
        lowerbound  % real, positive, scalar, lower bound for shift parameter criterion.
        upperbound  % real, positive, scalar, upper bound for shift parameter criterion.
        dim         % state space dimension of the underlying ODE.
    end

    methods

        function obj = mess_equation_glyap_so2(opt, M, D, K, B, lowerbound, upperbound)
            % Create instance from given arguments.
            %   eqn = MESS_EQUATION_GLYAP_SO2(opt, M, D, K, B, lowerbound, upperbound) returns
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

            % check matrix D
            assert(~isempty(D),'D is empty');
            assert(ismatrix(D),'D is no matrix');
            assert(isreal(D),'D is not real');
            assert(isnumeric(D),'D is not numeric');
            assert(issparse(D),'D is not sparse');
            assert(numel(unique(size(D)))==1,'D is not square');

            % check matrix K
            assert(~isempty(K),'K is empty');
            assert(ismatrix(K),'K is no matrix');
            assert(isreal(K),'K is not real');
            assert(isnumeric(K),'K is not numeric');
            assert(issparse(K),'K is not sparse');
            assert(numel(unique(size(K)))==1,'K is not square');

            %% assert same size
            assert(all(size(M)==size(D)),'M and D have different size');
            assert(all(size(D)==size(K)),'D and K have different size');

            % lowerbound and upperbound checks
            assert(isnumeric(lowerbound), 'lowerbound is not numeric');
            assert(isnumeric(upperbound), 'upperbound is not numeric');
            assert(isreal(lowerbound), 'lowerbound is not real');
            assert(isreal(upperbound), 'upperbound is not real');
            assert(lowerbound < upperbound, 'lowerbound is not smaller than upperbound');


            % check matrix B
            assert(~isempty(B),'B is not empty');
            assert(ismatrix(B),'B is no matrix');
            assert(isreal(B),'B is no real');
            assert(isnumeric(B),'B is not numeric');
            assert(~issparse(B),'B is not dense');

            % set dimension
            obj.dim = length(M);

            % set fields
            obj.M = M;
            obj.D = D;
            obj.K = K;
            obj.B = B;
            obj.lowerbound = lowerbound;
            obj.upperbound = upperbound;

        end

    end
end



