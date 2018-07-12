%MESS_SYLVESTER_DENSE Dense Sylvester Equation Solver
%
%   MESS_SYLVESTER_DENSE solves the dense Sylvester equations:
%
%   A X   +   X H + M = 0         (1)
%   A X F + E X H + M = 0         (2)
%
%   A is real or complex, dense n-by-n matrix
%   E is real or complex, dense n-by-n matrix
%   H is real or complex, dense m-by-m matrix
%   F is real or complex, dense m-by-m matrix
%   M is real or complex, dense n-by-m matrix
%
%   X = MESS_SYLVESTER_DENSE(A,H,M) solves the dense Sylvester equation (1).
%
%   X = MESS_SYLVESTER_DENSE(A,F,E,H,M) solves the dense Sylvester equation (2).
%
%   Examples and Tests:
%       <a href="matlab:help example_sylvester_dense">example_sylvester_dense</a>, <a href="matlab:edit example_sylvester_dense">example_sylvester_dense.m</a>
%       <a href="matlab:help test_sylvester_dense">test_sylvester_dense</a>, <a href="matlab:edit test_sylvester_dense">test_sylvester_dense.m</a>
%
%   See also  RES2_SYL.
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

function X = mess_sylvester_dense(A,F,E,H,M)
    switch nargin
        case 3
            X = mess_call('mess_sylvester_dense',A,F,E);
        case 5
            X = mess_call('mess_sylvester_dense',A,F,E,H,M);
        otherwise
            error('invalid calling sequence');
    end
end

