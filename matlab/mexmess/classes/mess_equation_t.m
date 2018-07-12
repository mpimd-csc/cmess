% MESS_EQUATION_T Enumeration for different Matrix equations.
%
%   See also EQUATION.
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


classdef (Sealed = true) mess_equation_t < int8
    enumeration
        MESS_EQN_NONE(0),       % No equation, mostly for initialization.
        MESS_EQN_LYAP(1),       % Lyapunov equation.
        MESS_EQN_GLYAP(2),      % Generalized Lyapunov equation.
        MESS_EQN_RICCATI(3),    % Riccati equation.
        MESS_EQN_GRICCATI(4),   % Generalized Riccati equation.
    end
end
