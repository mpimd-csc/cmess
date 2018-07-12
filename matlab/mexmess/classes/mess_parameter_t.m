% MESS_PARAMETER_T Enumeration for different shift parameter heuristics.
%
%   See also MESS_OPTIONS_ADI_SHIFTS.
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


classdef (Sealed = true) mess_parameter_t < int8
    enumeration
        MESS_LRCFADI_PARA_MINMAX(0),        % Compute complex parameters using the heuristic for the min-max problem.
        MESS_LRCFADI_PARA_MINMAXREAL(1),    % Compute real parameters using the heuristic for the min-max problem.
        MESS_LRCFADI_PARA_WACHSPRESS(2),    % Compute Wachspress shift parameters.
        MESS_LRCFADI_PARA_ADAPTIVE_V(3),    % Adaptive shift parameters, using V as projection space.
        MESS_LRCFADI_PARA_ADAPTIVE_Z(4),    % Adaptive shift parameters, using Z as projection space.
    end
end
