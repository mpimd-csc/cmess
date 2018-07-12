% MESS_MEMUSAGE_T Enumeration to control the memory usage
%
%   See also MESS_OPTIONS_ADI.
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


classdef (Sealed = true) mess_memusage_t < int8
    enumeration
        MESS_MEMORY_LOW(0),     % Low memory usage.
        MESS_MEMORY_MID(1),     % Medium memory usage.
        MESS_MEMORY_HIGH(2),    % High memory usage.
    end
end
