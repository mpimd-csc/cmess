% MESS_RESIDUAL_T Enumeration for different residual methods.
%
%   See also MESS_OPTIONS.
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


classdef (Sealed = true) mess_residual_t < int8
    enumeration
        MESS_RESIDUAL_INDEFINITE(0),    % Residual via Residual Factor Updates,
        MESS_RESIDUAL_SPECTRAL(1),      % Compute Largest Eigenvalue iteratively.
    end
end
