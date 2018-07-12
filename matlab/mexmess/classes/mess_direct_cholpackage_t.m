% MESS_DIRECT_CHOLPACKAGE_T Enumeration to select the Cholesky solver.
%
%   See also MESS_DIRECT_LUPACKAGE_T, MESS_MULTIDIRECT_T, MESS_DIRECT_CHOL_SELECT.
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


classdef (Sealed = true) mess_direct_cholpackage_t < int8
    enumeration
        MESS_DIRECT_DEFAULT_CHOLESKY(0),    % Default cholesky solver.
        MESS_DIRECT_LAPACK_CHOLESKY(1),     % LAPACK based cholesky solver.
        MESS_DIRECT_CSPARSE_CHOLESKY(2),    % CSPARSE based cholesky solver.
        MESS_DIRECT_CHOLMOD_CHOLESKY(3)     % CHOLMOD based cholesky solver.
    end
end
