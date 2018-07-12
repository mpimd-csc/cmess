% MESS_DIRECT_LUPACKAGE_T Enumeration to select the direct LU solver.
%
%   See also MESS_DIRECT_CHOLPACKAGE_T, MESS_MULTIDIRECT_T, MESS_DIRECT_LU_SELECT.
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


classdef (Sealed = true) mess_direct_lupackage_t < int8
    enumeration
        MESS_DIRECT_DEFAULT_LU(0),      % Default LU solver.
        MESS_DIRECT_SPARSE_LU(1),       % Internal LU solver.
        MESS_DIRECT_LAPACK_LU(2),       % LAPACK based LU solver.
        MESS_DIRECT_UMFPACK_LU(3),      % UMFPACK based LU solver.
        MESS_DIRECT_SUPERLU_LU(4),      % SUPERLU based LU solver.
        MESS_DIRECT_CSPARSE_LU(5),      % CSPARSE based LU solver.
        MESS_DIRECT_BANDED_LU(6),       % Banded matrices based LU solver.
        MESS_DIRECT_MKLPARDISO_LU(7)    % PARDISO based LU solver.
    end
end
