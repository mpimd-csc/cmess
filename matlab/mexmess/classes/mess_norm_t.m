% MESS_NORM_T Enumeration for different norm types.
%
%   See also MESS_DENSE_NM_GMPARE.
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


classdef (Sealed = true) mess_norm_t < int8
    enumeration
        MESS_2_NORM(0),             % The 2-norm.
        MESS_FROBENIUS_NORM(1),     % The Frobenius-norm.
        MESS_1_NORM(2),             % The 1-norm, maximum absolute column sum norm.
        MESS_INF_NORM(3)            % The inf-norm, maximum absolute row sum norm.
    end
end
