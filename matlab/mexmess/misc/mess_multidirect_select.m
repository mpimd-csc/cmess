%MESS_MULTIDIRECT_SELECT selects the solver package which should be used for Multi LU decompositions.
%
%   MESS_MULTIDIRECT_SELECT(multipackage) selects the multipackage as a default Multi LU decomposition
%   of shifted matrices (A+pE). multipackage must be an instance of <a href="matlab:help mess_multidirect_t">mess_multidirect_t</a>,
%
%   See also MESS_DIRECT_LU_SELECT, MESS_DIRECT_CHOL_SELECT.
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

function mess_multidirect_select(multipackage)
    mess_call('mess_multidirect_select',multipackage);
end



