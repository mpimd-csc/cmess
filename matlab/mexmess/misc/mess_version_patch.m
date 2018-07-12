%MESS_VERSION_PATCH returns the patch version of C-M.E.S.S.
%
%   patch = MESS_VERSION_PATCH() returns the patch version of C-M.E.S.S. as
%   scalar double value.
%
%   See also MESS_VERSION, MESS_VERSION_VERBOSE, MESS_VERSION_MAJOR, MESS_VERSION_MINOR.
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

function ret = mess_version_patch()
    ret = mess_call('mess_version_patch');
end

