%MESS_GIT_BRANCH returns the name of the git branch which was used to configure the build.
%
%   branch=MESS_GIT_BRANCH() returns the name of the git branch which
%   was used to configure the build. branch is an char array.
%
%   See also MESS_GIT_ID.
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

function branch=mess_git_branch()
    branch=mess_call('mess_git_branch');
end

