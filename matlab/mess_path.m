% Add all required directories to the MATLAB path.
% Run this script to add all required functions and directories to the.
% MATLAB path in order to run MEX-M.E.S.S. functions and demos.
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

s=pwd;
addpath(s);
addpath(genpath(sprintf('%s/mexmess',s)));
addpath(genpath(sprintf('%s/mexmess/classes',s)));

addpath(genpath(sprintf('%s/mexmess/equations',s)));
addpath(genpath(sprintf('%s/mexmess/equations/dae1',s)));
addpath(genpath(sprintf('%s/mexmess/equations/dae2',s)));
addpath(genpath(sprintf('%s/mexmess/equations/so1',s)));
addpath(genpath(sprintf('%s/mexmess/equations/so2',s)));
addpath(genpath(sprintf('%s/mexmess/equations/std',s)));
addpath(genpath(sprintf('%s/mexmess/resdiual',s)));
addpath(genpath(sprintf('%s/mexmess/sylvester',s)));

addpath(genpath(sprintf('%s/mexmess/tests',s)));
addpath(genpath(sprintf('%s/mexmess/tests/conversions',s)));
addpath(genpath(sprintf('%s/mexmess/tests/easy',s)));
addpath(genpath(sprintf('%s/mexmess/tests/lradi',s)));
addpath(genpath(sprintf('%s/mexmess/tests/lrnm',s)));
addpath(genpath(sprintf('%s/mexmess/tests/sylvester',s)));
addpath(genpath(sprintf('%s/mexmess/tests/callback',s)));

addpath(genpath(sprintf('%s/mexmess/examples',s)));
addpath(genpath(sprintf('%s/mexmess/mtx',s)));


clear s;
