%MESS_CALL  Gate to C-M.E.S.S. in general internal use only.
%
%   [...] = MESS_CALL(FUNC,...) calls the function FUNC. FUNC must be
%   the function name as string.
%
%   MESS_CALL('mess_version') prints information about the C-M.E.S.S. build.
%
%   minor=MESS_CALL('mess_version_minor') returns the minor version number
%   of C-M.E.S.S.
%
%   op2 = MESS_CALL('mex_test_mess_operation_t',op1) calls the internal function
%   'mex_test_mess_operation_t' with the mess_operation_t op1 and returns
%   op2. The 'mex_test_mess_operation_t' functions simply converts a
%   mess_operation_t from MEX-M.E.S.S. to C-M.E.S.S. and back.
%   Therefore the results op2 should be equal to op1.
%
%   MESS_CALL() prints all callable functions C-M.E.S.S. via MESS_CALL.
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
