% MESS_STATUS Class contains status informations about <a href="matlab:help mess_lradi">mess_lradi</a> and <a href="matlab:help mess_lrnm">mess_lrnm</a>.
%
%   MESS_STATUS properties:
%       res2_norms          - Absolute 2-norm in every step.
%       rel_changes         - Relative change of the Iterate in every step in Frobenius norm for
%                             <a href="matlab:help mess_lradi">mess_lradi</a> / in 2-norm residual for <a href="matlab:help mess_lrnm">mess_lrnm</a>.
%       res2_norm           - Final absolute 2-norm residual.
%       res2_0              - 2-norm of right hand side.
%       it                  - Number of iteration the algorithm took.
%       res2_change         - Final relative 2-norm change of the absolute residual.
%       rel_change          - Final relative change of solution factor Z.
%       stop_res2           - Number of iteration the algorithm stops with respect to relative 2-norm residual.
%       stop_res2c          - Number of iteration the algorithm stops with respect to  2-norm relative change.
%       stop_rel            - Number of iteration the algorithm stops with respect to  2-norm relative change.
%       stop_user           - Flag if iteration was cancelled.
%       time_all            - Overall time.
%       time_adi            - Time for <a href="matlab:help mess_lradi">mess_lradi</a>.
%       unstable            - True if instable ritz values occur.
%       n_internal_status   - Number of internal MESS_STATUS instance in internal_status.
%       internal_status     - List of MESS_STATUS instances.
%
%   See also MESS_LRADI, MESS_LRNM.
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


classdef mess_status < mess_handle

    properties
        res2_norms              % Absolute 2-norm in every step.
        rel_changes             % Relative change of the Iterate in every step in Frobenius norm for
                                % <a href="matlab:help mess_lradi">mess_lradi</a> / in 2-norm residual for <a href="matlab:help mess_lrnm">mess_lrnm</a>.
        res2_norm               % Final absolute 2-norm residual.
        res2_0                  % 2-norm of right hand side.
        it                      % Number of iteration the algorithm took.
        res2_change             % Final relative 2-norm change of the absolute residual.
        rel_change              % Final relative change of solution factor Z.
        stop_res2               % Number of iteration the algorithm stops with respect to relative 2-norm residual.
        stop_res2c              % Number of iteration the algorithm stops with respect to  2-norm relative change.
        stop_rel                % Number of iteration the algorithm stops with respect to  2-norm relative change.
        stop_user               % Flag if iteration was cancelled.
        time_all                % Overall time.
        time_adi                % Time for <a href="matlab:help mess_lradi">mess_lradi</a>.
        unstable                % True if unstable ritz values occur.
        n_internal_status       % Number of internal MESS_STATUS instance in internal_status.
        internal_status         % List of MESS_STATUS instances.
    end

    methods

        function obj = mess_status()
            % Create instance with default values.
            %   opt = MESS_STATUS() returns instance with default values.
            obj.res2_norms          = zeros(0,1);
            obj.rel_changes         = zeros(0,1);
            obj.res2_norm           = 0;
            obj.res2_0              = 0;
            obj.it                  = 0;
            obj.res2_change         = 0;
            obj.rel_change          = 0;
            obj.stop_res2           = false;
            obj.stop_res2c          = false;
            obj.stop_rel            = false;
            obj.stop_user           = false;
            obj.time_all            = 0;
            obj.time_adi            = 0;
            obj.unstable            = false;
            obj.n_internal_status   = 0;
            obj.internal_status     = cell(0,1);
        end
    end
end
