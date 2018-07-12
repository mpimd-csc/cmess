% MESS_OPTIONS_NEWTON Class to control the options for <a href="matlab:help mess_lrnm">mess_lrnm</a>.
%
%   MESS_OPTIONS_NEWTON properties:
%       maxit           - Maximal number of iterations in Newton method.
%       res2_tol        - Stopping criteria, 2-norm of relative residual.
%       gpStep          - Galerkin Projection frequency.
%       output          - Turn output in Newton method on/off.
%       singleshifts    - Turn singleshifts in Newton method on/off.
%       linesearch      - Turn linesearch in Newton method on/off.
%       K0              - Initial Feedback.
%
%   MESS_OPTIONS_NEWTON methods:
%       eq              - Compare two instances.
%
%   See also MESS_OPTIONS_ADI, MESS_OPTIONS_ADI_SHIFTS, MESS_OPTIONS_NEWTON, MESS_LRADI, MESS_LRNM.
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


classdef mess_options_newton < mess_handle

    properties
        maxit           % Maximal number of iterations in Newton method.
        res2_tol        % Stopping criteria, 2-norm of relative residual.
        gpStep          % Galerkin Projection frequency.
        output          % Turn output in Newton method on/off.
        singleshifts    % Turn singleshifts in Newton method on/off.
        linesearch      % Turn linesearch in Newton method on/off.
        K0              % Initial Feedback.
    end

    methods

        function obj = mess_options_newton()
            % Create instance with default values.
            %   opt = MESS_OPTIONS_NEWTON() returns instance with default values.

            obj.maxit           = 30;
            obj.res2_tol        = 1e-8;
            obj.gpStep          = 0;
            obj.output          = false;
            obj.singleshifts    = false;
            obj.linesearch      = false;
            obj.K0              = [];
        end

        % setter routines
        function obj = set.maxit(obj,in)
            % set corresponding field of class instance
            if isnumeric(in) && isscalar(in) && in > 0

                % check if value is integer
                if ( abs(round(in)-in)>0 )
                    warning(['%e Input is no integer value.' ...
                        'Value will be rounded.'],in) ;
                end

                obj.maxit = round(in);

            else
                error('Cannot set value to maxit.');
            end
        end

        function obj = set.res2_tol(obj,in)
            % set corresponding field of class instance
            if isnumeric(in) && isscalar(in) && in > 0
                obj.res2_tol = in;
            else
                error('Cannot set value to res2_tol.');
            end
        end

        function obj = set.gpStep(obj,in)
            % set corresponding field of class instance
            if isnumeric(in) && isscalar(in) && in >= 0

                % check if value is integer
                if ( abs(round(in)-in)>0 )
                    warning(['%e Input is no integer value.' ...
                        'Value will be rounded.'],in) ;
                end

                obj.gpStep = round(in);

            else
                error('Cannot set value to gpStep.');
            end
        end

        function obj = set.output(obj,in)
            % set corresponding field of class instance
            if islogical(in)
                obj.output = in;
            else
                error('Cannot set value to output.');
            end
        end

        function obj = set.singleshifts(obj,in)
            % set corresponding field of class instance
            if islogical(in)
                obj.singleshifts = in;
            else
                error('Cannot set value to singleshifts.');
            end
        end

        function obj = set.linesearch(obj,in)
            % set corresponding field of class instance
            if islogical(in)
                obj.linesearch = in;
            else
                error('Cannot set value to singleshifts.');
            end
        end

        function obj = set.K0(obj,in)
            % set corresponding field of class instance
            if isnumeric(in) && ismatrix(in) && isreal(in) && ~issparse(in)
                obj.K0 = in;
            else
                error('Cannot set value to K0');
            end

        end

        function bol = eq(obj1,obj2)
           % Compare two class instances.
            if ~strcmp(class(obj1),class(obj2))
                error('Objects are not of the same class')
            end
            bol = (obj1.maxit == obj2.maxit);
            bol = bol && (abs(obj1.res2_tol - obj2.res2_tol) < eps);
            bol = bol && (obj1.gpStep == obj2.gpStep);
            bol = bol && (obj1.output == obj2.output);
            bol = bol && (obj1.singleshifts == obj2.singleshifts);
            bol = bol && (obj1.linesearch == obj2.linesearch);
            bol = bol && (norm(obj1.K0 - obj2.K0)<sqrt(eps));
            bol = bol && (isreal(obj1.K0)   == isreal(obj2.K0));
            bol = bol && (issparse(obj1.K0) == issparse(obj2.K0));
            bol = bol && (isempty(obj1.K0)  == isempty(obj2.K0));
        end
    end
end
