% MESS_OPTIONS_ADI Class to control the options for <a href="matlab:help mess_lradi">mess_lradi</a>.
%
%   MESS_OPTIONS_ADI properties:
%       maxit           - Maximal number of iteration in <a href="matlab:help mess_lradi">mess_lradi</a>.
%       res2_tol        - Relative 2-norm stopping tolerance in <a href="matlab:help mess_lradi">mess_lradi</a>.
%       res2c_tol       - Relative 2-norm change tolerance in <a href="matlab:help mess_lradi">mess_lradi</a>.
%       rel_change_tol  - Relative Frobenius norm change tolerance in <a href="matlab:help mess_lradi">mess_lradi</a>.
%       output          - Turn output of <a href="matlab:help mess_lradi">mess_lradi</a> on/off.
%       memory_usage    - Set memory usage, instance of <a href="matlab:help mess_memusage_t">mess_memusage_t</a>.
%       shifts          - Options for shift parameters, instance of <a href="matlab:help mess_options_adi_shifts">mess_options_adi_shifts</a>.
%
%   MESS_OPTIONS_ADI methods:
%       eq              - Compare two instances.
%
%   See also MESS_OPTIONS, MESS_OPTIONS_ADI_SHIFTS, MESS_LRADI.
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


classdef mess_options_adi < mess_handle

    properties
        maxit           % Maximal number of iteration in <a href="matlab:help mess_lradi">mess_lradi</a>.
        res2_tol        % Relative 2-norm stopping tolerance in <a href="matlab:help mess_lradi">mess_lradi</a>.
        res2c_tol       % Relative 2-norm change tolerance in <a href="matlab:help mess_lradi">mess_lradi</a>.
        rel_change_tol  % Relative Frobenius norm change tolerance in <a href="matlab:help mess_lradi">mess_lradi</a>.
        output          % Turn output of <a href="matlab:help mess_lradi">mess_lradi</a> on/off.
        memory_usage    % Set memory usage, instance of <a href="matlab:help mess_memusage_t">mess_memusage_t</a>.
        shifts          % Options for shift parameters, instance of <a href="matlab:help mess_options_adi_shifts">mess_options_adi_shifts</a>.
    end

    methods

        function obj = mess_options_adi()
            % Create instance with default values.
            %   opt = MESS_OPTIONS_ADI() returns instance with default values.

            obj.maxit           = 500;
            obj.res2_tol        = 1e-10;
            obj.res2c_tol       = 1e-11;
            obj.rel_change_tol  = 1e-10;
            obj.output          = true;
            obj.memory_usage    = mess_memusage_t.MESS_MEMORY_MID;
            obj.shifts          = mess_options_adi_shifts();

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

        function obj = set.res2c_tol(obj,in)
            % set corresponding field of class instance
            if isnumeric(in) && isscalar(in) && in > 0
                obj.res2c_tol = in;
            else
                error('Cannot set value to res2c_tol.');
            end
        end

        function obj = set.rel_change_tol(obj,in)
            % set corresponding field of class instance
            if isnumeric(in) && isscalar(in) && in > 0
                obj.rel_change_tol = in;
            else
                error('Cannot set value to rel_change_tol.');
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

        function obj = set.memory_usage(obj,in)
            % set corresponding field of class instance
            if isa(in, 'mess_memusage_t')
                obj.memory_usage = in;
            else
                error('Cannot set value to memory_usage.');
            end
        end

        function obj = set.shifts(obj,in)
            % set corresponding field of class instance
            if isa(in,'mess_options_adi_shifts')
                obj.shifts = in;
            else
                error('Cannot set value to shifts.');
            end
        end

        function bol = eq(obj1,obj2)
           % Compare two class instances.
            if ~strcmp(class(obj1),class(obj2))
                error('Objects are not of the same class')
            end
            bol = (obj1.maxit == obj2.maxit);
            bol = bol && (abs(obj1.res2_tol - obj2.res2_tol) < eps);
            bol = bol && (abs(obj1.res2c_tol - obj2.res2c_tol) < eps);
            bol = bol && (abs(obj1.rel_change_tol - obj2.rel_change_tol) < eps);
            bol = bol && (obj1.output == obj2.output);
            bol = bol && (obj1.memory_usage == obj2.memory_usage);
            bol = bol && (obj1.shifts == obj2.shifts);
        end

    end
end
