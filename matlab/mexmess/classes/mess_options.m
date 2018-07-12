% MESS_OPTIONS Class to control the options <a href="matlab:help mess_lradi">mess_lradi</a> and <a href="matlab:help mess_lrnm">mess_lrnm</a>.
%
%   MESS_OPTIONS properties:
%       type            - Type of equation, instance of <a href="matlab:help mess_operation_t">mess_operation_t</a>.
%       residual_method - Method for residual computation, instance of <a href="matlab:help mess_residual_t">mess_residual_t</a>.
%       nm              - Options for Newton method, instance of <a href="matlab:help mess_options_newton">mess_options_newton</a>.
%       adi             - Options for ADI method, instance of <a href="matlab:help mess_options_adi">mess_options_adi</a>.
%
%   MESS_OPTIONS methods:
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


classdef mess_options < mess_handle

    properties
        type                % Type of equation, instance of <a href="matlab:help mess_operation_t">mess_operation_t</a>.
        residual_method     % Method for residual computation, instance of <a href="matlab:help mess_residual_t">mess_residual_t</a>.
        nm                  % Options for Newton method, instance of <a href="matlab:help mess_options_newton">mess_options_newton</a>.
        adi                 % Options for ADI method, instance of <a href="matlab:help mess_options_adi">mess_options_adi</a>.
    end

    methods

        function obj = mess_options()
            % Create instance with default values.
            %   opt = MESS_OPTIONS() returns instance with default values.

            obj.type            = mess_operation_t.MESS_OP_NONE;
            obj.residual_method = mess_residual_t.MESS_RESIDUAL_INDEFINITE;
            obj.nm              = mess_options_newton();
            obj.adi             = mess_options_adi();
        end

        % setter routines
        function obj = set.type(obj,in)
            % set corresponding field of class instance
            if isa(in, 'mess_operation_t')
                obj.type = in;
            else
                error('Cannot set value to type.');
            end
        end

        function obj = set.residual_method(obj,in)
            % set corresponding field of class instance
            if isa(in, 'mess_residual_t')
                obj.residual_method = in;
            else
                error('Cannot set value to type.');
            end
        end

        function obj = set.nm(obj,in)
            % set corresponding field of class instance
            if isa(in, 'mess_options_newton')
                obj.nm = in;
            else
                error('Cannot set value to nm.');
            end
        end

        function obj = set.adi(obj,in)
            % set corresponding field of class instance
            if isa(in, 'mess_options_adi')
                obj.adi = in;
            else
                error('Cannot set value to adi.');
            end
        end

        function bol = eq(obj1,obj2)
           % Compare two class instances.
            if ~strcmp(class(obj1),class(obj2))
                error('Objects are not of the same class')
            end
            bol =   (obj1.type == obj2.type)                        && ...
                    (obj1.residual_method == obj2.residual_method)  && ...
                    (obj1.nm   == obj2.nm)                          && ...
                    (obj1.adi  == obj2.adi);
        end
    end
end
