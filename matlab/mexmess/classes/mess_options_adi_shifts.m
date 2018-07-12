% MESS_OPTIONS_ADI_SHIFTS Class to control the shift parameter options for <a href="matlab:help mess_lradi">mess_lradi</a>.
%
%   MESS_OPTIONS_ADI_SHIFTS properties:
%       p           - Precomputed shifts.
%       arp_p       - Number of steps of Arnoldi process w.r.t. operator A.
%       arp_m       - Number of steps of Arnoldi process w.r.t. operator inv(A).
%       paratype    - Shift parameter type.
%       l0          - Number of shifts.
%       b0          - Initial vector for Arnoldi process.
%
%   MESS_OPTIONS_ADI_SHIFTS methods:
%       eq          - Compare two instances.
%
%   See also MESS_OPTIONS, MESS_OPTIONS_ADI, MESS_LRADI.
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


classdef mess_options_adi_shifts < mess_handle

    properties
        p           % Precomputed shifts.
        arp_p       % Number of steps of Arnoldi process w.r.t. operator A.
        arp_m       % Number of steps of Arnoldi process w.r.t. operator inv(A).
        paratype    % Shift parameter type.
        l0          % Number of shifts.
        b0          % Initial vector for Arnoldi process.
    end

    methods

        function obj = mess_options_adi_shifts()
            % Create instance with default values.
            %   opt = MESS_OPTIONS_ADI_SHIFTS() returns instance with default values.

            obj.p         = zeros(0,1);
            obj.arp_p     = 48;
            obj.arp_m     = 32;
            obj.paratype  = mess_parameter_t.MESS_LRCFADI_PARA_MINMAX;
            obj.l0        = 16;
            obj.b0        = zeros(0,1);
        end

        % setter routines
        function obj = set.p(obj,in)
            % set corresponding field of class instance
            if isnumeric(in) && ismatrix(in)

                rowvec = isrow(in);
                colvec = iscolumn(in);
                empty  = isempty(in);

                % check if in is row/colum vector
                if ~(rowvec || colvec) && ~empty
                    error('Cannot set value. Input is no vector.')
                end

                % check real part of shifts
                if any(real(in)>=0)
                    error(['Cannot set value.' ...
                        'Input has entries with nonegative realpart.']);
                end

                % set p as column vector
                if(rowvec)
                    obj.p = in';
                else
                    obj.p = in;
                end

            else
                error('Cannot set value to p');
            end
        end

        function obj = set.arp_p(obj,in)
            % set corresponding field of class instance
            if isnumeric(in) && isscalar(in) && in > 0

                % check if value is integer
                if ( abs(round(in)-in)>0 )
                    warning(['%e Input is no integer value.' ...
                        'Value will be rounded.'],in) ;
                end

                obj.arp_p = round(in);

            else
                error('Cannot set value to arp_p.');
            end
        end

        function obj = set.arp_m(obj,in)
            % set corresponding field of class instance
            if isnumeric(in) && isscalar(in) && in > 0

                % check if value is integer
                if ( abs(round(in)-in)>0 )
                    warning(['%e Input is no integer value.' ...
                        'Value will be rounded.'],in) ;
                end

                obj.arp_m = round(in);

            else
                error('Cannot set value to arp_m.');
            end
        end

        function obj = set.paratype(obj,in)
            % set corresponding field of class instance
            if isa(in, 'mess_parameter_t')
                obj.paratype = in;
            else
                error('Cannot set value to paratype.');
            end
        end

        function obj = set.l0(obj,in)
            % set corresponding field of class instance
            if isnumeric(in) && isscalar(in) && in > 0

                % check if value is integer
                if ( abs(round(in)-in)>0 )
                    warning(['%e Input is no integer value.' ...
                        'Value will be rounded.'],in) ;
                end

                obj.l0 = round(in);
            else
                error('Cannot set value to l0.');
            end
        end

        function obj = set.b0(obj,in)
            % set corresponding field of class instance
            if isnumeric(in) && ismatrix(in)

                rowvec = isrow(in);
                colvec = iscolumn(in);
                empty  = isempty(in);

                % check if in is row/colum vector
                if ~(rowvec || colvec) && ~empty
                    error('Cannot set value. Input is no vector.')
                end

                % set p as column vector
                if(rowvec)
                    obj.b0 = in';
                else
                    obj.b0 = in;
                end

            else
                error('Cannot set value to b0');
            end
        end

        function bol = eq(obj1,obj2)
            % Compare two class instances.
            if ~strcmp(class(obj1),class(obj2))
                error('Objects are not of the same class')
            end
            bol = (obj1.arp_p == obj2.arp_p);
            bol = bol && (obj1.arp_m == obj2.arp_m);
            bol = bol && (obj1.paratype == obj2.paratype);
            bol = bol && (obj1.l0 == obj2.l0);

            bol = bol && (isreal(obj1.p)   == isreal(obj2.p));
            bol = bol && (issparse(obj1.p) == issparse(obj2.p));
            bol = bol && (isempty(obj1.p)  == isempty(obj2.p));

            if(~isempty(obj1.p))
                bol = bol && (norm(obj1.p - obj2.p)<sqrt(eps));
            end

            bol = bol && (isreal(obj1.b0)   == isreal(obj2.b0));
            bol = bol && (issparse(obj1.b0) == issparse(obj2.b0));
            bol = bol && (isempty(obj1.b0)  == isempty(obj2.b0));

            if(~isempty(obj1.b0))
                bol = bol && (norm(obj1.b0 - obj2.b0)<sqrt(eps));
            end

        end
    end
end
