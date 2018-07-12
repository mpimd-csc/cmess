% MESS_HANDLE Handle Class which hides unwanted methods from handle.
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


classdef mess_handle < handle
    methods (Hidden)
        %% Hide unwanted methods for intheritage and call superclass methods
        function obj = eq(varargin)
            obj = eq@handle(varargin{:});
        end
        function ret = findobj(varargin)
            ret = findobj@handle(varargin{:});
        end
        function prop = findprop(varargin)
            prop = findprop@handle(varargin{:});
        end
        function obj = ge(varargin)
            obj = ge@handle(varargin{:});
        end
        function obj = gt(varargin)
            obj = gt@handle(varargin{:});
        end
        function obj = addlistener(varargin)
            obj = addlistener@handle(varargin{:});
        end
        function delete(varargin)
            delete@handle(varargin{:});
        end
        function obj = le(varargin)
            obj = le@handle(varargin{:});
        end
        function obj = lt(varargin)
            obj = lt@handle(varargin{:});
        end
        function obj = ne(varargin)
            obj = ne@handle(varargin{:});
        end
        function notify(varargin)
            notify@handle(varargin{:});
        end
    end
end
