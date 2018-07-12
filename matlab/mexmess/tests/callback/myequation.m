% MYEQUATION A simple implementation of <a href="matlab:help equation">equation</a> and demonstation of callback functionality.
%
%   MYEQUATION is a simple implementation of <a href="matlab:help equation">equation</a>
%   and shows how to use the callback functionality of MEX-M.E.S.S.
%   MYEQUATION implements the case when a Lyapunov or Riccati equation arises from
%   a standard linear time invariant system.
%
%     E x' = A x + B u,
%       y  = C x.
%
%   We use the properties of MYEQUATION to store all necessary data we need to apply the operations
%
%       - Matrix-Vector product with A and E.
%       - Solve a linear System with A and E.
%       - Matrix-Vector product with (A+pE) for a scalar p.
%       - Solve a linear System with (A+pE) for a scalar p.
%
%   The operations described above should be implemented using the methods
%   derived from the superclass <a href="matlab:help equation">equation</a>.
%   The methods named generate can be used for precomputing stuff,
%   for example compute a decomposition of A and use this decomposition in the
%   corresponding apply methods.
%
%   MYEQUATION properties:
%
%       dim         - state space dimension.
%       eqn_type    - instance of <a href="matlab:help mess_equation_t">mess_equation_t</a>.
%       RHS         - real, dense matrix factored right hand side of Matrix Equation
%       B           - real, dense matrix control to state matrix.
%       C           - real, dense matrix state to output matrix.
%       A           - real, sparse, stable, n-by-n matrix.
%       E           - real, sparse, regular, n-by-n matrix.
%       debug       - logical, if true more information are printed.
%
%
%   MYEQUATION Methods:
%       myequation          - Constructor.
%
%       AX_generate         - Dummy, print a message.
%       AX_apply            - Perform y=A*x or y=A'*x.
%       AX_clear            - Dummy, print a message.
%
%       EX_generate         - Dummy, print a message.
%       EX_apply            - Perform y=E*x or y=E'*x.
%       EX_clear            - Dummy, print a message.
%
%       AINV_generate       - Dummy, print a message.
%       AINV_apply          - Perform y=inv(A)*x or y=inv(A)'*x.
%       AINV_clear          - Dummy, print a message.
%
%       EINV_generate       - Dummy, print a message.
%       EINV_apply          - Perform y=inv(E)*x or y=inv(E)'*x.
%       EINV_clear          - Dummy, print a message.
%
%       ApEX_generate       - Dummy, print a message.
%       ApEX_apply          - Perform y=(A+pE)*x or y=(A+pE)'*x.
%       ApEX_clear          - Dummy, print a message.
%
%       ApEINV_generate     - Dummy, print a message.
%       ApEINV_apply        - Perform y=inv(A+pE)*x or y=inv(A+pE)'*x.
%       ApEINV_clear        - Dummy, print a message.
%
%       parameter           - Dummy, we return empty and let C-M.E.S.S. compute the shift parameters.
%
%   See also MESS_EQUATION_GRICCATI, MESS_LRADI, MESS_OPTIONS, TEST_CALLBACK_LYAPUNOV, TEST_CALLBACK_RICCATI.
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

classdef myequation <  equation

    properties
        dim             % state space dimension.
        eqn_type        % instance of <a href="matlab:help mess_equation_t">mess_equation_t</a>.
        RHS             % real, dense matrix factored right hand side of Matrix Equation
        B               % real, dense matrix control to state matrix.
        C               % real, dense matrix state to output matrix.
        A               % real, sparse, stable, n-by-n matrix.
        E               % real, sparse, regular, n-by-n matrix.
        debug = false;  % logical, if true more information are printed.
    end

    methods

        function obj = myequation(eqn_type,opt,A,E,B,C) %#ok<INUSL>
        %MYEQUATION  Contructor
        %
        %   eqn = MYEQUATION(eqn_type,opt,A,E,B,C) construct a MYEQUATION instance.
        %   Argments:
        %       eqn_type    - instance of <a href="matlab:help mess_equation_t">mess_equation_t</a>.
        %       opt         - instance of <a href="matlab:help mess_options_t">mess_options_t</a>.
        %       A           - real, sparse, stable, n-by-n matrix.
        %       E           - real, sparse, regular, n-by-n matrix.
        %       B           - real, dense matrix control to state matrix.
        %       C           - real, dense matrix state to output matrix.
        %
        %   Examples:
        %
        %   eqn = MYEQUATION(mess_equation_t.MESS_EQN_GLYAP,opt,A,E,B,[]) creates a generalized Lyapunov
        %   equation with right hand side BB'.
        %
        %   eqn = MYEQUATION(mess_equation_t.MESS_EQN_LYAP,opt,A,[],B,[]) creates a standard Lyapunov
        %   equation with right hand side BB'.
        %
        %   eqn = myequation(mess_equation_t.MESS_EQN_GRICCATI,test.opt,test.A,test.E,test.B,test.C)
        %   creates a generalized Riccati equation.
        %
        %   eqn = myequation(mess_equation_t.MESS_EQN_GRICCATI,test.opt,test.A,test.E,test.B,test.C)
        %   creates a standard Riccati equation.
        %
            if eqn_type == mess_equation_t.MESS_EQN_LYAP || eqn_type == mess_equation_t.MESS_EQN_GLYAP
                obj.A = A;
                obj.E = E;
                obj.RHS = B;
                obj.dim = length(A);
                obj.B = [];
                obj.C = [];
                obj.eqn_type = eqn_type;
            elseif eqn_type == mess_equation_t.MESS_EQN_RICCATI || eqn_type == mess_equation_t.MESS_EQN_GRICCATI
                obj.A = A;
                obj.E = E;
                obj.RHS = [];
                obj.dim = length(A);
                obj.B = B;
                obj.C = C;
                obj.eqn_type = eqn_type;
            else
                error('eqn_type %s is not allowed.\n',char(eqn_type));
            end
        end

        function AX_generate(obj,opt) %#ok<INUSD>
            if obj.debug, fprintf('call AX_generate\n'); end
        end

        function y=AX_apply(obj,opt,op,x) %#ok<INUSL>
            if obj.debug, fprintf('call AX_apply\n'); end
            if op == mess_operation_t.MESS_OP_NONE
                y = obj.A*x;
            else
                y = obj.A'*x;
            end
        end
        function AX_clear(obj)
            if obj.debug, fprintf('call AX_clear\n'); end
        end

        function EX_generate(obj,opt) %#ok<INUSD>
            if obj.debug, fprintf('call EX_generate\n'); end
        end

        function y=EX_apply(obj,opt,op,x) %#ok<INUSL>
            if obj.debug, fprintf('call EX_apply\n'); end
            if isempty(obj.E)
                % identity assumed
                y = x;
                return;
            end

            if op == mess_operation_t.MESS_OP_NONE
                y = obj.E*x;
            else
                y = obj.E'*x;
            end
        end

        function EX_clear(obj)
            if obj.debug, fprintf('call EX_clear\n'); end
        end

        function AINV_generate(obj,opt) %#ok<INUSD>
            if obj.debug, fprintf('call AINV_generate\n'); end
        end

        function y=AINV_apply(obj,opt,op,x) %#ok<INUSL>
            if obj.debug, fprintf('call AINV_apply\n'); end
            if op == mess_operation_t.MESS_OP_NONE
                y = obj.A\x;
            else
                y = obj.A'\x;
            end
        end

        function AINV_clear(obj)
            if obj.debug, fprintf('call AINV_clear\n'); end
        end

        function EINV_generate(obj,opt) %#ok<INUSD>
            if obj.debug, fprintf('call EINV_generate\n'); end
        end

        function y=EINV_apply(obj,opt,op,x) %#ok<INUSL>
            if obj.debug, fprintf('call EINV_apply\n'); end
            if isempty(obj.E)
                % identity assumed
                y = x;
                return;
            end
            if op == mess_operation_t.MESS_OP_NONE
                y = obj.E\x;
            else
                y = obj.E'\x;
            end
        end

        function EINV_clear(obj)
            if obj.debug, fprintf('call EINV_clear\n'); end
        end

        function ApEX_generate(obj,opt) %#ok<INUSD>
            if obj.debug, fprintf('call ApEX_generate\n'); end
        end

        function y=ApEX_apply(obj,opt, op, p, id, x) %#ok<INUSL>
            if obj.debug, fprintf('call ApEX_apply\n'); end
            if isempty(obj.E)
                myE = speye(obj.dim);
            else
                myE = obj.E;
            end
            if op == mess_operation_t.MESS_OP_NONE
                y = (obj.A+p*myE)*x;
            else
                y = (obj.A+p*myE)'*x;
            end
        end

        function ApEX_clear(obj)
            if obj.debug, fprintf('call ApEX_clear\n'); end
        end

        function ApEINV_generate(obj,opt,p) %#ok<INUSD>
            if obj.debug, fprintf('call ApE_generate\n'); end
        end

        function y=ApEINV_apply(obj,opt,op,p,idx,x) %#ok<INUSL>
            if obj.debug, fprintf('call ApEINV_apply\n'); end
            if isempty(obj.E)
                myE = speye(obj.dim);
            else
                myE = obj.E;
            end
            if op == mess_operation_t.MESS_OP_NONE
                y = (obj.A+p*myE)\x;
            else
                y = (obj.A+p*myE)'\x;
            end
        end

        function ApEINV_clear(obj)
            if obj.debug, fprintf('call ApEINV_clear\n'); end
        end

        function p = parameter(obj, arp_p, arp_m, B, K) %#ok<INUSD>
            if obj.debug, fprintf('call parameter\n'); end
            % return empty to make cmess compute shifts
            p = [];
            return;
        end

    end


end
