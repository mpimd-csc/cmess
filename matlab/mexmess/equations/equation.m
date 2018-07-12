% EQUATION Class represents an abstract interface for user defined equations.
%
%   EQUATION represents an abstract interface for user defined equations
%   using callback functionality.
%
%   EQUATION properties:
%       dim         - int, state space dimension.
%       eqn_type    - Type of Equation.
%       B           - Control to state matrix.
%       C           - State to output matrix.
%       RHS         - Matrix, right hand side of equation.
%
%   EQUATION methods:
%
%       AX_generate     -   Prepares the application of A to a right hand side vector or
%                           matrix. It is called onces before x = Ay is called the first time.
%       AX_apply        -   Applies A to a right hand side. It has to return the result of y = Ax.
%       AX_clear        -   Finalizes the application of the A operator. It is called after y = Ax is applied the last time.
%
%       EX_generate     -   Prepares the application of E to a right hand side vector or
%                           matrix. It is called onces before x = Ey is called the first time.
%       EX_apply        -   Applies E to a right hand side. It has to return the result of y = Ex.
%       EX_clear        -   Finalizes the application of the E operator. It is called after y = Ex is applied the last time.
%
%       AINV_generate   -   Prepares the application of the inverse of the inv(A) operator on a
%                           given right hand side vector or matrix. It is called onces before y = inv(A)x is called the first time.
%       AINV_apply      -   Applies the inverse of A to a right hand side. It has to return the solution of Ax = y.
%       AINV_clear      -   Finalizes the application of the inv(A) operator. It is called after y = inv(A)x is applied the last time.
%
%       EINV_generate   -   Prepares the application of the inverse of the E operator on a
%                           given right hand side vector or matrix. It is called onces before y = inv(E)x is called the first time.
%       EINV_apply      -   Applies the inverse of E to a right hand side. It has to return the solution of Ex = y.
%       EINV_clear      -   Finalizes the application of the inv(E) operator. It is called after y = inv(E)x is applied the last time.
%
%       ApEX_generate   -   Prepares the application of the A+pE operator. It is called onces before ApEX_apply is called the first time.
%       ApEX_apply      -   Applies function the operator of A+pE. It has to return y=(A+pE)x.
%       ApEX_clear      -   Finalizes the application of the inverse of the A+pE operator. It is called after ApEX_apply is used the last time.
%
%       ApEINV_generate -   Prepares the application of the inv(A+pE) operator.
%                           It is called onces before ApEINV_apply is called the first time.
%       ApEINV_apply    -   Applies the inverse of inv(A+pE) to a right hand side. It has to return the solution of y = inv(A+pE)x.
%       ApEINV_clear    -   Finalizes the application of the inverse of the inv(A+pE) operator. It is called after ApEINV_apply is used the last time.
%
%       parameter       -   The parmeter function has to return the shift parameters. If [] is returned shift paramter will be automatically determined.
%                           The Shift parameter strategy is determined by the options structure.
%   Examples / Tests:
%       <a href="matlab:help my_equation">my_equation</a>, <a href="matlab:edit my_equation">my_equation.m</a>
%       <a href="matlab:help test_callback_lyapunov">test_callback_lyapunov</a>, <a href="matlab:edit test_callback_lyapunov">test_callback_lyapunov.m</a>
%       <a href="matlab:help test_callback_riccati">test_callback_riccati</a>, <a href="matlab:edit test_callback_riccati">test_callback_riccati.m</a>
%
%   See also MESS_EQUATION_T, MESS_LRADI, MESS_LRNM, MESS_OPTIONS, MYEQUATION, TEST_CALLBACK_LYAPUNOV, TEST_CALLBACK_RICCATI.
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


classdef (Abstract) equation < base_equation

    properties(Abstract)
        dim         % int, state space dimension
        eqn_type    % Type of Equation
        B           % control to state matrix
        C           % state to output matrix
        RHS         % matrix, right hand side of equation
    end

    methods(Abstract)
        AX_generate(obj,opt)
        AX_apply(obj,opt,op,x)
        AX_clear(obj)

        EX_generate(obj,opt)
        EX_apply(obj,opt,op,x)
        EX_clear(obj)

        AINV_generate(obj,opt)
        AINV_apply(obj,opt,op,x)
        AINV_clear(obj)

        EINV_generate(obj,opt)
        EINV_apply(obj,opt,op,x)
        EINV_clear(obj)

        ApEX_generate(obj,opt)
        ApEX_apply(obj,opt,op,p, id, x)
        ApEX_clear(obj)

        ApEINV_generate(obj,opt,p)
        ApEINV_apply(obj,opt,op,p,idx,x)
        ApEINV_clear(obj)

        parameter(obj, arp_p, arp_m, B, K)
    end
end
