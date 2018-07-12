%EXAMPLE_GLYAP  Example for the Lyapunov Equation solver <a href="matlab:help mess_glyap">mess_glyap</a>,
%
%   EXAMPLE_GLYAP() is an example to demostrate how to solve a Lyapunov Equation.
%   This example shows the following steps:
%
%       - generate random matrices A, E, Y
%       - call <a href="matlab:help mess_glyap">mess_glyap</a> to solve Lyapunov equation
%       - call <a href="matlab:help mess_glyap">mess_glyap</a> to solve Lyapunov equation again using Schur decomposition
%       - print the relative and absolute 2-Norm error
%
%   See also EXAMPLE_GSTEIN, EXAMPLE_LYAP.
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


function ret = example_glyap()

    %% generate random matrices
    ret = 0;
    rng(0);
    tol = sqrt(eps);
    n = 40;
    A=rand(n);
    E=rand(n);
    Xs = ones(n);
    Y = -(A*Xs*E'+E*Xs*A');
    YT = -(A'*Xs*E+E'*Xs*A);

    %% Non Transpose
    X1 = mess_glyap(A,E,Y);
    [X2,Ahat,Ehat,QA,QE] = mess_glyap(A,E,Y);
    X3 = mess_glyap(A,E,Ahat,Ehat,QA,QE,Y);

    rel_res1=norm(X1-Xs,1)/norm(Xs,1);
    rel_res2=norm(X2-Xs,1)/norm(Xs,1);
    rel_res3=norm(X3-Xs,1)/norm(Xs,1);

    fprintf('tol                                    = %e\n',tol);
    fprintf('relative error (1, Non Transpose)      = %e\n',rel_res1);
    fprintf('relative error (2, Non Transpose)      = %e\n',rel_res2);
    fprintf('relative error (3, Non Transpose)      = %e\n',rel_res3);

    ret = ret + (tol<rel_res1);
    ret = ret + (tol<rel_res2);
    ret = ret + (tol<rel_res3);

    %% Non Transpose 2
    X1 = mess_glyap(A,E,Y,mess_operation_t.MESS_OP_NONE);
    [X2,Ahat,Ehat,QA,QE] = mess_glyap(A,E,Y,mess_operation_t.MESS_OP_NONE);
    X3 = mess_glyap(A,E,Ahat,Ehat,QA,QE,Y,mess_operation_t.MESS_OP_NONE);

    rel_res1=norm(X1-Xs,1)/norm(Xs,1);
    rel_res2=norm(X2-Xs,1)/norm(Xs,1);
    rel_res3=norm(X3-Xs,1)/norm(Xs,1);

    fprintf('tol                                    = %e\n',tol);
    fprintf('relative error (1, Non Transpose 2)    = %e\n',rel_res1);
    fprintf('relative error (2, Non Transpose 2)    = %e\n',rel_res2);
    fprintf('relative error (3, Non Transpose 2)    = %e\n',rel_res3);

    ret = ret + (tol<rel_res1);
    ret = ret + (tol<rel_res2);
    ret = ret + (tol<rel_res3);

    %% Transpose
    X1 = mess_glyap(A,E,YT,mess_operation_t.MESS_OP_TRANSPOSE);
    [X2,Ahat,Ehat,QA,QE] = mess_glyap(A,E,YT,mess_operation_t.MESS_OP_TRANSPOSE);
    X3 = mess_glyap(A,E,Ahat,Ehat,QA,QE,YT,mess_operation_t.MESS_OP_TRANSPOSE);

    rel_res1=norm(X1-Xs,1)/norm(Xs,1);
    rel_res2=norm(X2-Xs,1)/norm(Xs,1);
    rel_res3=norm(X3-Xs,1)/norm(Xs,1);

    fprintf('tol                                    = %e\n',tol);
    fprintf('relative error (1, Non Transpose 2)    = %e\n',rel_res1);
    fprintf('relative error (2, Non Transpose 2)    = %e\n',rel_res2);
    fprintf('relative error (3, Non Transpose 2)    = %e\n',rel_res3);

    ret = ret + (tol<rel_res1);
    ret = ret + (tol<rel_res2);
    ret = ret + (tol<rel_res3);

end
