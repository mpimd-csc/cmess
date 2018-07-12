%EXAMPLE_LRNM_RAIL  Example for the Riccati Equation solver <a href="matlab:help mess_lrnm">mess_lrnm</a>,
%
%   EXAMPLE_LRNM_RAIL() is an example to demostrate how to solve a Riccati Equation using <a href="matlab:help mess_lrnm">mess_lrnm</a>,
%   This example shows the following steps:
%
%       - read matrices A, E, B, C from a file
%
%       - create a <a href="matlab:help mess_options">mess_options</a> instance
%       - modifiy options
%       - use <a href="matlab:help mess_equation_griccati">mess_equation_griccati</a> to create a equation
%       - call <a href="matlab:help mess_lrnm">mess_lrnm</a> to solve A' X E + E' X A - E' X B B' X E = -C' C
%       - print information from <a href="matlab:help mess_status">mess_status</a> instance
%
%   See also EXAMPLE_CARE, EXAMPLE_LYAP, EXAMPLE_LRADI_RAIL.
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


function ret = example_lrnm_rail()

    %% read matrices
    A = mmread('@CMAKE_SOURCE_DIR@/tests/data/Rail/A.mtx');
    B = full(mmread('@CMAKE_SOURCE_DIR@/tests/data/Rail/B.mtx'));
    C = full(mmread('@CMAKE_SOURCE_DIR@/tests/data/Rail/C.mtx'));
    E = mmread('@CMAKE_SOURCE_DIR@/tests/data/Rail/E.mtx');

    %% create options
    opt = mess_options();

    %% create riccati equation
    opt.type = mess_operation_t.MESS_OP_TRANSPOSE;
    opt.adi.res2_tol = 1e-10;
    opt.nm.res2_tol = 1e-5;
    opt.adi.shifts.paratype = mess_parameter_t.MESS_LRCFADI_PARA_ADAPTIVE_V;
    opt.adi.output = false;
    opt.nm.output = true;
    eqn = mess_equation_griccati(opt,A,E,B,C);

    %% solve generalized riccati equation
    [Z,status] = mess_lrnm(eqn,opt);

    %% print residual
    fprintf('Size of Low Rank Solution Factor Z: %d-x-%d\n',size(Z));
    fprintf('absolute/relative residual = %e \t %e \n',status.res2_norm,status.res2_norm/status.res2_0);

    %% return value
    ret = ~(status.res2_norm<eps^(1/4) && status.res2_norm/status.res2_0<eps^(1/4));

    %% call again to test the short calling sequence, and compare the solutions
    opt.adi.output = false;
    opt.nm.output = false;
    Z2 = mess_lrnm(eqn,opt);
    ret = ret || ~(res2_diffZ(Z,Z2)<eps^(1/2) );

end
