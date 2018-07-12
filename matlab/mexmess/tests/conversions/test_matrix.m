% TEST_MATRIX test conversion of matrices between MEX-M.E.S.S. and C-M.E.S.S.
%
%   TEST_MATRIX calls mex_test_matrix function using
%   <a href=matlab:help mess_call>mess_call</a>.
%   Matrices are converted from MEX-M.E.S.S. to C-M.E.S.S. and back.
%   Therefore the input matrix should be the same as the output matrix.
%   Random matrices are used.
%
%   TEST_MATRIX properties:
%       rows            - vector, contains different number of rows of test matrices.
%       cols            - vector, contains different number of cols of test matrices.
%       cpxs            - vector, contains 0 or 1. 1 means test with complex matrices.
%       sparseA         - vector, contains 0 or 1. 1 means test with sparse matrices.
%       output          - logical, turn output on/off.
%
%   TEST_MATRIX methods:
%       mess_test_matrix    - conversion of matrices
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

classdef test_matrix < matlab.unittest.TestCase

    properties
        rows = 0:20;
        cols = 0:20;
        cpxs = 0:1;
        sparseA = 0:1;
        output = false;
    end


    methods(Test)

        function mess_test_matrix(test)

            rng(1);
            for row=test.rows
                for col=test.cols
                    for cpx = test.cpxs
                        for sp = test.sparseA

                            if sp
                                mat = sprand(row,col,0.5)+cpx*1i*sprand(row,col,0.5);
                            else
                                mat = rand(row,col)+cpx*1i*rand(row,col);
                            end
                            mat2 = mess_call('mex_test_matrix',mat);

                            if test.output
                                fprintf('class C:%s\tMEX:%s\n',class(mat2),class(mat));
                                fprintf('empty C:%d\tMEX:%d\n',isempty(mat2),isempty(mat));
                                fprintf('size C: ');   fprintf('%d ',size(mat2)); fprintf('\n');
                                fprintf('size MEX: '); fprintf('%d ',size(mat2)); fprintf('\n');
                            end

                            if(~isempty(mat))
                                nrm = norm(full(mat-mat2));
                                if test.output
                                    fprintf('norm=%e\n',nrm);
                                end

                                verifyLessThanOrEqual(test,nrm, eps);
                            end

                            if test.output
                                fprintf('class C:%s\tMEX:%s\n',class(mat2),class(mat));
                            end
                            verifyTrue(test,strcmp(class(mat2),class(mat)));
                        end
                    end
                end
            end
        end
    end
end





