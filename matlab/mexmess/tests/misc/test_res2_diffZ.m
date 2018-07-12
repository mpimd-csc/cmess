% TEST_RES2_DIFFZ test the <a href="matlab:help res2_diffZ">res2_diffZ</a> function.
%
%   TEST_RES2_DIFFZ tests if the norm||Z1*Z1' - Z2*Z2'||_2, is computed correctly
%
%   TEST_RES2_DIFFZ properties:
%       rows            - vector, contains different number of rows of test matrices.
%       cpxZ1           - vector, contains 0 or 1. 1 means test with complex factor Z1.
%       cpxZ2           - vector, contains 0 or 1. 1 means test with complex factor Z2.
%       sparseZ1        - vector, contains 0 or 1. 1 means test with sparse factor Z1.
%       sparseZ2        - vector, contains 0 or 1. 1 means test with sparse factor Z2.
%       output          - logical, turn output on/off.
%
%   TEST_RES2_DIFFZ methods:
%       mess_test_res2_diffZ    - test function <a href="matlab:help res2_diffZ">res2_diffZ</a>
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

classdef test_res2_diffZ < matlab.unittest.TestCase

    properties
        rows = 10:10:100;
        cpxZ1 = 0:1;
        cpxZ2 = 0:1;
        sparseZ1 = 0:1;
        sparseZ2 = 0:1;
        output = false;
    end


    methods(Test)

        function mess_test_res2_diffZ(test)

            rng(1);
            for row=test.rows
                for cpx1 = test.cpxZ1
                    for cpx2 = test.cpxZ2
                        for sparse1 = test.sparseZ1
                            for sparse2 = test.sparseZ2

                                % create Z1 and Z2
                                colsZ1 = randi([1,ceil(row/4)]);
                                colsZ2 = randi([1,ceil(row/4)]);
                                if sparse1
                                    Z1 = sprand(row,colsZ1,0.5) + cpx1*1i*sprand(row,colsZ1,0.5);
                                else
                                    Z1 = rand(row,colsZ1) + cpx1*1i*rand(row,colsZ1);
                                end

                                if sparse2
                                    Z2 = sprand(row,colsZ2,0.5) + cpx2*1i*sprand(row,colsZ2,0.5);
                                else
                                    Z2 = rand(row,colsZ2) + cpx2*1i*rand(row,colsZ2);
                                end

                                % compute difference naive way
                                res2_val1 = norm(full(Z1*Z1' - Z2*Z2'));

                                % compute using the function
                                res2_val2 = res2_diffZ(Z1,Z2);

                                % print some output
                                if test.output
                                    fprintf('Z1: %d-by-%d real=%d, sparse=%d\n',size(Z1,1),size(Z1,2),isreal(Z1),issparse(Z1));
                                    fprintf('Z2: %d-by-%d real=%d, sparse=%d\n',size(Z2,1),size(Z2,2),isreal(Z2),issparse(Z2));
                                    fprintf('naive way, ||Z1*Z1''-Z2*Z2''||_2 = %e \n',res2_val1);
                                    fprintf('res2_diffZ  =  %e \n',res2_val2);
                                    fprintf('relative error = %e \n',abs(res2_val1-res2_val2)/abs(res2_val1));
                                    fprintf('absolute error = %e \n',abs(res2_val1-res2_val2));
                                end

                                % compare absolute and relative error
                                verifyLessThanOrEqual(test,abs(res2_val1-res2_val2), sqrt(eps));
                                verifyLessThanOrEqual(test,abs(res2_val1-res2_val2)/abs(res2_val1), sqrt(eps));

                            end
                        end
                    end
                end
            end
        end
    end
end





