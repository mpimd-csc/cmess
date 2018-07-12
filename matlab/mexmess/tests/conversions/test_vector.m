% TEST_VECTOR test conversion of vectors between MEX-M.E.S.S. and C-M.E.S.S.
%
%   TEST_VECTOR calls mex_test_vector function using
%   <a href=matlab:help mess_call>mess_call</a>.
%   Vectors are converted from MEX-M.E.S.S. to C-M.E.S.S. and back.
%   Therefore the input vector should be the same as the output vector.
%   Random matrices are used.
%
%   TEST_VECTOR properties:
%       dim             - vector, different dimension of test vectors.
%       ts              - vector, contains 0,1,2, indicates if non-transposed, transposed, complex transposed is used.
%       cpxs            - vector, contains 0 or 1, indicates if real or complex matrices are used.
%       output          - logical, turn output on/off.
%
%   TEST_VECTOR methods:
%       mess_test_vector    - conversion of vectors
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

classdef test_vector < matlab.unittest.TestCase

    properties
        dims = 0:20;
        ts = 0:2;
        cpxs = 0:1;
        output = false;
    end


    methods(Test)
        function mess_test_vector(test)

            rng(1);
            %% convert vector to C and back and compare
            for dim=test.dims
                for t=test.ts
                    for cpx = test.cpxs
                        if cpx
                            v = rand(dim,1)+1i*rand(dim,1);
                        else
                            v = rand(dim,1);
                        end

                        if t==2
                            v = transpose(v);
                        elseif t==3
                            v = ctranspose(v);
                        end

                        v2 = mess_call('mex_test_vector',v);

                        if t==2
                            v = transpose(v);
                        elseif t==3
                            v = ctranspose(v);
                        end

                        if test.output
                            fprintf('class C:%s\tMEX:%s\n',class(v2),class(v));
                            fprintf('empty C:%d\tMEX:%d\n',isempty(v2),isempty(v));
                            fprintf('size C: ');   fprintf('%d ',size(v2)); fprintf('\n');
                            fprintf('size MEX: '); fprintf('%d ',size(v2)); fprintf('\n');
                        end

                        if(~isempty(v))
                            nrm = norm(v-v2);
                            if test.output
                                fprintf('norm=%e\n',nrm);
                            end
                            verifyLessThanOrEqual(test,nrm,eps);
                        end

                        if test.output
                            fprintf('class C:%s\tMEX:%s\n',class(v2),class(v));
                        end

                        verifyTrue(test,strcmp(class(v2),class(v)));

                    end
                end
            end
        end
    end
end

