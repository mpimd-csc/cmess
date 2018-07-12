% TEST_OPTIONS test conversion of <a href=matlab:help mess_options>mess_options</a> between MEX-M.E.S.S. and C-M.E.S.S.
%
%   TEST_OPTIONS calls mex_test_options function using
%   <a href=matlab:help mess_call>mess_call</a>.
%   <a href=matlab:help mess_options>mess_options</a> are converted from MEX-M.E.S.S. to C-M.E.S.S. and back.
%   Therefore the input instance should be the same as the output instance.
%   Input instances are filled with random consistent data.
%
%   TEST_OPTIONS properties:
%       trials      - scalar, number of iterations
%       paras       - vector contains all enums of <a href=matlab:help mess_parameter_t>mess_parameter_t</a>
%       ops         - vector contains all enums of <a href=matlab:help mess_operation_t>mess_operation_t</a>
%       mems        - vector contains all enums of <a href=matlab:help mess_memusage_t>mess_memusage_t</a>
%       output      - logical, turn output on/off.
%
%   TEST_OPTIONS methods:
%       prepare_test        - fill TEST_OPTIONS instance with data
%       mess_test_options   - perform conversion of <a href=matlab:help mess_options>mess_options</a> instance
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

classdef test_options < matlab.unittest.TestCase

    properties
       trials
       paras
       ops
       mems
       output
    end

    methods(TestMethodSetup)
        function prepare_test(test)
            % number of trials
            test.trials = 30;

            % output
            test.output = false;

            % get all members of all enum types
            [test.paras,~]      = enumeration('mess_parameter_t');
            test.paras          = test.paras';
            [test.ops,~]        = enumeration('mess_operation_t');
            test.ops            = test.ops';
            [test.mems,~]       = enumeration('mess_memusage_t');
            test.mems           = test.mems';

            % random seed
            rng(1);

        end

    end


    methods(Test)

        function mess_test_options(test)

            testnum =1;
            % convert optionstype to C and back and compare
            for trial = test.trials
                for para = test.paras
                    for op = test.ops
                        for mem = test.mems

                            % set some field in options
                            o = mess_options();
                            o.type  = op;

                            % nm
                            o.nm.maxit          = randi(200);
                            o.nm.res2_tol       = rand();
                            o.nm.gpStep         = randi(200);
                            o.nm.output         = logical(rand()>0.5);
                            o.nm.singleshifts   = logical(rand()>0.5);
                            o.nm.linesearch     = logical(rand()>0.5);
                            if logical(rand()>0.5)
                                o.nm.K0 = [];
                            else
                                o.nm.K0 = rand(randi(100),randi(100));
                            end

                            % adi
                            o.adi.maxit             = randi(100);
                            o.adi.res2_tol          = rand();
                            o.adi.res2c_tol         = rand();
                            o.adi.rel_change_tol    = rand();
                            o.adi.output            = logical(rand()>0.5);
                            o.adi.memory_usage      = mem;

                            % adi shifts
                            if logical(rand()>0.5)
                                o.adi.shifts.p = zeros(1,0);
                            else
                                o.adi.shifts.p = -rand(randi(100),1);
                            end

                            o.adi.shifts.arp_p      = randi(100);
                            o.adi.shifts.arp_m      = randi(100);
                            o.adi.shifts.paratype   = para;
                            o.adi.shifts.l0         = randi(100);
                            if logical(rand()>0.5)
                                o.adi.shifts.b0 = zeros(randi(100),0);
                            else
                                o.adi.shifts.b0 = rand(randi(100),1);
                            end

                            % transfert to cmess and back
                            o2 = mess_call('mex_test_options',o);

                            % test equality
                            comp = (o == o2);
                            if test.output
                                fprintf('test=%d\t o==o2 -> %d\n',testnum,int8(comp));
                            end
                            if(~comp)
                                fprintf('FAILED\n');
                                fprintf('o=\n');
                                disp(o);
                                disp(o.nm);
                                disp(o.adi);
                                disp(o.adi.shifts);
                                fprintf('o2=\n');
                                disp(o2);
                                disp(o2.nm);
                                disp(o2.adi);
                                disp(o2.adi.shifts);
                                error('test failed options are not assumed to be equal');
                            end
                            verifyTrue(test,comp);

                            testnum = testnum+1;
                        end
                    end
                end
            end

        end

    end
end

