% TEST_STATUS test conversion  <a href=matlab:help mess_status>mess_status</a>  between C-M.E.S.S. and MEX-M.E.S.S.
%
%   TEST_STATUS calls mex_test_status function using
%   <a href=matlab:help mess_call>mess_call</a>.
%   It receives an instances of <a href=matlab:help mess_status>mess_status</a> with
%   predefined values from C-M.E.S.S. Therefore it simply checks if the received values
%   are correct.
%
%   TEST_STATUS properties:
%       output      - logical, turn output on/off.
%
%   TEST_STATUS methods:
%       mess_test_status    - receive <a href=matlab:help mess_status>mess_status</a> instance and tests its values
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

classdef test_status < matlab.unittest.TestCase

    properties
        output = false;
    end

    methods(Test)

        function mess_test_status(test)
            stat = mess_call('mex_test_status');

            if test.output
                disp(stat);
                disp(stat.internal_status{1});
                disp(stat.internal_status{2});
            end


            % compare with predefined values in C
            bol = true;

            bol = bol && stat.res2_norms(1)==1;
            bol = bol && stat.res2_norms(2)==2;
            bol = bol && stat.res2_norms(3)==3;

            bol = bol && stat.rel_changes(1)==-1;
            bol = bol && stat.rel_changes(2)==-2;
            bol = bol && stat.rel_changes(3)==-3;

            bol = bol && stat.it                 == 2;
            bol = bol && stat.res2_change        == 2;
            bol = bol && stat.res2_norm           == 2;
            bol = bol && stat.res2_0             == 2;
            bol = bol && stat.rel_change         == 2;
            bol = bol && stat.stop_res2          == true;
            bol = bol && stat.stop_res2c         == true;
            bol = bol && stat.stop_rel           == true;
            bol = bol && stat.stop_user          == true;
            bol = bol && stat.time_all           == 2;
            bol = bol && stat.time_adi           == 2;
            bol = bol && stat.unstable           == 1;
            bol = bol && stat.n_internal_status  == 2;


            bol = bol && stat.internal_status{1}.it                  == 0;
            bol = bol && stat.internal_status{1}.res2_change         == 0;
            bol = bol && stat.internal_status{1}.res2_norm           == 0;
            bol = bol && stat.internal_status{1}.res2_0              == 0;
            bol = bol && stat.internal_status{1}.rel_change          == 0;
            bol = bol && stat.internal_status{1}.stop_res2           == false;
            bol = bol && stat.internal_status{1}.stop_res2c          == false;
            bol = bol && stat.internal_status{1}.stop_rel            == false;
            bol = bol && stat.internal_status{1}.stop_user           == false;
            bol = bol && stat.internal_status{1}.time_all            == 0;
            bol = bol && stat.internal_status{1}.time_adi            == 0;
            bol = bol && stat.internal_status{1}.unstable            == 0;
            bol = bol && stat.internal_status{1}.n_internal_status   == 0;


            bol = bol && stat.internal_status{2}.it                  == 1;
            bol = bol && stat.internal_status{2}.res2_change         == 1;
            bol = bol && stat.internal_status{2}.res2_norm           == 1;
            bol = bol && stat.internal_status{2}.res2_0              == 1;
            bol = bol && stat.internal_status{2}.rel_change          == 1;
            bol = bol && stat.internal_status{2}.stop_res2           == true;
            bol = bol && stat.internal_status{2}.stop_res2c          == true;
            bol = bol && stat.internal_status{2}.stop_rel            == true;
            bol = bol && stat.internal_status{2}.stop_user           == true;
            bol = bol && stat.internal_status{2}.time_all            == 1;
            bol = bol && stat.internal_status{2}.time_adi            == 1;
            bol = bol && stat.internal_status{2}.unstable            == 1;
            bol = bol && stat.internal_status{2}.n_internal_status   == 0;

            if ~bol
                fprintf('FAILED\n');
                if test.output
                    disp(stat);
                    disp(stat.internal_status{1});
                    disp(stat.internal_status{2});
                end
            else
                fprintf('PASSED\n');
            end
            verifyTrue(test,bol);

        end
    end
end
