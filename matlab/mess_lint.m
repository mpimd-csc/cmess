% MESS_LINT Call generates a lint report for MEX-M.E.S.S.
%
%   This function has internal use only.
%   Adapt blacklist cell array in mess_lint to exclude files.
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


function ret = mess_lint()

    %% check version
    if verLessThan('matlab','9.2')
        error('mess_lint needs at least MATLAB 2017a');
    end

    %% exclude files from lint checking
    blacklist = {'mtx/*','mess_lint\.m'};

    %% get all m-files in subdirectories
    mess_lint_dir = fileparts(mfilename('fullpath'));
    mfiles = dir(fullfile(mess_lint_dir,'**','*.m'));

    %% create full file name
    mfilesfull = cell(0);
    k=1;
    for i=1:numel(mfiles)

        % get file name
        file = fullfile(mfiles(i).folder,mfiles(i).name);

        %% check if file is blacklisted and add if not
        onblacklist = false;
        for j = 1:numel(blacklist)
            if ~isempty(regexpi(file,blacklist{j}))
                % mark this file as blacklisted
                onblacklist = true;
                fprintf('mess_lint: File on blacklist, not included for report: %s\n',file);
            end
        end

        if ~onblacklist
            mfilesfull{k} = file;
            k=k+1;
        end
    end

    %% create lint report as printable string
    report = checkcode(mfilesfull,'-string');

    %% create lint report as cell array
    report_data = checkcode(mfilesfull);

    %% show lint report
    fprintf('============= %s: checkcode report for MEX-M.E.S.S=============\n',datestr(now));
    fprintf('\n');
    fprintf('%s',report);

    %% check report for errors / warnings
    ret = 0;
    for i = 1:numel(report_data)
        if ~isempty(report_data{i})
            ret = 1;
            return
        end
    end
end
