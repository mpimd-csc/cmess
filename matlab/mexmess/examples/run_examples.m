%RUN_EXAMPLES  Runs examples for MEX-M.E.S.S.
%
%   RUN_EXAMPLES() runs all examples for MEX-M.E.S.S.
%
%   See also EXAMPLE_CARE, EXAMPLE_LRADI_RAIL, EXAMPLE_LRNM_RAIL, EXAMPLE_LYAP, EXAMPLE_SYLVESTER_SPARSEDENSE, EXAMPLE_GLYAP, EXAMPLE_GSTEIN,
%   EXAMPLE_SYLVESTER_DENSE.
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

function ret = run_examples()

    fprintf('RUN EXAMPLE: example_lyap.m\n\n');
    ret = example_lyap();
    fprintf('-------------------------------------------------------------\n\n');
    if ret, return; end

    fprintf('RUN EXAMPLE: example_care.m\n\n');
    ret = example_care();
    fprintf('-------------------------------------------------------------\n\n');
    if ret, return; end

    fprintf('RUN EXAMPLE: example_lradi_rail.m\n\n');
    ret = example_lradi_rail();
    fprintf('-------------------------------------------------------------\n\n');
    if ret, return; end

    fprintf('RUN EXAMPLE: example_lrnm_rail.m\n\n');
    ret = example_lrnm_rail();
    fprintf('-------------------------------------------------------------\n\n');
    if ret, return; end

    fprintf('RUN EXAMPLE: example_sylvester_sparsedense.m\n\n');
    ret = example_sylvester_sparsedense();
    fprintf('-------------------------------------------------------------\n\n');
    if ret, return; end

    fprintf('RUN EXAMPLE: example_glyap.m\n\n');
    ret = example_glyap();
    fprintf('-------------------------------------------------------------\n\n');
    if ret, return; end

    fprintf('RUN EXAMPLE: example_gstein.m\n\n');
    ret = example_gstein();
    fprintf('-------------------------------------------------------------\n\n');
    if ret, return; end

    fprintf('RUN EXAMPLE: example_sylvester_dense.m\n\n');
    ret = example_sylvester_dense();
    fprintf('-------------------------------------------------------------\n\n');
    if ret, return; end

    fprintf('RUN EXAMPLE: example_dense_nm_gmpare.m\n\n');
    ret = example_dense_nm_gmpare();
    fprintf('-------------------------------------------------------------\n\n');
    if ret, return; end

end
