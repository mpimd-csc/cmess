%RES2_DIFFZ  computes the 2-norm residual of a difference of two low-rank factors.
%
%   RES2_DIFFZ computes the 2-norm of the difference of 2 low-rank factors
%   Z1 and Z2, i.e. it computs || Z1*Z1' - Z2*Z2'||_2.
%
%   Arguments:
%       Z1  - real or complex, n-by-z1 matrix with z1<=n
%       Z2  - real or complex, n-by-z2 matrix with z2<=n
%
%   res2 = RES2_DIFFZ(Z1,Z2) comptes the difference of 2 low-rank factors
%   Z1 and Z2, i.e. it computes || Z1*Z1' -Z2*Z2'||_2.
%
%   Examples and Tests:
%       <a href="matlab:help test_res2_diffZ">test_res2_diffZ</a>, <a href="matlab:edit test_res2_diffZ">test_res2_diffZ.m</a>
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

function res2 = res2_diffZ(Z1,Z2)

    assert(ismatrix(Z1) && isnumeric(Z1),'Z1 is no numeric matrix.');
    assert(ismatrix(Z2) && isnumeric(Z1),'Z2 is no numeric matrix.');
    assert(size(Z1,1)>size(Z1,2),'Z1 is oversized.');
    assert(size(Z2,1)>size(Z2,2),'Z2 is oversized.');
    assert(size(Z1,1)==size(Z2,1),'Z1 and Z2 have different number of rows.');

    % compute norm(Z1*Z1' - Z2*Z2') and exploit that the matrix is hermitian,
    % it is the largest absolute value
    eopts.issym = true;
    eopts.isreal = isreal(Z1) && isreal(Z2);
    res2 = abs(eigs(@(x) Z1*(Z1'*x)-Z2*(Z2'*x),size(Z1,1),1,'LM',eopts));
end

