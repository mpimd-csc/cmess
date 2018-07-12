%MESS_CITE print function to include references in the documentation of
% MEX-M.E.S.S.
%
%   MESS_CITE is a function which returns full reference information
%   for a given paper/book/preprint ect.
%   MESS_CITE is only used internal for providing references to papers
%
%   MESS_CITE() returns available keys.
%
%   MESS_CITE(KEY) return references for a KEY.
%   KEY must be a string. If KEY is not a dictory an error is returned.
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

function ref = mess_cite(key)

    %% create a container and add references
    refs = containers.Map();

    refs('Wei16')=sprintf([
        'Weichelt, Heiko K.','\n', ...
        'Numerical aspects of flow stabilization by Riccati feedback','\n',...
        'Dissertation, 2016'
    ]);

    refs('morGugSW13')=sprintf([
        'Gugercin, S. and Stykel, T. and Wyatt, S.','\n', ...
        'Model Reduction of Descriptor Systems by Interpolatory Projection Methods','\n', ...
        '2013, 10.1137/130906635'
    ]);

    refs('BenSSetal13')=sprintf([
        'Benner, P. and Saak, J. and Stoll, M. and Weichelt, H.','\n', ...
        'Efficient Solution of Large-Scale Saddle Point',' ',...
        'Systems Arising in Riccati-Based Boundary Feedback',' ',...
        'Stabilization of Incompressible Stokes Flow','\n', ...
        '2013, 10.1137/120881312'
    ]);

    refs('BaeBSetal15')=sprintf([
        'Baensch, E. and Benner, P. and Saak, J. and Weichelt, H.','\n', ...
        'Riccati-based boundary feedback stabilization of incompressible',' ', ...
        'Navier-Stokes flows','\n',...
        '2015, 10.1137/140980016'
    ]);

    refs('BenS10')=sprintf([
        'Benner, P. and Saak, J.','\n', ...
        'A Galerkin-Newton-ADI Method for Solving Large-Scale Algebraic Riccati Equations','\n', ...
        'Preprint, 2010', ', ',...
        '<a href="http://www.am.uni-erlangen.de/home/spp1253/wiki/index.php/Preprints">online available</a>'
    ]);

    refs('BenHSetal16')=sprintf([
        'Benner, Peter and Heinkenschloss, Matthias and Saak, Jens and Weichelt, Heiko','\n', ...
        'An inexact low-rank Newton-ADI method for large-scale algebraic Riccati equations','\n', ...
        '2016', ', ',...
        '<a href="http://www.sciencedirect.com/science/article/pii/S0168927416300824">online available</a>'
    ]);

    refs('Ham82a')=sprintf([
        'S. J. Hammarling','\n', ...
        'Newtons Method for Solving the Algebraic Riccati Equation','\n', ...
        '1982'
    ]);

    refs('BenKS12')=sprintf([
        'P. Benner and P. Kuerschner and J. Saak','\n', ...
        'Avoiding complex arithmetic in the low-rank ADI method efficiently','\n', ...
        '2012, 10.1002/pamm.201210308'
    ]);

    refs('BenKS13a')=sprintf([
        'Benner, P. and Kuerschner, P. and Saak, J.','\n', ...
        'A Reformulated Low-Rank ADI Iteration with Explicit Residual Factors','\n', ...
        '2013, 10.1002/pamm.201310273'
    ]);

    refs('BenKS13')=sprintf([
        'Benner, P. and Kuerschner, P. and Saak, J.','\n', ...
        'Efficient Handling of Complex Shift Parameters in the Low-Rank Cholesky Factor ADI Method','\n', ...
        '2013, 10.1007/s11075-012-9569-7'
    ]);

    refs('BenKS14b')=sprintf([
        'Benner, P. and Kuerschner, P. and Saak, J.','\n', ...
        'Self-Generating and Efficient Shift Parameters in',' ', ...
        'ADI Methods for Large Lyapunov and Sylvester Equations','\n', ...
        '2014'
    ]);

    refs('Kue16')=sprintf([
        'Patrick Kuerschner','\n', ...
        'Efficient Low-Rank Solution of Large-Scale Matrix Equations','\n', ...
        'ADI Methods for Large Lyapunov and Sylvester Equations','\n', ...
        '2016',', ',...
        '<a href="http://hdl.handle.net/11858/00-001M-0000-0029-CE18-2">online available</a>'
    ]);


    refs('Saa09')=sprintf([
        'J. Saak','\n', ...
        'Efficient Numerical Solution of Large Scale Algebraic Matrix Equations in',' ', ...
        'PDE Control and Model Order Reduction','\n', ...
        'Dissertation, 2009, ', ...
        '<a href="http://nbn-resolving.de/urn:nbn:de:bsz:ch1-200901642">online available</a>'
    ]);


    refs('morBenKS11')=sprintf([
        'Benner, P. and Koehler, M. and Saak, J.','\n', ...
        'Sparse-Dense Sylvester Equations in H_2-Model Order Reduction','\n', ...
        'Preprint, 2011'
    ]);


    refs('KoeS16a')=sprintf([
        'Koehler, M. and Saak, J.','\n',...
        'On BLAS Level-3 Implementations of Common Solvers for (Quasi-) Triangular Generalized Lyapunov Equations','\n', ...
        'ACM Transactions on Mathematical Software, 2016'
    ]);


    %% perform operation
    if nargin == 0
        ref = refs.keys();
        return;
    end

    %% check if key is valid and available
    assert(ischar(key),'key is no character array.');
    assert(refs.isKey(key),'key is not contained in the database.');

    %% get reference and return
    ref = refs(key);

end

