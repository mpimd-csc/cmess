% MEX-M.E.S.S. (Matlab Interface to C-M.E.S.S.)
% This is MEX-M.E.S.S. the Matlab Interface to C-M.E.S.S.
% Version: @MESS_MAJOR_VERSION@.@MESS_MINOR_VERSION@.@MESS_PATCH_VERSION@
%
% Overview about available routines, solvers and classes:
%
% Enumeration Types
%   mess_direct_cholpackage_t           - available Cholesky based solvers
%   mess_direct_lupackage_t             - available LU based solver
%   mess_equation_t                     - available matrix equations types
%   mess_memusage_t                     - control memory usage
%   mess_multidirect_t                  - available Multi-LU based solvers
%   mess_operation_t                    - represents transposed and hermitian transposed operations
%   mess_parameter_t                    - available shift parameter heuristics
%   mess_residual_t                     - available residual computation methods
%   mess_norm_t                         - available norm types
%
% Classes
%   mess_handle                         - internal usage only
%   mess_options                        - options for matrix equations solvers
%   mess_options_adi                    - options for mess_lradi
%   mess_options_adi_shifts             - options for shift parameter heuristics
%   mess_options_newton                 - options for mess_lrnm
%   mess_status                         - contains status informations from matrix equations solvers
%
% Matrix Equations
%   mess_equation_glyap_dae1            - Lyapunov Equation arising from Index 1 DAE System
%   mess_equation_griccati_dae1         - Riccati Equation arising from Index 1 DAE System
%   mess_equation_glyap_dae2            - Lyapunov Equation arising from Index 2 DAE System
%   mess_equation_griccati_dae2         - Riccati Equation arising from Index 2 DAE System
%   mess_equation_glyap_so1             - Lyapunov Equation arising from Second Order System (Variant 1)
%   mess_equation_griccati_so1          - Riccati Equation arising from Second Order System (Variant 1)
%   mess_equation_glyap_so2             - Lyapunov Equation arising from Second Order System (Variant 2)
%   mess_equation_griccati_so2          - Riccati Equation arising from Second Order System (Variant 2)
%   mess_equation_glyap                 - standard/generalized Lyapunov Equation from Standard System
%   mess_equation_griccati              - standard/generalized Riccati Equation from Standard System
%   equation                            - Abstact Equation class, only for user defined callback equations
%   base_equation                       - Base Equation class, internal usage only
%
% Matrix Equations Solvers:
%   mess_care                           - Dense and Sparse Riccati Equation solver
%   mess_lyap                           - Dense and Sparse Lyapunov Equation solver
%   mess_lradi                          - Sparse Lyapunov Equation solver, Low-Rank ADI method
%   mess_lrnm                           - Sparse Riccati Equation solver, Low-Rank ADI-Newton method
%   mess_sylvester_sparsedense          - Sparse Dense Sylvester Equation solver
%   mess_sylvester_dense                - Dense Sylvester Equation solver
%   mess_glyap                          - Dense BLAS Level 3 Lyapunov equation solver
%   mess_gstein                         - Dense BLAS Level 3 Stein equation solver
%   mess_dense_nm_gmpare                - Dense Newton Based positive and negative Algebraic Riccati Equation solver
%
% Residual routines:
%   res2_lyap                           - 2-norm residual of standard/generalized Lyapunov equation
%   res2_ric                            - 2-norm residual of standard/generalized Riccati equation
%   res2_syl                            - 2-norm residual of sparse-dense or dense Sylvester equation
%   res2_stein                          - 2-norm residual of standard/generalized Stein equation
%   res2_gmpare                         - 2-norm residual of standard/generalized positive/negative Riccati equation
%
% Misc routines:
%   mess_call                           - gate to C-M.E.S.S. can call some function of C-M.E.S.S., in general internal use only
%   mess_version                        - print information about C-M.E.S.S.
%   mess_version_verbose                - print detailed information about C-M.E.S.S.
%   mess_version_major                  - returns the major version of C-M.E.S.S.
%   mess_version_minor                  - returns the minor version of C-M.E.S.S.
%   mess_version_patch                  - returns the patch version of C-M.E.S.S.
%   mess_git_id                         - returns the git id which was used to configure C-M.E.S.S. build
%   mess_git_branch                     - returns the git branch which was used to configure C-M.E.S.S. build
%   mess_direct_lu_select               - select the LU decomposition solver
%   mess_direct_chol_select             - select the Cholesky decomposition solver
%   mess_multidirect_select             - select the Multi-LU solver for shifted linear Systems (A+pE)
%
% Examples:
%   example_care                        - solve a genrealized Riccati equation, compute 2-norm residual
%   example_lyap                        - solve a generalized Lyapunov equation, compute 2-norm residual
%   example_lradi_rail                  - solve a generalized Lyapunov equation, compute 2-norm residual
%   example_lrnm_rail                   - solve a generalizd Riccati equation, compute 2-norm residual
%   example_sylvester_sparsedense       - solve a sparse dense Sylvester equation, compute 2-norm residual
%   example_sylvester_dense             - solve a dense Sylvester equation, compute 2-norm residual
%   example_glyap                       - solve a generalized Lyapunov equation using mess_glyap, compute 2-norm residual
%   example_gstein                      - solve a generalized Stein equation using mess_gstein, compute 2-norm residual
%   example_dense_nm_gmpare             - solve a standard / generalized positive Riccati equation using mess_dense_nm_gmpare, compute residual
%   run_examples                        - runs all defined examples
%
% Unit Tests:
%   run_tests                           - execute all Unit tests
%
%   Callback:
%       myequation                      - user defined callback equation
%       test_callback_lyapunov          - Lyapunov Equation solver (mess_lradi) test of user defined callback equation myequation
%       test_callback_riccati           - Riccati Equation solver (mess_lrnm) test of user defined callback equation myequation
%
%   Conversions:
%       test_enums                      - proper enum type conversion
%       test_matrix                     - conversion of matrix types
%       test_options                    - conversion of mess_options class to C-M.E.S.S.
%       test_status                     - conversion of mess_status class to C-M.E.S.S.
%       test_vector                     - conversion of vector types
%
%   Easy Interface:
%       test_care                       - mess_care tests
%       test_lyap                       - mess_lyap tests
%       test_dense_nm_gmpare            - mess_dense_nm_gmpare tests
%
%   Sylvester Solvers:
%       test_sylvester_sparsedense      - mess_sylvester_sparsedense tests
%       test_sylvester_dense            - mess_sylvester_dense tests
%
%   Glyap3 based solvers:
%       test_glyap                      - mess_glyap tests
%       test_gstein                     - mess_gstein tests
%
%   Low-Rank ADI:
%       test_lradi_filter               - mess_lradi tests for standard/generalized Lyapunov equations
%       test_lradi_nse_dae2             - mess_lradi tests for Lyapunov equation arising from Index 2 DAE system
%       test_lradi_triple_chain_so1     - mess_lradi tests for Lyapunov equation arising from Second Order system (variant 1)
%       test_lradi_triple_chain_so2     - mess_lradi tests for Lyapunov equation arising from Second Order system (variant 2)
%       test_lradi_ww_dae1              - mess_lradi tests for Lyapunov equation arising from Index 1 DAE system
%
%   Low-Rank Newton Method:
%       test_lrnm_filter                - mess_lrnm tests for standard/generalized Riccati equations
%       test_lrnm_nse_dae2              - mess_lrnm tests for Riccati equation arising from Index 2 DAE system
%       test_lrnm_triple_chain_so1      - mess_lrnm tests for Riccati equation arising from Second Order system (variant 1)
%       test_lrnm_triple_chain_so2      - mess_lrnm tests for Riccati equation arising from Second Order system (variant 2)
%       test_lrnm_ww_dae1               - mess_lrnm tests for Riccati equation arising from Index 1 DAE system
%
%   Miscellaneous:
%       test_mess_have                  - test of miscellaneous routines
%       test_res2_diffZ                 - test of norm computation routine
%
%
% Others:
%   mess_path                           - add MEX-M.E.S.S. to the Matlab path
%   mess_lint                           - Matlab lint test for MEX-M.E.S.S. code
%

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
