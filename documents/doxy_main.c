//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, see <http://www.gnu.org/licenses/>.
//
// Copyright (C) Peter Benner, Martin Koehler, Jens Saak and others
//               2009-2018
//

/**
  @mainpage M.E.S.S. - Matrix Equation Sparse Solver
  @section main_intro Welcome to M.E.S.S.

  @mess, the Matrix Equation Sparse Solver library, is the successor to the
  Lyapack Toolbox for @matlab.  It is available as a @matlab toolbox as well as
  a @c C library.  It is intended for the solution of symmetric linear anyd
  quadratic, differential and algebraic matrix equations with real, large, and
  sparse (or sparse + low-rank) coefficient matrices         in the linear
  terms, low-rank quadratic and constant terms, and linearow-rank solutions.
  The @c C version provides a large set of auxillary subroutines for sparse
  matrix computations and efficient usage of modern multicore workstations.


  @section main_features Features
  \li Solving large scale (generalized) Lyapunov and Riccati Equations.
  \li Model order reduction for 1st and 2nd order systems via Balanced Truncation.
  \li \f$ \mathcal{H}_2  \f$  model order reduction via the IRKA and TSIA algorithms.
  \li Tools for LQR problems.
  \li A uniform inferface to handle sparse and dense matrices.
  \li Shared memory parallel algorithms using @openmp and @pthreads.
  \li A solver for special sparse-dense Sylvester equations.
  \li A large set of direct and iterative solvers for linear systems.
  \li Various helper subroutines for thread pools, configuration files, interfaces to @python.

  @section main_overview Overview of the Functionality
  The libary is organized in different functional groups:
  \li \ref general "General Type Definitions and Macros"
  \li \ref matrix "Matrix Operations"
  \li \ref vector "Vector Operations"
  \li \ref direct "Direct Solvers for Linear Systems"
  \li \ref multisolvers "Solvers for Sets of Linear Systems"
  \li \ref itsolver "Iterative Solvers for Linear Systems"
  \li \ref dynsys "Dynamical Systems"
  \li \ref eigenvalues "Eigenvalue and Related Problems"
  \li \ref lrcfadi "Low-Rank Methods for Matrix Equations"
  \li \ref easyfrontend "Easy Frontends to Main Functions"
  \li \ref graph "Graph Algorithms"
  \li \ref matgen "Matrix Generators"
  \li \ref error   "Error Handling"
  \li \ref interfaces "Interfaces to other Programming Languages and Libraries"
  \li \ref plot "Plotting Interface"
  \li \ref misc "Miscellaneous Helper Functions"
  \li \ref tutorials "Tutorials"
  \li \ref test "Software Testing"

  @section supported_meq Supported Matrix Equations
  @subsection supported_meq_ale Algebraic Lyapunov equations
  \f[
     A X E^T + E X A^T + B B^T = 0
  \f]
  \f[
     A^T X E + E^T X A + C^T C = 0
  \f]
  These equations are solved via the low-rank (\f$X = Z Z^T\f$, or \f$X = L D
  L^T\f$) alternating directions implicit iteration (LR-ADI) using the
  function \ref mess_lradi. For the case where \f$E\f$, \f$A\f$ are real sparse matrices we
  also provide a lyap-like function \ref mess_lyap and \ref mess_glyap.

  @subsection supported_meq_are Continuous Time Algebraic Riccati equations
  \f[
     A X E^T + E X A^T - E X B B^T X E^T + C^T C = 0
  \f]
  \f[
    A^T X E + E^T X A - E^T X C^T C X E + B B^T = 0
  \f]
  For these standard control and filter equations we provide
  several solvers. Similar to \ref mess_lyap above we also provide a
  quick-start interface mess_care using a similar syntax as
  care. For experts and for use with the non-standard equations,
  we implemented a low-rank Newton ADI method with linesearch in
  \ref mess_lrnm.


  @section Installation
  For installation details please look at the \ref install "Installation" page.

  @section Contact
  The main website of @mess is located at  http://www.mpi-magdeburg.mpg.de/projects/mess.
  */

