..
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   You should have received a copy of the GNU General Public License
   along with this program; if not, see <http://www.gnu.org/licenses/>.
   Copyright (C) Peter Benner, Martin Koehler, Jens Saak and others
                 2009-2018

Interfaced Methods
==================


Easy Interface
..............

lyap
----

.. currentmodule:: pymess.easy
.. autosummary::
    :nosignatures:

    lyap

.. autofunction:: lyap

care
----

.. currentmodule:: pymess.easy
.. autosummary::
    :nosignatures:

    care

.. autofunction:: care

sylvester_sparsedense
---------------------

.. currentmodule:: pymess.easy
.. autosummary::
    :nosignatures:

    sylvester_sparsedense

.. autofunction:: sylvester_sparsedense

positive_riccati
----------------

.. currentmodule:: pymess.dense_nm_gmpare
.. autosummary::
    :nosignatures:

    dense_nm_gmpare

.. autofunction:: dense_nm_gmpare


GLYAP3 Interface
................

glyap
-----

.. currentmodule:: pymess.glyap3
.. autosummary::
    :nosignatures:

    glyap

.. autofunction:: glyap

gstein
------

.. currentmodule:: pymess.glyap3
.. autosummary::
    :nosignatures:

    gstein

.. autofunction:: gstein


Low Rank Methods
.................

lradi
-----
.. currentmodule:: pymess.lradi
.. autosummary::
    :nosignatures:

    lradi

.. autofunction:: lradi


lrnm
----
.. currentmodule:: pymess.lrnm
.. autosummary::
    :nosignatures:

    lrnm

.. autofunction:: lrnm

residual
--------

.. currentmodule:: pymess.residual

.. autosummary::
    :nosignatures:

    res2_lyap
    res2_glyap
    res2_gstein
    res2_ric
    res2_sylv_sd
    res2_gmpare

.. autofunction:: res2_lyap
.. autofunction:: res2_glyap
.. autofunction:: res2_gstein
.. autofunction:: res2_ric
.. autofunction:: res2_sylv_sd
.. autofunction:: res2_gmpare


Other Methods
..............

.. currentmodule:: pymess.misc

.. autosummary::
    :nosignatures:

    eps
    mess_version
    mess_version_verbose
    mess_version_major
    mess_version_minor
    mess_version_patch
    mess_is_debug
    mess_have_bzip2
    mess_have_zlib
    mess_have_umfpack
    mess_have_amd
    mess_have_colamd
    mess_have_cholmod
    mess_have_csparse
    mess_have_superlu
    mess_have_mklpardiso
    mess_have_arpack
    mess_have_matio
    mess_have_openmp
    mess_have_mess64
    mess_git_id
    mess_git_branch
    set_errorlevel

.. autofunction:: eps
.. autofunction:: mess_version
.. autofunction:: mess_version_verbose
.. autofunction:: mess_version_major
.. autofunction:: mess_version_minor
.. autofunction:: mess_version_patch
.. autofunction:: mess_is_debug
.. autofunction:: mess_have_bzip2
.. autofunction:: mess_have_zlib
.. autofunction:: mess_have_umfpack
.. autofunction:: mess_have_amd
.. autofunction:: mess_have_colamd
.. autofunction:: mess_have_cholmod
.. autofunction:: mess_have_csparse
.. autofunction:: mess_have_superlu
.. autofunction:: mess_have_mklpardiso
.. autofunction:: mess_have_arpack
.. autofunction:: mess_have_matio
.. autofunction:: mess_have_openmp
.. autofunction:: mess_have_mess64
.. autofunction:: mess_git_id
.. autofunction:: mess_git_branch
.. autofunction:: set_errorlevel



Select (Multi)Direct Solvers
.............................

.. currentmodule:: pymess.direct_select

.. autosummary::
    :nosignatures:

    direct_select
    direct_chol_select
    multidirect_select

.. autofunction:: direct_select
.. autofunction:: direct_chol_select
.. autofunction:: multidirect_select





