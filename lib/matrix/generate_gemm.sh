#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>.
#
# Copyright (C) Peter Benner, Martin Koehler, Jens Saak and others
#               2009-2018
#

cpp gemm_dense_sparse_v2.f90 -o gemm_dense_sparse_pp.f90
sed -i -e '/^#/d' -e '/^$/d' gemm_dense_sparse_pp.f90
cpp gemm_sparse_dense_v2.f90 -o gemm_sparse_dense_pp.f90
sed -i -e '/^#/d' -e '/^$/d' gemm_sparse_dense_pp.f90

