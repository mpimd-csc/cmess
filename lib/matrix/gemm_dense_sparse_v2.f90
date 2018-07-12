!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program; if not, see <http:!www.gnu.org/licenses/>.
!
! Copyright (C) Peter Benner, Martin Koehler, Jens Saak and others
!               2009-2018
!

#define DD 0
#define DZ 1
#define ZD 2
#define ZZ 3

! C = A * B
#define TYPE DD
#define DATATYPE_A real*8
#define DATATYPE_B real*8
#define DATATYPE_C real*8
#include "gemm_dense_sparse_nn.f90"

#define TYPE DZ
#define DATATYPE_A real*8
#define DATATYPE_B complex*16
#define DATATYPE_C complex*16
#include "gemm_dense_sparse_nn.f90"

#define TYPE ZD
#define DATATYPE_A complex*16
#define DATATYPE_B real*8
#define DATATYPE_C complex*16
#include "gemm_dense_sparse_nn.f90"

#define TYPE ZZ
#define DATATYPE_A complex*16
#define DATATYPE_B complex*16
#define DATATYPE_C complex*16
#include "gemm_dense_sparse_nn.f90"

! C = A * B^T
#define TYPE DD
#define DATATYPE_A real*8
#define DATATYPE_B real*8
#define DATATYPE_C real*8
#include "gemm_dense_sparse_nt.f90"

#define TYPE DZ
#define DATATYPE_A real*8
#define DATATYPE_B complex*16
#define DATATYPE_C complex*16
#include "gemm_dense_sparse_nt.f90"

#define TYPE ZD
#define DATATYPE_A complex*16
#define DATATYPE_B real*8
#define DATATYPE_C complex*16
#include "gemm_dense_sparse_nt.f90"

#define TYPE ZZ
#define DATATYPE_A complex*16
#define DATATYPE_B complex*16
#define DATATYPE_C complex*16
#include "gemm_dense_sparse_nt.f90"

! C = A * B^H
#define TYPE DD
#define DATATYPE_A real*8
#define DATATYPE_B real*8
#define DATATYPE_C real*8
#include "gemm_dense_sparse_nh.f90"

#define TYPE DZ
#define DATATYPE_A real*8
#define DATATYPE_B complex*16
#define DATATYPE_C complex*16
#include "gemm_dense_sparse_nh.f90"

#define TYPE ZD
#define DATATYPE_A complex*16
#define DATATYPE_B real*8
#define DATATYPE_C complex*16
#include "gemm_dense_sparse_nh.f90"

#define TYPE ZZ
#define DATATYPE_A complex*16
#define DATATYPE_B complex*16
#define DATATYPE_C complex*16
#include "gemm_dense_sparse_nh.f90"



! C = A^T * B
#define TYPE DD
#define DATATYPE_A real*8
#define DATATYPE_B real*8
#define DATATYPE_C real*8
#include "gemm_dense_sparse_tn.f90"

#define TYPE DZ
#define DATATYPE_A real*8
#define DATATYPE_B complex*16
#define DATATYPE_C complex*16
#include "gemm_dense_sparse_tn.f90"

#define TYPE ZD
#define DATATYPE_A complex*16
#define DATATYPE_B real*8
#define DATATYPE_C complex*16
#include "gemm_dense_sparse_tn.f90"

#define TYPE ZZ
#define DATATYPE_A complex*16
#define DATATYPE_B complex*16
#define DATATYPE_C complex*16
#include "gemm_dense_sparse_tn.f90"

! C = A^T * B^T
#define TYPE DD
#define DATATYPE_A real*8
#define DATATYPE_B real*8
#define DATATYPE_C real*8
#include "gemm_dense_sparse_tt.f90"

#define TYPE DZ
#define DATATYPE_A real*8
#define DATATYPE_B complex*16
#define DATATYPE_C complex*16
#include "gemm_dense_sparse_tt.f90"

#define TYPE ZD
#define DATATYPE_A complex*16
#define DATATYPE_B real*8
#define DATATYPE_C complex*16
#include "gemm_dense_sparse_tt.f90"

#define TYPE ZZ
#define DATATYPE_A complex*16
#define DATATYPE_B complex*16
#define DATATYPE_C complex*16
#include "gemm_dense_sparse_tt.f90"

! C = A^T * B^H
#define TYPE DD
#define DATATYPE_A real*8
#define DATATYPE_B real*8
#define DATATYPE_C real*8
#include "gemm_dense_sparse_th.f90"

#define TYPE DZ
#define DATATYPE_A real*8
#define DATATYPE_B complex*16
#define DATATYPE_C complex*16
#include "gemm_dense_sparse_th.f90"

#define TYPE ZD
#define DATATYPE_A complex*16
#define DATATYPE_B real*8
#define DATATYPE_C complex*16
#include "gemm_dense_sparse_th.f90"

#define TYPE ZZ
#define DATATYPE_A complex*16
#define DATATYPE_B complex*16
#define DATATYPE_C complex*16
#include "gemm_dense_sparse_th.f90"

! C = A^H * B
#define TYPE DD
#define DATATYPE_A real*8
#define DATATYPE_B real*8
#define DATATYPE_C real*8
#include "gemm_dense_sparse_hn.f90"

#define TYPE DZ
#define DATATYPE_A real*8
#define DATATYPE_B complex*16
#define DATATYPE_C complex*16
#include "gemm_dense_sparse_hn.f90"

#define TYPE ZD
#define DATATYPE_A complex*16
#define DATATYPE_B real*8
#define DATATYPE_C complex*16
#include "gemm_dense_sparse_hn.f90"

#define TYPE ZZ
#define DATATYPE_A complex*16
#define DATATYPE_B complex*16
#define DATATYPE_C complex*16
#include "gemm_dense_sparse_hn.f90"

! C = A^H * B^T
#define TYPE DD
#define DATATYPE_A real*8
#define DATATYPE_B real*8
#define DATATYPE_C real*8
#include "gemm_dense_sparse_ht.f90"

#define TYPE DZ
#define DATATYPE_A real*8
#define DATATYPE_B complex*16
#define DATATYPE_C complex*16
#include "gemm_dense_sparse_ht.f90"

#define TYPE ZD
#define DATATYPE_A complex*16
#define DATATYPE_B real*8
#define DATATYPE_C complex*16
#include "gemm_dense_sparse_ht.f90"

#define TYPE ZZ
#define DATATYPE_A complex*16
#define DATATYPE_B complex*16
#define DATATYPE_C complex*16
#include "gemm_dense_sparse_ht.f90"

! C = A^H * B^H
#define TYPE DD
#define DATATYPE_A real*8
#define DATATYPE_B real*8
#define DATATYPE_C real*8
#include "gemm_dense_sparse_hh.f90"

#define TYPE DZ
#define DATATYPE_A real*8
#define DATATYPE_B complex*16
#define DATATYPE_C complex*16
#include "gemm_dense_sparse_hh.f90"

#define TYPE ZD
#define DATATYPE_A complex*16
#define DATATYPE_B real*8
#define DATATYPE_C complex*16
#include "gemm_dense_sparse_hh.f90"

#define TYPE ZZ
#define DATATYPE_A complex*16
#define DATATYPE_B complex*16
#define DATATYPE_C complex*16
#include "gemm_dense_sparse_hh.f90"

