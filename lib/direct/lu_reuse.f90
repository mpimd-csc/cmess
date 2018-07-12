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

subroutine LURCSRC(rows, cols, values, colptr, rowptr, lcolptr, lrowptr, ucolptr, urowptr, havep, p, haveq, q, lvalues, uvalues, w)
 implicit none
 integer, intent(in) :: rows, cols
 integer, dimension(0:*), intent(in) :: colptr, rowptr
 integer, dimension(0:*), intent(in) :: lcolptr, lrowptr
 integer, dimension(0:*), intent(in) :: ucolptr, urowptr
 integer, dimension(0:*), intent(in) :: p,q
 integer, intent(in) :: havep, haveq
 complex*16,dimension(0:*), intent(in) :: values
 complex*16,dimension(0:*), intent(out) :: lvalues, uvalues
 complex*16,dimension(0:*), intent(inout) :: w

 integer :: i,j,k,col, row
 complex*16 :: alpha

 if ( havep ==0 .AND. haveq == 0) then
         row = 0
         do i = rowptr(row), rowptr(row+1)-1
                ! col = colptr(i)
                w(colptr(i)) = values(i)
         end do
         do i = urowptr(0), urowptr(1)-1
                uvalues(i) = w (ucolptr(i))
         end do
         do i = 0, cols-1
                w(i) = 0
         end do
         lvalues(0) = 1.0

         ! main loop
         do i = 1, rows-1
                ! copy row to working array
                row = i
                ! if ( havep == 1) row = p(row)
                do j = rowptr(row), rowptr(row+1)-1
                        ! col = colptr(j)
                        ! if ( haveq == 1) col = q(col)
                        w(colptr(j)) = values(j)
                end do
                !modify the row
                do j = lrowptr(i), lrowptr(i+1)-2
                        col = lcolptr(j)
                        alpha = w(col)/uvalues(urowptr(col))
                        do k = urowptr(col), urowptr(col+1)-1
                                w (ucolptr(k)) = w(ucolptr(k)) - alpha *uvalues(k)
                        end do
                        w(col) = 0
                        lvalues(j) = alpha
                end do
                lvalues(lrowptr(i+1)-1) = 1.0
                do j =  urowptr(i), urowptr(i+1)-1
                        uvalues(j) = w(ucolptr(j))
                        w(ucolptr(j)) = 0
                end do
         end do
 else
         if ( havep == 1) then
                row = p(0)
         else
                 row = 0
         end if
         do i = rowptr(row), rowptr(row+1)-1
                col = colptr(i)
                if ( haveq == 1) col = q(col)
                w(col) = values(i)
         end do
         do i = urowptr(0), urowptr(1)-1
                uvalues(i) = w (ucolptr(i))
         end do
         do i = 0, cols-1
                w(i) = 0
         end do
         lvalues(0) = 1.0

         ! main loop
         do i = 1, rows-1
                ! copy row to working array
                row = i
                if ( havep == 1) row = p(row)
                do j = rowptr(row), rowptr(row+1)-1
                        col = colptr(j)
                        if ( haveq == 1) col = q(col)
                        w(col) = values(j)
                end do
                !modify the row
                do j = lrowptr(i), lrowptr(i+1)-2
                        col = lcolptr(j)
                        alpha = w(col)/uvalues(urowptr(col))
                        do k = urowptr(col), urowptr(col+1)-1
                                w (ucolptr(k)) = w(ucolptr(k)) - alpha *uvalues(k)
                        end do
                        w(col) = 0
                        lvalues(j) = alpha
                end do
                lvalues(lrowptr(i+1)-1) = 1.0
                do j =  urowptr(i), urowptr(i+1)-1
                        uvalues(j) = w(ucolptr(j))
                        w(ucolptr(j)) = 0
                end do
         end do
 end if

end subroutine

subroutine LURCSRR(rows, cols, values, colptr, rowptr, lcolptr, lrowptr, ucolptr, urowptr, havep, p, haveq, q, lvalues, uvalues, w)
 implicit none
 integer, intent(in) :: rows, cols
 integer, dimension(0:*), intent(in) :: colptr, rowptr
 integer, dimension(0:*), intent(in) :: lcolptr, lrowptr
 integer, dimension(0:*), intent(in) :: ucolptr, urowptr
 integer, dimension(0:*), intent(in) :: p,q
 integer, intent(in) :: havep, haveq
 real*8,dimension(0:*), intent(in) :: values
 real*8,dimension(0:*), intent(out) :: lvalues, uvalues
 real*8,dimension(0:*), intent(inout) :: w

 integer :: i,j,k,col, row
 real*8  :: alpha

 if ( havep == 1) then
         row = p(0)
 else
         row = 0
 end if
 do i = rowptr(row), rowptr(row+1)-1
        col = colptr(i)
        if ( haveq == 1) col = q(col)
        w(col) = values(i)
 end do
 do i = urowptr(0), urowptr(1)-1
        uvalues(i) = w (ucolptr(i))
 end do
 do i = 0, cols-1
        w(i) = 0
 end do
 lvalues(0) = 1.0
! write (*,*) lrowptr(0:9)
 ! main loop
 do i = 1, rows-1
        ! copy row to working array
        row = i
        if ( havep == 1) row = p(row)
        do j = rowptr(row), rowptr(row+1)-1
                col = colptr(j)
                if ( haveq == 1) col = q(col)
                w(col) = values(j)
        end do
        !modify the row
        do j = lrowptr(i), lrowptr(i+1)-2
                col = lcolptr(j)
                alpha = w(col)/uvalues(urowptr(col))
                do k = urowptr(col), urowptr(col+1)-1
                        w (ucolptr(k)) = w(ucolptr(k)) - alpha *uvalues(k)
                end do
                w(col) = 0
                lvalues(j) = alpha
        end do
        ! if ( i == 1 ) then
                ! write (*,*) lrowptr(i+1), lrowptr(i+1)-1
        ! end if
        lvalues(lrowptr(i+1)-1) = 1.0
        do j =  urowptr(i), urowptr(i+1)-1
                uvalues(j) = w(ucolptr(j))
                w(ucolptr(j)) = 0
        end do
 end do

end subroutine


