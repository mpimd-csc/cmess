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
 * @file include/mess/solvername.h
 * @brief Structure of solvernames.
 * @author @koehlerm
 */

#ifndef SOLVERNAME_H
#define SOLVERNAME_H

#define SET_SOLVERNAME(X,Y) do { \
    mess_try_alloc(X, char *, (strlen(Y)+1) * sizeof(char)); \
    strncpy(X, Y, strlen(Y)+1); \
} while(0)

#endif
