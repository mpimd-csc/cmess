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
 * @file lib/misc/misc.c
 * @brief Miscellaneous functions (print memory, round value).
 */


#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/errors.h"
#include "mess/matrix.h"
#include "mess/vector.h"
#include "mess/error_macro.h"
#include <complex.h>


/**
 * @brief Print @mess version.
 * @param[in] verbose input if nonzero the output is more verbose.
 *
 * The @ref  mess_versionx prints out the @mess version and the @ref mess_int_t data type.
 * If @p verbose is nonzero the output is more verbose.
 *
 * @sa mess_version
 * @sa mess_version_verbose
 *
 */
void mess_versionx(int verbose){
    MSG_PRINT("C-M.E.S.S. - The Matrix Equations Sparse Solver Library:\n" );

    MSG_PRINT("     Information about C-M.E.S.S. build:\n");
    MSG_PRINT("         - Version                   = %d.%d.%d\n",MESS_MAJOR_VERSION,MESS_MINOR_VERSION,MESS_PATCH_VERSION);
    MSG_PRINT("         - CMake Version             = %s\n",MESS_CMAKE_VERSION);
    MSG_PRINT("         - Build Config Time (UTC)   = %s\n",MESS_CONFIG_TIME);
    MSG_PRINT("         - git branch                = %s\n",MESS_GIT_BRANCH);
    MSG_PRINT("         - git commit hash           = %s\n",MESS_GIT_COMMIT_HASH);
    MSG_PRINT("         - Build Type                = %s\n",MESS_CMAKE_BUILD_TYPE);
    MSG_PRINT("         - C-Compiler                = %s\n",MESS_CMAKE_C_COMPILER);
    MSG_PRINT("         - CXX-Compiler              = %s\n",MESS_CMAKE_CXX_COMPILER);
    MSG_PRINT("         - Fortran-Compiler          = %s\n",MESS_CMAKE_Fortran_COMPILER);
    MSG_PRINT("         - Host System               = %s\n",MESS_CMAKE_HOST_SYSTEM);
    MSG_PRINT("         - BLAS Vendor               = %s\n",MESS_BLA_VENDOR);
    MSG_PRINT("         - BLAS Libaries             = %s\n",MESS_BLAS_LIBRARIES);
    MSG_PRINT("         - LAPACK Libaries           = %s\n",MESS_LAPACK_LIBRARIES);
    MSG_PRINT("         - sizeof(mess_int_t)        = %d",sizeof(mess_int_t));
    if (sizeof(mess_int_t) == 4 ) {
        MSG_PRINT(" (32 bit integer)\n");
    } else if(sizeof(mess_int_t) == 8) {
        MSG_PRINT(" (64 bit integer)\n");
    }else{
        MSG_PRINT("\n");
    }
#ifdef MESS_DEBUG
    MSG_PRINT("         - MESS_DEBUG                = YES\n");
#endif


    MSG_PRINT("\n");
    MSG_PRINT("     Available features:\n");
#ifdef MESS_HAVE_ZLIB
    MSG_PRINT("         - ZLIB\n");
#endif
#ifdef MESS_HAVE_BZ2
    MSG_PRINT("         - BZip2\n");
#endif
#ifdef MESS_HAVE_SUPERLU
    MSG_PRINT("         - SUPERLU\n");
#endif
#ifdef MESS_HAVE_UMFPACK
    MSG_PRINT("         - UMFPACK\n");
#endif
#ifdef MESS_HAVE_AMD
    MSG_PRINT("         - AMD\n");
#endif
#ifdef  MESS_HAVE_COLAMD
    MSG_PRINT("         - COLAMD\n");
#endif
#ifdef MESS_HAVE_CSPARSE
    MSG_PRINT("         - CSPARSE\n");
#endif
#ifdef MESS_HAVE_SUPERLU
    MSG_PRINT("         - SUPERLU\n");
#endif
#ifdef MESS_HAVE_PARDISO
    MSG_PRINT("         - PARDISO\n");
#endif
#ifdef MESS_HAVE_ARPACK
    MSG_PRINT("         - ARPACK\n");
#endif
#ifdef MESS_HAVE_MATIO
    MSG_PRINT("         - MATIO\n");
#endif
#ifdef MESS_HAVE_X11
    MSG_PRINT("         - X11\n");
#endif
#ifdef MESS_HAVE_XPM
    MSG_PRINT("         - X11_XPM\n");
#endif
#ifdef MESS_HAVE_OPENMP
    MSG_PRINT("         - OpenMP\n");
#endif
#ifdef MESS_HAVE_BACKTRACE
    MSG_PRINT("         - BACKTRACE\n");
#endif

#ifdef MESS_HAVE_MATIO
    MSG_PRINT("         - MATIO\n");
#endif
    if (verbose){
        MSG_PRINT("\n");
        MSG_PRINT("     Detailed Library and Compiler Flags Information:\n");
        MSG_PRINT("         - SUITESPARSE Libaries              = \n%s\n",MESS_SUITESPARSE_LIBRARIES);
        MSG_PRINT("         - SUITESPARSE Includes              = \n%s\n",MESS_SUITESPARSE_INCLUDE_DIR);
        MSG_PRINT("         - Arpack      Libaries              = %s\n",MESS_ARPACK_LIBRARIES);
        MSG_PRINT("         - SUPERLU     Includes              = %s\n",MESS_SUPERLU_INCLUDE_DIR);
        MSG_PRINT("         - SUPERLU     Libraries             = %s\n",MESS_SUPERLU_LIBRARIES);
        MSG_PRINT("         - Libraries                         = \n%s\n",MESS_LIBS);
        MSG_PRINT("         - Includes                          = \n%s\n",MESS_INCLUDE_DIR);
        MSG_PRINT("         - C-Compiler Flags                  = %s\n",MESS_CMAKE_C_FLAGS);
        MSG_PRINT("         - C-Compiler Flags Debug            = %s\n",MESS_CMAKE_C_FLAGS_DEBUG);
        MSG_PRINT("         - C-Compiler Flags Release          = %s\n",MESS_CMAKE_C_FLAGS_RELEASE);
        MSG_PRINT("         - CXX-Compiler Flags                = %s\n",MESS_CMAKE_CXX_FLAGS);
        MSG_PRINT("         - CXX-Compiler Flags Debug          = %s\n",MESS_CMAKE_CXX_FLAGS_DEBUG);
        MSG_PRINT("         - CXX-Compiler Flags Release        = %s\n",MESS_CMAKE_CXX_FLAGS_RELEASE);
        MSG_PRINT("         - Fortran-Compiler Flags            = %s\n",MESS_CMAKE_Fortran_FLAGS);
        MSG_PRINT("         - Fortran-Compiler Flags Debug      = %s\n",MESS_CMAKE_Fortran_FLAGS_DEBUG);
        MSG_PRINT("         - Fortran-Compiler Flags Release    = %s\n",MESS_CMAKE_Fortran_FLAGS_RELEASE);
    }
}


/**
 * @brief Print @mess version.
 *
 * The @ref  mess_version prints out the @mess version and the @ref mess_int_t data type.
 *
 */
void mess_version(){
    mess_versionx(0);
}

/**
 * @brief Print @mess version information.
 *
 * The @ref  mess_version_verbose prints defailed out the @mess version and the @ref mess_int_t data type.
 *
 */
void mess_version_verbose(){
    mess_versionx(1);
}


/**
 * @brief Return the git id which was used to configure the @mess build.
 * @return git id SHA-1 hash value of length 40.
 *
 * The @ref mess_git_id returns the git id (SHA-1 hash value of length 40) which was used to configure the
 * @mess build.
 *
 * @sa mess_git_branch
 *
 */
const char* mess_git_id(){
    return MESS_GIT_COMMIT_HASH;
}


/**
 * @brief Return the name of the git branch which was used to configure the @mess build.
 * @return name of git branch.
 *
 * The @ref mess_git_branch returns the name of the git branch which was used to configure the  @mess build.
 *
 */
const char* mess_git_branch(){
    return MESS_GIT_BRANCH;
}


/**
 * @brief Print a memory size.
 * @param[in] size input size to print
 *
 * The @ref mess_print_bytes prints out a memory size with the correct unit.
 *
 */
void mess_print_bytes(mess_int_t size){
    if ( size >= 1024*1024*1024){
        MSG_PRINT("%g Gbytes", (double)size/(1024.0*1024.0*1024.0));
        return;
    }
    if ( size >= 1024*1024){
        MSG_PRINT("%g Mbytes", (double)size/(1024.0*1024.0));
        return;
    }
    if ( size >= 1024){
        MSG_PRINT("%g kbytes", (double)size/(1024.0));
        return;
    }
    MSG_PRINT(MESS_PRINTF_INT " bytes", size);
    return ;
}


/**
 * @brief Round a value.
 * @param[in] x input value
 * @return integer value of \f$ x \f$
 *
 * The @ref nearbyint function rounds a value \f$ x \f$ to an integral value and returns it.
 *
 */
#ifndef MESS_HAVE_NEARBYINT
double nearbyint(double x){
    if ( x > 0 )
        return (mess_int_t)(x+0.5);
    else
        return (mess_int_t)(x-0.5);
}
#endif


int mess_version_major(void){return MESS_MAJOR_VERSION;}
int mess_version_minor(void){return MESS_MINOR_VERSION;}
int mess_version_patch(void){return MESS_PATCH_VERSION;}

int mess_is_debug(void){
#ifdef MESS_DEBUG
    return 1;
#endif
    return 0;
}

int mess_have_zlib(void){
#ifdef MESS_HAVE_ZLIB
    return 1;
#endif
    return 0;
}

int mess_have_bzip2(void){
#ifdef MESS_HAVE_BZIP2
    return 1;
#endif
    return 0;
}

int mess_have_umfpack(void){
#ifdef MESS_HAVE_UMFPACK
    return 1;
#endif
    return 0;
}

int mess_have_amd(void){
#ifdef MESS_HAVE_AMD
    return 1;
#endif
    return 0;
}

int mess_have_colamd(void){
#ifdef MESS_HAVE_COLAMD
    return 1;
#endif
    return 0;
}


int mess_have_cholmod(void){
#ifdef MESS_HAVE_CHOLMOD
    return 1;
#endif
    return 0;
}



int mess_have_csparse(void){
#ifdef MESS_HAVE_CSPARSE
    return 1;
#endif
    return 0;
}

int mess_have_superlu(void){
#ifdef MESS_HAVE_SUPERLU
    return 1;
#endif
    return 0;
}

int mess_have_mklpardiso(void){
#ifdef MESS_HAVE_MKLPARDISO
    return 1;
#endif
    return 0;
}

int mess_have_arpack(void){
#ifdef MESS_HAVE_ARPACK
    return 1;
#endif
    return 0;
}

int mess_have_matio(void){
#ifdef MESS_HAVE_MATIO
    return 1;
#endif
    return 0;
}

int mess_have_mess64(void){
#ifdef MESS64
    return 1;
#endif
    return 0;
}

int mess_have_openmp(void){
#ifdef MESS_HAVE_OPENMP
    return 1;
#endif
    return 0;
}


