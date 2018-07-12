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
 * @addtogroup test_io
 * @{
 * @file tests/io/check_fputs_block.c
 * @brief Check writing of large blocks.
 * @author @koehlerm
 * @test
 * This function checks if an
 * input filestream can be opened and if the reading and writing functions are working.
 * The tested functions belong to @cscutils package.
 *
 * @}
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <string.h>

#include "../call_macro.h"
#include "mess/mess.h"
#include "cscutils/io.h"


int main (int argc, char ** argv) {
    mess_init();
    char * buffer1, *buffer2;
    size_t len = CSC_IO_FILE_BUFFER_SIZE * 4 - 125 ;
    int ret =0;
    char *fn ;
    int i;
    csc_io_file_t *f;
    if ( argc != 2 )
        fn = "test.mtx";
    else
        fn = argv[1];

    buffer1 = (char *) malloc ( (len+1) *sizeof(char));
    buffer2 = (char *) malloc ( (len+1) *sizeof(char));
    buffer1[len]='\0';
    for ( i = 0; i < len; i++) {
        buffer1[i] = 'A' + i % 60;
    }

    f = csc_io_open(fn,CSC_IO_FILE_WRITE);
    if ( !f ) {
        fprintf(stderr,"Can not open FILE");
        free(buffer1);
        free(buffer2);
        return 1;
    } ;
    if ( csc_io_puts(buffer1, f) <= 0 ) {
        free(buffer1);
        free(buffer2);
        csc_io_close(f);
        return -1;
    }
    csc_io_close(f);

    f = csc_io_open(fn,CSC_IO_FILE_READ);
    if ( !f ) {
        fprintf(stderr,"Can not open FILE");
        free(buffer1);
        free(buffer2);
        csc_io_close(f);
        return 1; } ;
    if ( csc_io_gets(buffer2, len+1, f) == NULL) {
        free(buffer1);
        free(buffer2);
        csc_io_close(f);
        return -1;
    }
    csc_io_close(f);

    // printf("%s\n\n", buffer1);
    // printf("%s\n\n", buffer2);
    if ( strncmp(buffer1,buffer2,len) != 0 ) {
        fprintf(stderr,"ERROR.\n");
        ret = 1;
    } else {
        fprintf(stderr,"OK\n");
        ret = 0;
    }

    free(buffer1);
    free(buffer2);

    return ret;

}

