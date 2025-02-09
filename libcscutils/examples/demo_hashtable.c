/*
 * Demo Program - Hashtable from libcscutils
 * Copyright (C) Martin Koehler, 2015
 *
 * This library is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 2.1 of the License, or
 * (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this library; if not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include "cscutils/ds.h"

typedef struct _entry_t {
    char key[128];
    char value[128];
} entry_t;


char *getkey(csc_ds_object_t e2) {
    entry_t * e = (entry_t *) e2;
    return e->key;
}

void freeobj(csc_ds_object_t e2) {
    entry_t * e = (entry_t *) e2;
    free(e);
    e = NULL;
}

size_t hash (const char *key, size_t size)
{
    size_t h = 0;
    size_t i = 0;
    size_t len = strlen(key);

    for (i = 0; i < len; i++) {
        h = ( h + key[i]) % size;
    }
    return h;
}


entry_t * create_entry(const char *k, const char *v) {
    entry_t *e = calloc(1, sizeof(entry_t));
    strcpy ( e->key, k);
    strcpy ( e->value, v);
    return e;
}

int main(int argc, char **argv)
{
    csc_ds_t * ds;
    entry_t *e;



    ds = csc_ds_hashtable(17, getkey, hash , freeobj);
    if ( ds == NULL ) {
        fprintf(stderr, "Error Create Hashtable\n");
        return -1;
    }

    e = create_entry("Test", "Wert 1");
    csc_ds_insert(ds, e);
    e = create_entry("TestXZ", "Wert 2");
    csc_ds_insert(ds, e);
    e = create_entry("Hallo", "Welt!");
    csc_ds_insert(ds, e);

    printf("Initial Table\n");
    csc_ds_dump(stdout, ds);


    e = create_entry("Hello", "world");
    csc_ds_insert(ds, e);

    printf("After Insert Table\n");
    csc_ds_dump(stdout, ds);

    csc_ds_remove(ds, "TextXZ");

    printf("After Remove Table\n");
    csc_ds_dump(stdout, ds);

    e = csc_ds_find(ds, "Hallo");
    printf("%s - %s\n", e->key, e->value);


    csc_ds_free(ds);
    return 0;
}
