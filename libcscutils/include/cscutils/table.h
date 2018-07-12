/*
* CSCUTILS - A collection of various software routines uses in CSC projects
* Copyright (C) 2015 Martin Koehler
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

#ifndef CSC_TABLE_H
#define CSC_TABLE_H

#ifdef  __cplusplus
extern "C" {
#endif
    /**
     @file libcscutils/include/cscutils/table.h
      @defgroup table  Table: ASCII and Tex Tables

      This part of the library contains routines to deal with ASCII printable tables. The tables
      can be print as classical ASCII tables in output files as well a been formated to Latex code.

      @addtogroup table
      @{
    */

    /** Maximum character length of an entry in a table cell. */
    #define CSC_TABLE_MAXLEN 256


    /**
     @brief Enumeration for the columns data type.

     The csc_table_value_t enumeration represents the possible values in a column.
     */
    typedef enum {
        CSC_TABLE_INTEGER = 0,  /**< The column contains integer values. */
        CSC_TABLE_FLOAT = 1,    /**< The column contains double precision values. */
        CSC_TABLE_STRING = 2    /**< The column contains strings. */
    } csc_table_value_t;

    /**
     @brief Enumeration for the alignment in a column.

     The csc_table_align_t enumeration represents the possible alignments inside a column.
     */
    typedef enum {
        CSC_TABLE_CENTER = 0,  /**< The entry gets centered. */
        CSC_TABLE_LEFT = 1 ,   /**< The entry gets aligned to the left. */
        CSC_TABLE_RIGHT = 2    /**< The enrty gets aligned to the right. */
    } csc_table_align_t;

     /**
       @brief The csc_table_comment_t structure represents a comment for a table.

       The csc_table_comment_t structure represents a comment for the table. Additionally,
       it stores the comment sign which is printed in front of a comment. By default this is
       initialized with "#" to obtain an output which is easily parsed by \b sed or \b awk.
      */
    typedef struct _csc_table_comment_t {
        char start[CSC_TABLE_MAXLEN];       /**< The comment start string. Printed in front of every comment line. */
        char **lines;                       /**< Array of strings containing the comments. */
        int len;                            /**< Number of elements in lines. */
    } csc_table_comment_t;


    /**
      @brief The csc_table_formater_t typedef declares a function pointer to format cell entries.

      The csc_table_formater_t typedef declares a function pointer to declare cell entries. It takes
      an output buffer and its length as as input as well a the data type of the column as arguments. The
      last argument contains the actual value which needs to be formated. The function does not compute
      the alignment inside the table cell.
     */
    typedef void (*csc_table_formater_t)(char *, int, csc_table_value_t, ...);

    /**
      @brief The csc_table_comment_t structure represents a column inside a csc_table_t object.

      The csc_table_comment_t structure represents a column inside a csc_table_t structure. Each column
      can only contain data of a single data type.
     */
    typedef struct _csc_table_column_t {
        csc_table_value_t type;         /**< Data type of the cell entries. */
        char name[CSC_TABLE_MAXLEN];    /**< Name of the column. Used as headline while printing. */
        union {
            long   *integer_values;     /**< Array containing the integer values of the cell entries. */
            double *float_values;       /**< Array containing the double precision values of the cell entries. */
            char  **string_values;      /**< Array containing the strings in the cell entries. */
            void *ptr;                  /**< Array for future extension and generic access to the other entries. */
        } v;                            /**< Union representing the column entries. */
        char format_str[CSC_TABLE_MAXLEN]; /**< Printf compatible format string for the table cells. */
        csc_table_formater_t formater; /**< Format function for the table cells. If it is set than the format_str is ignored. */
        int *set;                      /**< Array which indicates which rows in the column are set.  */
        int len;                       /**< Number of elments in the column. */
        int width;                     /**< Width of the column. Recomputed everytime an entry is added. */
        csc_table_align_t align;       /**< Alignment of the column entries.  */
    } csc_table_column_t;

    /**
      @brief The csc_table_t structure represents a table.

      The csc_table_t structure represents a table which can be printed on the screen, to an ascii file,
      or exported as a Latex document. It is the main structure used by the csc_table module.
     */
    typedef struct _csc_table_t {
        int number_of_columns;             /**< Number of column in the table. */
        int number_of_rows;                /**< Number of rows in the table. */
        int current_row;                   /**< Current row. This is the index of the row all set entry function will use. */
        csc_table_column_t *columns;       /**< Array of the columns. */
        int cp;                            /**< Indicator if continous print is used. If true that last row is printed when a new row is started. */
        csc_table_comment_t * comment;     /**< Comments for the table. These are printed before the table starts.  */
    } csc_table_t;


    csc_table_t * csc_table_new(int continous_print);
    void csc_table_destroy(csc_table_t * t);
    void csc_table_print_ascii(FILE *stream, csc_table_t *t, const char *colsep);
    void csc_table_print_fortran(csc_table_t *t, const char *colsep);

    int csc_table_save_ascii(const char *filename, csc_table_t *t, const char *colsep);
    int csc_table_save_latex(const char *filename, csc_table_t *t, int standalone);

    csc_table_t * csc_table_new_from_table(csc_table_t *table);

    int  csc_table_add_column(csc_table_t *t, const char *name, csc_table_value_t type, csc_table_align_t align);
    void csc_table_column_destroy(csc_table_column_t col);
    int  csc_table_new_row(csc_table_t * t);
    int  csc_table_set_entry(csc_table_t *t, int column, ...);
    void  csc_table_set_entry_integer(csc_table_t *t, int column, int val);
    void  csc_table_set_entry_float(csc_table_t *t, int column, double val);
    void  csc_table_set_entry_string(csc_table_t *t, int column, char *val);

    int  csc_table_column_set_format(csc_table_t *t, int column, const char *fmt);
    int  csc_table_column_set_formater(csc_table_t *t, int column, csc_table_formater_t fmt);
    int csc_table_append_row(csc_table_t *t, csc_table_t *tab, int row);


    int csc_table_comment_printf(csc_table_t *t, const char * comment, ...);
    int csc_table_comment_text(csc_table_t *t, const char * text);
    int csc_table_comment_sign(csc_table_t *t, const char *sign);
    void csc_table_comment_date(csc_table_t *t );
    void csc_table_comment_cmd(csc_table_t *t, int argc, char ** argv);

    /* Table Operations   */
    int csc_table_max_row(csc_table_t * t, int column);
    int csc_table_min_row(csc_table_t * t, int column);


    /* Table comments internal management */
    csc_table_comment_t * csc_table_new_comment(void);
    void csc_table_destroy_comment(csc_table_comment_t * c);
    void csc_table_comment_start(csc_table_comment_t *c, const char *start);
    int  csc_table_comment_add(csc_table_comment_t *c, const char * fmt, ...);
    int  csc_table_comment_add_va(csc_table_comment_t *c, const char * fmt, va_list ap);
    void csc_table_comment_print(FILE *stream, csc_table_comment_t *c);
    void csc_table_comment_clear(csc_table_comment_t *c);



    /** @}
     *
     * @defgroup table_formaters  Different Formaters for Numeric Values
     * @ingroup table
     *
     * This group contains different formaters for the values inside the columns.
     *
     * @addtogroup table_formaters
     * @{
     **/

     void csc_table_formater_integer(char *out, int len, csc_table_value_t type, ...);
     void csc_table_formater_integer_latex(char *out, int len, csc_table_value_t type, ...);

 /** @} */
#ifdef  __cplusplus
}
#endif
#endif
