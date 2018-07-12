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
 * @file include/mess/matrix.h
 * @brief Interface to handle matrices.
 * @author @koehlerm
 */

#ifndef matrix_H_
#define matrix_H_

#include "mess/vector.h"

#ifdef __cplusplus
extern "C" {
#endif
    /** @addtogroup matrix_ds
     * @{ */
    /**
     * @brief Enumeration to represent the storage format of a matrix.
     *
     * The @ref mess_storage_t enumeration is used to represent the storage type of a @ref mess_matrix object. \n
     * It is used for the store_type filed inside \ref mess_matrix_st.
     **/
    typedef enum {
        /** Matrix storage type is not set yet.  */
        MESS_UNKNOWN = 0,
        /** Matrix is stored in Compressed Sparse Row (@ref MESS_CSR) storage.  */
        MESS_CSR = 1,
        /** Matrix is stored in Compressed Sparse Column (@ref MESS_CSC) storage.  */
        MESS_CSC = 2,
        /** Matrix is stored in dense Fortran storage (@ref MESS_DENSE).  */
        MESS_DENSE = 3,
        /** Matrix is storage in sparse coordinate storage (@ref MESS_COORD) */
        MESS_COORD = 4
    } mess_storage_t;

    /**
     * @brief Macro to check if a matrix is stored in @ref MESS_CSR storage.
     *
     * The @ref MESS_IS_CSR checks if a matrix is stored in the @ref MESS_CSR storage.
     **/
#define MESS_IS_CSR(c)  ((*c).store_type == MESS_CSR)
    /**
     * @brief Macro to check if a matrix is stored in @ref MESS_CSC storage.
     *
     * The @ref MESS_IS_CSC checks if a matrix is stored in the @ref MESS_CSC storage.
     **/
#define MESS_IS_CSC(c)  ((*c).store_type == MESS_CSC)
    /**
     * @brief Macro to check if a matrix is stored in @ref MESS_DENSE storage.
     *
     * The @ref MESS_IS_DENSE checks if a matrix is stored in the @ref MESS_DENSE storage.
     **/
#define MESS_IS_DENSE(c)    ((*c).store_type == MESS_DENSE)
    /**
     * @brief Macro to check if a matrix is stored in @ref MESS_COORD storage.
     *
     * The @ref MESS_IS_COORD checks if a matrix is stored in the @ref MESS_COORD storage.
     **/
#define MESS_IS_COORD(c)    ((*c).store_type == MESS_COORD)
    /**
     * @brief Macro to check if a matrix is stored in some sparse  storage.
     *
     * The @ref MESS_IS_SPARSE checks if a matrix is stored in some sparse format. It
     * is the combination of \ref MESS_IS_CSR, \ref MESS_IS_CSC and \ref MESS_IS_COORD.
     **/
#define MESS_IS_SPARSE(c) (MESS_IS_CSR(c) || MESS_IS_CSC(c) || MESS_IS_COORD(c))


    /**
     * @brief Enumeration to represent the symmetry of a matrix.
     *
     * The @ref mess_symmetry_t enumeration is used to represent the symmetry of a matrix. It is used
     * for the symmetry flag inside the \ref mess_matrix_st structure. \n
     * At the moment only \ref MESS_GENERAL is used inside @mess. The other ones exist for further extensions.
     */
    typedef enum {
        /** The matrix is a general one without any symmetry.  */
        MESS_GENERAL = 0,
        /** The matrix is symmetric. */
        MESS_SYMMETRIC = 1,
        /** The matrix is skew-symmetric */
        MESS_SKEWSYMMETRIC = 2,
        /** The matrix is hermitian symmetric.*/
        MESS_HERMITIAN = 3
    } mess_symmetry_t;

    /**
     * @brief Check if a matrix is a general one, that means it has no special symmetry.
     *
     * The @ref MESS_IS_GENERAL macro checks if the matrix has no special symmetric structure.
     * */
#define MESS_IS_GENERAL(m)  ((*m).symmetry == MESS_GENERAL)
    /**
     * @brief Check if a matrix is a general one.
     *
     * The @ref MESS_IS_UNSYMMETRIC macro is used as a wrapper for \ref MESS_IS_GENERAL.
     * */
#define MESS_IS_UNSYMMETRIC(m)  ((*m).symmetry == MESS_GENERAL)
    /**
     * @brief Check if a matrix is stored symmetric.
     *
     * The @ref MESS_IS_SYMMETRIC macro checks if a matrix is stored symmetric. That means
     * the algorithms have to take care about implicit symmetry.\n
     * At the moment this is not used in @mess.
     */
#define MESS_IS_SYMMETRIC(m)    ((*m).symmetry == MESS_SYMMETRIC)
    /**
     * @brief Check if a matrix is stored skew-symmetric.
     *
     * The @ref MESS_IS_SKEWSYMMETRIC macro checks if a matrix is stored skew symmetric. That means
     * the algorithms have to take care about implicit skew symmetry.\n
     * At the moment this is not used in @mess.
     */
#define MESS_IS_SKEWSYMMETRIC(m)    ((*m).symmetry == MESS_SKEWSYMMETRIC)
    /**
     * @brief Check if a matrix is stored symmetric.
     *
     * The @ref MESS_IS_HERMITIAN macro checks if a matrix is stored hermitian symmetric. That means
     * the algorithms have to take care about implicit hermitian symmetry. \n
     * At the moment this is not used in @mess.
     */
#define MESS_IS_HERMITIAN(m)    ((*m).symmetry == MESS_HERMITIAN)


    /**
     * @brief The mess_matrix structure to handle dense and sparse matrices in a common context.
     *
     * This structure represents a matrix which is stored as compressed row, compressed column or dense storage.
     * The used storage stored matrix is determined by the store_type property. In case of dense storage the matrix is stored
     * in the Fortran format. Otherwise a 0 based indirect addressing format like @ref MESS_CSR, @ref MESS_CSC or @ref MESS_COORD is used.
     *
     * The data_type flag indicates the data type of the matrix. It is set to one of the values of the \ref mess_datatype_t
     * enumeration. If the matrix is real the values array is used to store the data. In case of a complex matrix
     * the values_cpx array is used. The unused one of those two arrays should be set to @c NULL.
     *
     * @sa mess_matrix_convert
     *
     * \attention
     * The symmetry flag is supposed to indicate whether the matrix is symmetric or not. Unfortunately we
     * do not use this at the moment. This is the reason why the symmetry is always set to \ref MESS_GENERAL
     * even if you read a symmetric matrix from a file.
     */
    typedef struct mess_matrix_st {
        /** Number of rows of the matrix.  */
        mess_int_t rows;
        /** Number of columns of the matrix.  */
        mess_int_t cols;
        /** Leading dimension of the matrix. This is only used for dense matrices.  */
        mess_int_t ld;
        /** Number of non zero elements in the matrix. This is only used for sparse matrices. */
        mess_int_t nnz;
        /** Array to store the row pointer or the row indices. @c NULL for dense matrices. */
        mess_int_t *rowptr;
        /** Array to store the column pointer or column indices. @c NULL for dense matrices.  */
        mess_int_t *colptr;
        /** Array containing the real matrix entries.  */
        double *values;
        /** Array containing the complex matrix entries.  */
        mess_double_cpx_t *values_cpx;
        /** Identify the matrix storage type. This can be one of the \ref mess_storage_t enumeration. */
        mess_storage_t store_type;
        /** Symmetry flag. This flag is always set to \ref MESS_GENERAL at the moment. Primarily used to support implicit symmetry. @sa mess_symmetry_t */
        mess_symmetry_t symmetry;
        /** Determine the data type of the matrix. This flag identifies whether values or values_cpx is used.  */
        mess_datatype_t data_type;
    } mess_matrix_st;

    /**
     * @brief Type definition to make \ref mess_matrix_st easier to use.
     *
     * This type definition creates a pointer out of the \ref mess_matrix_st structure. \n
     * This pointer is the data type which should be used for the matrices. Most functions take a mess_matrix
     * object directly. Only \ref mess_matrix_init and \ref mess_matrix_clear need a pointer to a mess_matrix,
     * that means a double pointer to a \ref mess_matrix_st structure.
     *
     * @sa mess_matrix_init
     * @sa mess_matrix_clear
     * */
    typedef struct mess_matrix_st*  mess_matrix;

    /**
     * @brief Type definition and structure for mess_matrix files.
     *
     * The @ref mess_matrix_file type definition and structure defines a @ref mess_matrix file.
     * */
    typedef struct mess_matrix_file_st {
        int8_t  magic[6];           /**< Mess_matrix binary file */
        int8_t  version;            /**< Binary format version */
        int8_t  sizemess_int_t;     /**< Size of mess_int_t in file */
        int64_t cols;               /**< Number of columns */
        int64_t rows;               /**< Number of rows */
        int64_t nnz;                /**< Number of non zero elements */
        int64_t ld;                 /**< Number of leading dimension */
        int8_t  data_type;          /**< Data type */
        int8_t  store_type;         /**< Storage type */
        int8_t  symmetry;           /**< Symmetry */
    } mess_matrix_file;

    /** Definition of a @ref mess_matrix binary file. */
#define MESS_MATRIX_FILE_MAGIC {0x63,0x6D,0x65,0x73,0x73,0x4D}
    /**  Macro to set a file to a valid @ref mess_matrix binary file.*/
#define MESS_MATRIX_FILE_MAGIC_SET(m) {m[0]=0x63; m[1]=0x6D; m[2]=0x65; m[3]=0x73 ; m[4]=0x73 ; m[5]=0x4D; }
    /**  Macro to check if a file is a valid @ref mess_matrix binary file. */
#define MESS_MATRIX_FILE_MAGIC_CHECK(m) ( (m[0] == 0x63 )&&(m[1]==0x6D) &&(m[2]==0x65)&&(m[3]==0x73)&&(m[4]==0x73)&&(m[5]==0x4D))
    /** Definition of the binary format version */
#define MESS_MATRIX_FILE_VERSION 0x01


    /**
     * @brief Check if a matrix is stored in the correct storage format and converts it if necessary.
     * @param[in]  input    Matrix object which storage formart should be checked
     * @param[out] output   Matrix object with the desired storage format
     * @param[out] status   Status of the conversion
     * @param[in]  fmt      input Desired storage format
     *
     * The @ref MESS_MATRIX_CHECKFORMAT macro checks if a given input matrix has the desired storage format fmt. If
     * the matrix is already in the correct format the output object will be set to the input object and stat will
     * be -1. If the input is not of the desired format output will be a newly created matrix object and contains a
     * copy of the input matrix in the desired storage format. If the conversion was successful  stat contains 0 otherwise
     * and positive error code. If stat is 0 the matrix output needs to be freed using \ref mess_matrix_clear.
     *
     * For example if one want to guarantee that the matrix \c A is in @ref MESS_CSR storage:
     * \code{.c}
     *  int status;
     *  mess_matrix helper;
     *  MESS_MATRIX_CHECKFORMAT(A, helper, status, MESS_CSR);
     *  ... Work with helper ...
     *  if ( status == 0 ) {
     *  mess_matrix_clear(&helper);
     *  }
     *  \endcode
     *
     */
#define MESS_MATRIX_CHECKFORMAT(input, output, status, fmt){                                                                            \
    if ( (input)->store_type == (fmt)){                                                                                                 \
        (output) = (input);                                                                                                             \
        (status) = -1;                                                                                                                  \
    } else {                                                                                                                            \
        mess_matrix_init(&output);                                                                                                      \
        MSG_WARN("convert matrix \"%s\" from %s to %s\n", #input, mess_storage_t_str((input)->store_type), mess_storage_t_str(fmt));    \
        (status) = mess_matrix_convert((input),(output), (fmt));                                                                        \
        FUNCTION_FAILURE_HANDLE((status), ((status)>0), mess_matrix_convert);                                                           \
    }                                                                                                                                   \
}

    /**
     * @brief Resets a mess_matrix object.
     * @param[in,out] matrix  Matrix object to reset
     *
     * The MESS_MATRIX_RESET is a wrapper around \ref mess_matrix_reset .
     */
#define MESS_MATRIX_RESET(matrix)  do { mess_matrix_reset((matrix)); } while(0)

    int mess_matrix_clear(mess_matrix *matrix);
    int mess_matrix_init(mess_matrix *matrix);
    /** @} */

    /** @addtogroup matrix_io
     * @{ */
    int mess_matrix_read(const char *filename, mess_matrix matrix);
    int mess_matrix_write(const char * filename, mess_matrix matrix);
    int mess_matrix_read_formated(const char *filename, mess_matrix matrix, mess_storage_t format );

    int mess_matrix_display(mess_matrix matrix);
    int mess_matrix_print(mess_matrix matrix);
    int mess_matrix_printshort(mess_matrix matrix);
    int mess_matrix_printdata(mess_matrix mat);
    int mess_matrix_printinfo(mess_matrix mat);

    /** @} */

    /** @addtogroup matrix_data
     * @{ */
    int mess_matrix_alloc(mess_matrix matrix, mess_int_t rows, mess_int_t cols, mess_int_t nnz, mess_storage_t storetype, mess_datatype_t datatype);
    void mess_matrix_alloc_setld(mess_int_t ld_m);
    int mess_matrix_realloc(mess_matrix matrix, mess_int_t rows, mess_int_t cols, mess_int_t nnz);
    int mess_matrix_resize(mess_matrix matrix, mess_int_t rows, mess_int_t cols);
    int mess_matrix_reset(mess_matrix matrix);

    int mess_matrix_convert(mess_matrix input, mess_matrix output, mess_storage_t outtype);
    int mess_matrix_sort(mess_matrix mat);
    int mess_matrix_cat(mess_matrix A, mess_matrix B, mess_matrix C, mess_matrix D, mess_storage_t otype, mess_matrix O);
    int mess_matrix_rowsub(mess_matrix input, mess_int_t srow, mess_int_t erow, mess_matrix out);
    int mess_matrix_colsub(mess_matrix input, mess_int_t scol, mess_int_t ecol, mess_matrix out);
    int mess_matrix_sub(mess_matrix in, mess_int_t rowS, mess_int_t rowE, mess_int_t colS, mess_int_t colE, mess_matrix out);
    int mess_matrix_lift(mess_matrix in, mess_int_t n, mess_matrix out);

    int mess_matrix_realpart(mess_matrix matrix, mess_matrix realpart);
    int mess_matrix_imagpart(mess_matrix matrix, mess_matrix imagpart);
    int mess_matrix_complex_from_parts(mess_matrix xr, mess_matrix xc, mess_matrix x);
    int mess_matrix_getelement(mess_matrix matrix, mess_int_t row, mess_int_t col, double *val, mess_double_cpx_t *valc);
    int mess_matrix_setelement_complex(mess_matrix matrix, mess_int_t row, mess_int_t col, mess_double_cpx_t value);
    int mess_matrix_setelement(mess_matrix matrix, mess_int_t row, mess_int_t col, double value);
    int mess_matrix_convert_csr_csc(mess_matrix inmatrix, mess_matrix outmatrix);
    int mess_matrix_convert_csc_csr(mess_matrix inmatrix, mess_matrix outmatrix);
    int mess_matrix_symfillup(mess_matrix mat) ;
    int mess_matrix_dense_from_farray(mess_matrix mat, mess_int_t rows, mess_int_t cols, mess_int_t ld, double *realv, mess_double_cpx_t *complexv);
    int mess_matrix_dense_from_carray(mess_matrix mat, mess_int_t rows, mess_int_t cols, double **realv, mess_double_cpx_t **complexv);
    int mess_matrix_csr(mess_matrix matrix, mess_int_t rows, mess_int_t cols, mess_int_t * rowptr, mess_int_t * colptr, double * values, mess_double_cpx_t *values_cpx);
    int mess_matrix_csc(mess_matrix matrix, mess_int_t rows, mess_int_t cols, mess_int_t * rowptr, mess_int_t * colptr, double * values, mess_double_cpx_t *values_cpx);
    int mess_matrix_coord(mess_matrix matrix, mess_int_t rows, mess_int_t cols, mess_int_t nnz, mess_int_t onebased, mess_int_t * rowptr, mess_int_t * colptr, double * values, mess_double_cpx_t *values_cpx);

    int mess_matrix_eye(mess_matrix A, mess_int_t rows, mess_int_t cols, mess_storage_t storetype);
    int mess_matrix_eyec(mess_matrix mat, mess_int_t rows, mess_int_t cols, mess_storage_t store);
    int mess_matrix_tridiag(mess_matrix matrix, mess_int_t rows, mess_int_t cols, mess_storage_t store_t, mess_datatype_t data_t, mess_double_cpx_t lower, mess_double_cpx_t diag, mess_double_cpx_t upper);
    int mess_matrix_rand_init(mess_int_t *seed);
    int mess_matrix_rand(mess_matrix mat, mess_int_t rows, mess_int_t cols, mess_storage_t storetype, mess_datatype_t dt, double p);
    int mess_matrix_rand_dense(mess_matrix mat, mess_int_t rows, mess_int_t cols, mess_datatype_t dt );
    int mess_matrix_rand_dense_uniform( mess_matrix mat, mess_int_t *seed , mess_int_t rows, mess_int_t cols);
    int mess_matrix_rand_csr(mess_matrix mat, mess_int_t rows, mess_int_t cols , double p, mess_datatype_t dt );
    int mess_matrix_rand_csc(mess_matrix mat, mess_int_t rows, mess_int_t cols , double p, mess_datatype_t dt );
    int mess_matrix_rand_coord(mess_matrix mat, mess_int_t rows, mess_int_t cols , double p, mess_datatype_t dt );
    int mess_matrix_dupl(mess_matrix mat);
    int mess_matrix_eliminate_zeros(mess_matrix mat);
    int mess_matrix_zeros(mess_matrix matrix);
    int mess_matrix_ones(mess_matrix matrix);
    int mess_matrix_one_value(mess_matrix matrix, double value);
    int mess_matrix_one_valuec(mess_matrix matrix, mess_double_cpx_t value);
    int mess_matrix_need_alloc(mess_matrix matrix, mess_int_t rows,mess_int_t cols,mess_int_t nnz, mess_storage_t storetype,mess_datatype_t datatype);

    int mess_matrix_copy(mess_matrix in, mess_matrix out);
    int mess_matrix_getcol(mess_matrix matrix, mess_int_t col, mess_vector c);
    int mess_matrix_getrow(mess_matrix matrix, mess_int_t row, mess_vector r);
    int mess_matrix_setcol(mess_matrix matrix, mess_int_t col, mess_vector colv);
    int mess_matrix_triu(mess_matrix mat, mess_int_t k);
    int mess_matrix_tril(mess_matrix mat, mess_int_t k);
    int mess_matrix_rowscalem(mess_vector r, mess_matrix mat);
    int mess_matrix_colscalem(mess_vector c, mess_matrix mat);
    int mess_matrix_setrow(mess_matrix matrix, mess_int_t row, mess_vector rowv);
    int mess_matrix_diag(mess_matrix matrix, mess_vector d);
    int mess_matrix_diag_from_vector(mess_vector v, mess_matrix diag);
    int mess_matrix_diagpos(mess_matrix matrix, mess_int_t *pos);
    int mess_matrix_addcols(mess_matrix matrix, mess_matrix toadd);
    int mess_matrix_addcols2p(mess_matrix Z , double s1, mess_matrix V1, double s2, mess_matrix V2);
    int mess_matrix_addcols1(mess_matrix Z , double s1, mess_matrix V1);
    int mess_matrix_toreal(mess_matrix m);
    int mess_matrix_tocomplex(mess_matrix m);
    int mess_matrix_totype(mess_matrix mat, mess_datatype_t dt);
    int mess_matrix_view(mess_matrix original, mess_int_t srow, mess_int_t scol, mess_int_t rows, mess_int_t cols, mess_matrix view);
    int mess_matrix_bandwidth(mess_matrix matrix, mess_int_t *kl, mess_int_t *ku);
    /** @}  */


    /*-----------------------------------------------------------------------------
     *  Norms
     *-----------------------------------------------------------------------------*/
    /** @addtogroup matrix_norm
     * @{ */


    int mess_matrix_norm1(mess_matrix A , double *nrm);
    int mess_matrix_norm2(mess_matrix A , double *nrm);
    int mess_matrix_norm2inv(mess_matrix A, double *nrm);
    int mess_matrix_normf(mess_matrix A , double *nrm);
    int mess_matrix_normf2(mess_matrix A, double *nrm );
    int mess_matrix_norminf(mess_matrix A, double *nrm );
    int mess_matrix_norm(mess_matrix A , mess_norm_t nrm_t, double *nrm);
    int mess_matrix_diffnorm(mess_matrix A, mess_matrix B, double *nrm );
    int mess_matrix_diffnormf(mess_matrix A, mess_matrix B, double *nrm );
    int mess_matrix_dynorm2(mess_matrix A, double *nrm );
    int mess_matrix_dynormf(mess_matrix A, double *nrm );
    int mess_matrix_dynormf2(mess_matrix A, double *nrm );
    int mess_matrix_indefinite_dynorm2(mess_matrix W, mess_matrix K, double *nrm);
    int mess_matrix_indefinite_dynormf(mess_matrix W, mess_matrix K, double *nrm);
    int mess_matrix_indefinite_dynormf2(mess_matrix W, mess_matrix K, double *nrm);
    int mess_matrix_mulnormf2(mess_operation_t op1, mess_matrix A, mess_operation_t op2, mess_matrix B, double *nrm);

    /** @} */

    /** @addtogroup matrix_op
     * @{ */
    int mess_matrix_ctranspose (mess_matrix A, mess_matrix AH);
    int mess_matrix_transpose (mess_matrix A, mess_matrix AT);
    int mess_matrix_xtranspose (mess_matrix A, mess_matrix AT, int hermitian);
    int mess_matrix_conj (mess_matrix A);

    int mess_matrix_proj_sym(mess_matrix matrix);
    /*-----------------------------------------------------------------------------
     *  rank and orthogonal subspaces
     *-----------------------------------------------------------------------------*/
    int mess_matrix_rank(mess_matrix A, mess_int_t *rank);
    int mess_matrix_orth(mess_matrix A, mess_matrix Q);
    int mess_matrix_null(mess_matrix A, mess_matrix Z);
    int mess_matrix_mgs ( mess_matrix A, mess_matrix Q, mess_matrix R );
    int mess_matrix_mgs_inplace ( mess_matrix Q );
    int mess_matrix_mgs_add ( mess_matrix Z, mess_matrix V, mess_matrix E );

    int mess_matrix_mvp(mess_operation_t op, mess_matrix A, mess_vector x, mess_vector y);
    int mess_matrix_gaxpy( mess_operation_t op, mess_matrix A, mess_vector x , mess_vector y );

    int mess_matrix_colsums(mess_matrix A, mess_vector x);
    int mess_matrix_rowsums(mess_matrix A, mess_vector x);

    int mess_matrix_add(double alpha, mess_matrix A, double beta, mess_matrix B);
    int mess_matrix_addc(mess_double_cpx_t alpha, mess_matrix A, mess_double_cpx_t beta, mess_matrix B);

    int mess_matrix_multiply(mess_operation_t opA, mess_matrix A, mess_operation_t opB, mess_matrix B, mess_matrix C);

    int mess_matrix_scale(double alpha, mess_matrix mat);
    int mess_matrix_scalec(mess_double_cpx_t alpha, mess_matrix mat);

    int mess_matrix_fro_inn(mess_operation_t op, mess_matrix A, mess_matrix B, void *fro);

    int mess_matrix_trace(mess_matrix A, double *tr);
    int mess_matrix_tracec(mess_matrix A, mess_double_cpx_t *tr);

    int mess_matrix_biorth(mess_matrix Vin, mess_matrix Win, mess_matrix Vout, mess_matrix Wout);
    int mess_matrix_gbiorth(mess_matrix E, mess_matrix Vin, mess_matrix Win, mess_matrix Vout, mess_matrix Wout);
    int mess_matrix_proj_mgs(mess_vector in, mess_matrix Q, mess_matrix QT, mess_int_t colcount, mess_vector out ,double *nrmout);
    int mess_matrix_proj_gmgs(mess_vector in, mess_matrix Q, mess_matrix QT, mess_matrix E, char *Eflag, mess_int_t colcount, mess_vector out, double *nrmout);

    int mess_matrix_map(mess_matrix mat, double (*f_real)(double), mess_double_cpx_t (*f_cpx)(mess_double_cpx_t), int full);
    int mess_matrix_map_abs(mess_matrix mat, int full);
    int mess_matrix_map_acos(mess_matrix mat, int full);
    int mess_matrix_map_acosh(mess_matrix mat, int full);
    int mess_matrix_map_arg(mess_matrix mat, int full);
    int mess_matrix_map_asin(mess_matrix mat, int full);
    int mess_matrix_map_asinh(mess_matrix mat, int full);
    int mess_matrix_map_atan(mess_matrix mat, int full);
    int mess_matrix_map_atanh(mess_matrix mat, int full);
    int mess_matrix_map_ceil(mess_matrix mat, int full);
    int mess_matrix_map_conj(mess_matrix mat, int full);
    int mess_matrix_map_cos(mess_matrix mat, int full);
    int mess_matrix_map_cosh(mess_matrix mat, int full);
    int mess_matrix_map_exp(mess_matrix mat, int full);
    int mess_matrix_map_expm1(mess_matrix mat, int full);
    int mess_matrix_map_floor(mess_matrix mat, int full);
    int mess_matrix_map_log(mess_matrix mat, int full);
    int mess_matrix_map_not(mess_matrix mat, int full);
    int mess_matrix_map_round(mess_matrix mat, int full);
    int mess_matrix_map_sin(mess_matrix mat, int full);
    int mess_matrix_map_sinh(mess_matrix mat, int full);
    int mess_matrix_map_sqrt(mess_matrix mat, int full);
    int mess_matrix_map_tan(mess_matrix mat, int full);
    int mess_matrix_map_tanh(mess_matrix mat, int full);
    int mess_matrix_map_isfinite(mess_matrix mat, int full);
    int mess_matrix_map_isinf(mess_matrix mat, int full);
    int mess_matrix_map_isnan(mess_matrix mat, int full);

    int mess_matrix_any(mess_matrix mat, mess_int_t (*f_real) (double), mess_int_t (*f_cpx) (mess_double_cpx_t), mess_int_t dim, mess_vector anyvec);


    /** @} */

    /** @addtogroup matrix_misc
     * @{ */

    /**
     * @brief Enumeration to collect the different matrix reordering.
     *
     * The @ref mess_reorder_t enumeration is used to select different matrix
     * reorderings. Mostly they are used to reduce the fill-in in direct
     * sparse solvers or to increase the data-locality in matrix-vector
     * products.
     *
     * @sa mess_matrix_reorder
     * @sa mess_matrix_reorder_amd
     * @sa mess_matrix_reorder_colamd
     * @sa mess_matrix_reorder_rcm
     * */
    typedef enum {
        /** Select the natural (identity) reordering */
        MESS_REORDER_NONE = 0,
        /** Select the Reverse Cuthill McKee (RCM) reordering. */
        MESS_REORDER_RCM  = 1,
        /** Select the AMD reordering from @suitesparse. */
        MESS_REORDER_AMD  = 2,
        /** Select the COLAMD reordering from @suitesparse. */
        MESS_REORDER_COLAMD = 3
    } mess_reorder_t;

    mess_int_t mess_matrix_memsize(mess_matrix matrix);
    mess_int_t mess_matrix_memsize_nnz(mess_int_t rows, mess_int_t cols, mess_int_t nnz, mess_int_t type, mess_int_t dtype);
    int mess_matrix_spy(mess_matrix matrix, const char* filename,int width, int height);
    int mess_matrix_reorder(mess_reorder_t reorder, mess_matrix matrix, mess_int_t *p, mess_int_t *q );
    int mess_matrix_reorder_amd(mess_matrix matrix, mess_int_t *p);
    int mess_matrix_reorder_colamd(mess_matrix matrix, mess_int_t *p);
    int mess_matrix_reorder_rcm(mess_matrix A, mess_int_t * perm);
    int mess_matrix_perm(mess_matrix matrix, mess_int_t *p, mess_int_t *q);
    int mess_matrix_permcopy(mess_matrix matrix, mess_int_t *p, mess_int_t *q, mess_matrix out);
    int mess_matrix_colperm(mess_matrix A, mess_int_t *perm);
    int mess_matrix_colpermcopy(mess_matrix A, mess_int_t *perm, mess_matrix B);
    int mess_matrix_joinpatterns(mess_matrix A, mess_matrix B);
    int mess_matrix_equalpattern(mess_matrix A, mess_matrix B);
    int mess_vector_tomatrix(mess_vector v, mess_matrix mat);
    int mess_vector_frommatrix(mess_matrix mat, mess_vector v);
    int mess_matrix_equal(mess_matrix mat1, mess_matrix mat2);
    int mess_matrix_equal_verbose(mess_matrix mat1, mess_matrix mat2 );
    int mess_matrix_decomp(mess_matrix A, mess_matrix Asym, mess_matrix Askewsym);

    /** @} */



    /** @addtogroup matrix_col
     * @{ */
    /*-----------------------------------------------------------------------------
     *  column operations
     *-----------------------------------------------------------------------------*/
    int mess_matrix_colnorm(mess_matrix Q, mess_int_t col, double *norm);
    int mess_matrix_colscale(mess_matrix Q, mess_int_t col, mess_double_cpx_t scale);
    int mess_matrix_coldot(mess_matrix Q, mess_int_t col1, mess_int_t col2, double *dot);
    int mess_matrix_coldotc(mess_matrix Q, mess_int_t col1, mess_int_t col2, mess_double_cpx_t *dot);
    int mess_matrix_colaxpy(mess_double_cpx_t alpha, mess_vector vect, mess_int_t col, mess_matrix Q);
    int mess_matrix_colaxpy2(mess_matrix Q, mess_double_cpx_t alpha, mess_int_t colc, mess_int_t col1);
    int mess_matrix_colvecdot(mess_matrix Q, mess_int_t col, mess_vector v, double *dot);
    int mess_matrix_colvecdotc(mess_matrix Q, mess_int_t col, mess_vector v, mess_double_cpx_t *dot);
    int mess_matrix_colvecaxpy(mess_double_cpx_t alpha, mess_int_t col, mess_matrix Q ,mess_vector v);
    int mess_matrix_coldotE(mess_matrix Q, mess_matrix E, mess_int_t col1, mess_int_t col2, double *dot );
    int mess_matrix_coldotcE(mess_matrix Q, mess_matrix E, mess_int_t col1, mess_int_t col2, mess_double_cpx_t *dot );
    int mess_matrix_colnormE(mess_matrix Q, mess_matrix E,mess_int_t col, double *dot);


    /** @}  */

    // Graph algorihtms
    /** @addtogroup graph
     *  @{ */
    int mess_graph_reach(mess_matrix G, mess_matrix B, mess_int_t k, mess_int_t *top, mess_int_t *xi, const mess_int_t *pinv);
    int mess_graph_dfs(mess_matrix G, mess_int_t j, mess_int_t* top, mess_int_t* xi, mess_int_t *pstack, const mess_int_t *pinv);
    int mess_graph_lblockpart(mess_matrix G, mess_int_t * perm, mess_int_t *nblocks, mess_int_t *blocks );
    /** @}  */

    /**
     * @addtogroup matrix_mvpcall
     * @{
     */


    /**
     * @brief Represent a generic matrix-vector-product.
     *
     * The @ref mess_mvpcall object represents a generic matrix-vector product \f$ y=Ax \f$
     * where \f$ A \f$ is not longer restricted to be a simple matrix. \n
     * This is necessary to use for example iterative solvers or iterative eigenvalue computations on
     * more complex objects than a simple matrix.
     * \see eigenvalues_sparse
     * \see itsolver_sol
     *
     * */
    typedef struct mess_mvpcall_st {
        /** Dimension of the matrix-vector product.  */
        mess_int_t dim;
        /** Data type of the matrix-vector product. \see mess_datatype_t.  */
        mess_datatype_t data_type;
        /** Function pointer which points to the implementation of the matrix-vector product.
         * It has to compute \f$ y = Ax \f$. */
        int (*mvp)(void *data, mess_operation_t op, mess_vector x, mess_vector y);
        /** Pointer to the data needed by the matrix-vector product.  */
        void *data;
    } *  mess_mvpcall;

    int mess_mvpcall_matrix(mess_mvpcall * mvpcall, mess_operation_t op, mess_matrix A);
    int mess_mvpcall_operator(mess_mvpcall *mvpcall, mess_int_t dim, mess_datatype_t data_type, int (*mvp)(void *data, mess_operation_t op, mess_vector x, mess_vector y), void *data);
    int mess_mvpcall_apply(mess_mvpcall mvpcall, mess_operation_t op, mess_vector x, mess_vector y);
    int mess_mvpcall_clear(mess_mvpcall *mvpcall);
    /** @}  */


    /** @addtogroup matgen_fdm
     * @{ */
    /**
     * @brief Type definition for the \f$ f_x \f$, \f$ f_y \f$ and \f$ g \f$  function in the FDM operator.
     *
     * The type definition creates functions \f$ f_x, f_y \f$ and \f$ g \f$ of the partial differential equation
     * (PDE)
     * \f[ \Delta u - f_x \frac{\mathrm{d}u}{\mathrm{d}x} - f_y \frac{\mathrm{d}u}{\mathrm{d}y} - gu= RHS \f]
     * on \f$ \Omega = (0,1) \times (0,1) \f$ with boundary condition
     * \f[ u=0  \f]
     * on \f$ d \Omega \f$.
     *
     * */
    typedef double (*mess_matgen_fdm_function)(double x, double y);

    int mess_matgen_fdmmatrix(mess_matrix A, mess_int_t n0, mess_matgen_fdm_function fnfx, mess_matgen_fdm_function fnfy, mess_matgen_fdm_function fng );
    int mess_matgen_fdmvector(mess_vector v, mess_int_t n0, mess_matgen_fdm_function func );
    int mess_matgen_fdmcolumn(mess_matrix B, mess_int_t n0, mess_matgen_fdm_function func );
    int mess_matgen_fdmrow(mess_matrix C, mess_int_t n0, mess_matgen_fdm_function func  );


    /** @}  */



#ifdef __cplusplus
}
#endif

#endif /*MESS_H_*/

/** \}@ */
