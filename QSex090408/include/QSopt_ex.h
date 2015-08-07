/****************************************************************************/
/*                                                                          */
/*  This file is part of QSopt_ex.                                          */
/*                                                                          */
/*  (c) Copyright 2006-2008 by David Applegate, William Cook, Sanjeeb Dash, */
/*  and Daniel Espinoza                                                     */
/*                                                                          */
/*  This code may be used under the terms of the GNU General Public License */
/*  (Version 2.1 or later) as published by the Free Software Foundation.    */
/*                                                                          */
/*  Alternatively, use is granted for research purposes only.               */ 
/*                                                                          */
/*  It is your choice of which of these two licenses you are operating      */
/*  under.                                                                  */
/*                                                                          */
/*  We make no guarantees about the correctness or usefulness of this code. */
/*                                                                          */
/****************************************************************************/

#ifndef __QSopt_ex__
#define __QSopt_ex__

#include <stdlib.h>
#include <gmp.h>

#ifndef __BASICDEFS__
#define __BASICDEFS__

/* algo */
#define PRIMAL_SIMPLEX 1
#define DUAL_SIMPLEX   2
#define PRIMAL_OR_DUAL 3


#ifndef __QS_BASIS__
#define __QS_BASIS__
typedef struct qsbasis
{
    int nstruct;
    int nrows;
    char *cstat;
    char *rstat;
}
QSbasis;
#endif

extern mpq_t mpq_ILL_MINDOUBLE;
extern mpq_t mpq_ILL_MAXDOUBLE;

#ifndef QS_DEFINITIONS
#define QS_DEFINITIONS
#define QS_MIN       (1)
#define QS_MAX       (-1)


/****************************************************************************/
/*                                                                          */
/*                 PARAMETERS THAT CAN BE SET BY setparam                   */
/*                                                                          */
/****************************************************************************/


#define QS_PARAM_PRIMAL_PRICING    0
#define QS_PARAM_DUAL_PRICING      2
#define QS_PARAM_SIMPLEX_DISPLAY   4
#define QS_PARAM_SIMPLEX_MAX_ITERATIONS 5
#define QS_PARAM_SIMPLEX_MAX_TIME  6
#define QS_PARAM_SIMPLEX_SCALING   7


/****************************************************************************/
/*                                                                          */
/*                     VALUES FOR PRICING PARAMETERS                        */
/*                                                                          */
/****************************************************************************/

#define QS_PRICE_PDANTZIG 1
#define QS_PRICE_PDEVEX 2
#define QS_PRICE_PSTEEP 3
#define QS_PRICE_PMULTPARTIAL 4

#define QS_PRICE_DDANTZIG 6
#define QS_PRICE_DSTEEP 7
#define QS_PRICE_DMULTPARTIAL 8
#define QS_PRICE_DDEVEX 9


/****************************************************************************/
/*                                                                          */
/*                         VALUES FOR BASIS STATUS                          */
/*                                                                          */
/****************************************************************************/


#define QS_COL_BSTAT_LOWER     '0'
#define QS_COL_BSTAT_BASIC     '1'
#define QS_COL_BSTAT_UPPER     '2'
#define QS_COL_BSTAT_FREE      '3'

#define QS_ROW_BSTAT_LOWER     '0'
#define QS_ROW_BSTAT_BASIC     '1'
#define QS_ROW_BSTAT_UPPER     '2'


/****************************************************************************/
/*                                                                          */
/*    Return Status for QSopt_primal, QSopt_dual, QSget_status              */
/*                                                                          */
/****************************************************************************/

#define QS_LP_OPTIMAL           1
#define QS_LP_INFEASIBLE        2
#define QS_LP_UNBOUNDED         3
#define QS_LP_ITER_LIMIT        4
#define QS_LP_TIME_LIMIT        5
#define QS_LP_UNSOLVED          6
#define QS_LP_ABORTED           7
#define QS_LP_MODIFIED        100
#endif


/**** error collection ****/
#define QS_DATA_ERROR            0
#define QS_DATA_WARN            1
#define QS_MPS_FORMAT_ERROR        2
#define QS_MPS_FORMAT_WARN        3
#define QS_LP_FORMAT_ERROR        4
#define QS_LP_FORMAT_WARN        5
#define QS_INPUT_NERROR        8

/* defs for phase I ratio test */
#define BBOUND    1
#define BATOLOWER 2
#define BATOUPPER 3
#define BBTOLOWER 4
#define BBTOUPPER 5
#define BSKIP     6

/* result of ratio test */
#define RATIO_UNBOUNDED 1
#define RATIO_NOBCHANGE 2
#define RATIO_BCHANGE   3
#define RATIO_FAILED    4
#define RATIO_NEGATIVE  5

#endif

/****************************************************************************/
/****************************************************************************/

#ifndef mpq___QS_QSTRUCT_H
#define mpq___QS_QSTRUCT_H

typedef struct mpq_qsdata {
	struct mpq_ILLlpdata *qslp;
	struct mpq_lpinfo *lp;
	struct mpq_price_info *pricing;
	struct mpq_ILLlp_basis *basis;
	struct mpq_ILLlp_cache *cache;
	char *name;
	int qstatus;          /* QS_LP_UNSOLVED or an opt status  */
	int factorok;         /* set to 0 if factorization is old */
	int simplex_display;  /* 0 off, 1 on */
	int simplex_scaling;  /* 0 off, 1 on */
} mpq_QSdata;

#endif /* mpq___QS_QSTRUCT_H */

/****************************************************************************/
/****************************************************************************/

#ifndef dbl___QS_QSTRUCT_H
#define dbl___QS_QSTRUCT_H

typedef struct dbl_qsdata {
	struct dbl_ILLlpdata *qslp;
	struct dbl_lpinfo *lp;
	struct dbl_price_info *pricing;
	struct dbl_ILLlp_basis *basis;
	struct dbl_ILLlp_cache *cache;
	char *name;
	int qstatus;          /* QS_LP_UNSOLVED or an opt status  */
	int factorok;         /* set to 0 if factorization is old */
	int simplex_display;  /* 0 off, 1 on */
	int simplex_scaling;  /* 0 off, 1 on */
} dbl_QSdata;

#endif /* dbl___QS_QSTRUCT_H */

/****************************************************************************/
/****************************************************************************/

#ifndef mpf___QS_QSTRUCT_H
#define mpf___QS_QSTRUCT_H

typedef struct mpf_qsdata {
	struct mpf_ILLlpdata *qslp;
	struct mpf_lpinfo *lp;
	struct mpf_price_info *pricing;
	struct mpf_ILLlp_basis *basis;
	struct mpf_ILLlp_cache *cache;
	char *name;
	int qstatus;          /* QS_LP_UNSOLVED or an opt status  */
	int factorok;         /* set to 0 if factorization is old */
	int simplex_display;  /* 0 off, 1 on */
	int simplex_scaling;  /* 0 off, 1 on */
} mpf_QSdata;

#endif /* mpf___QS_QSTRUCT_H */


/****************************************************************************/
/****************************************************************************/

#ifndef dbl___QS_QSOPT_H
#define dbl___QS_QSOPT_H


#ifdef WIN32

#ifdef QSLIB_EXPORTS
#define dbl_QSLIB_INTERFACE __declspec(dllexport)
#else
#define dbl_QSLIB_INTERFACE __declspec(dllimport)
#endif

#else
#define dbl_QSLIB_INTERFACE extern
#endif

#ifdef WIN32
typedef struct dbl_QSLIB_INTERFACE dbl_qsdata *dbl_QSprob;
typedef struct dbl_QSLIB_INTERFACE qsbasis *dbl_QSbas;
#else
typedef struct dbl_qsdata *dbl_QSprob;
typedef struct qsbasis *dbl_QSbas;
#endif

/****************************************************************************/
/*                                                                          */
/*                      QSopt Library Functions                             */
/*                                                                          */
/****************************************************************************/
#ifdef  __cplusplus
extern "C"
{
#endif

#ifdef WIN32
/* 
 *  in WINDOWS we make 
 *     dbl_solver_main/dbl_reader_main part of DLL
 */
    dbl_QSLIB_INTERFACE int dbl_solver_main (int argc, char **argv);
    dbl_QSLIB_INTERFACE int dbl_reader_main (int argc, char **argv);
#endif

    dbl_QSLIB_INTERFACE void dbl_QSfree (void *ptr),
    dbl_QSfree_prob (dbl_QSprob p),
    dbl_QSfree_basis (dbl_QSbas B),
    dbl_QSset_precision (const unsigned prec),
    dbl_QSstart (void), 
    dbl_QSend (void); 
    dbl_QSLIB_INTERFACE int dbl_QSopt_primal (dbl_QSprob p, int *status),
    dbl_QSopt_dual (dbl_QSprob p, int *status),
    dbl_QSopt_pivotin_col (dbl_QSprob p, int ccnt, int *clist),
    dbl_QSopt_pivotin_row (dbl_QSprob p, int rcnt, int *rlist),
    dbl_QSopt_strongbranch (dbl_QSprob p, int ncand, int *candidatelist,
        double * xlist, double * down_vals, double * up_vals, int iterations,
        double objbound),
    dbl_QSchange_objsense (dbl_QSprob p, int newsense),
    dbl_QSget_objsense (dbl_QSprob p, int *newsense),
    dbl_QSnew_col (dbl_QSprob p, double obj, double lower, double upper,
        const char *name),
    dbl_QSadd_cols (dbl_QSprob p, int num, int *cmatcnt, int *cmatbeg,
        int *cmatind, double * cmatval, double * obj, double * lower,
        double * upper, const char **names),
    dbl_QSadd_col (dbl_QSprob p, int cnt, int *cmatind, double * cmatval,
        double obj, double lower, double upper, const char *name),
    dbl_QSnew_row (dbl_QSprob p, double rhs, int sense, const char *name),
    dbl_QSadd_rows (dbl_QSprob p, int num, int *rmatcnt, int *rmatbeg,
        int *rmatind, double * rmatval, double * rhs, char *sense,
        const char **names),
    dbl_QSadd_row (dbl_QSprob p, int cnt, int *rmatind, double * rmatval,
        double * rhs, int sense, const char *name),
    dbl_QSdelete_rows (dbl_QSprob p, int num, int *dellist),
    dbl_QSdelete_row (dbl_QSprob p, int rowindex),
    dbl_QSdelete_setrows (dbl_QSprob p, int *flags),
    dbl_QSdelete_named_row (dbl_QSprob p, const char *rowname),
    dbl_QSdelete_named_rows_list (dbl_QSprob p, int num,
        const char **rownames),
    dbl_QSdelete_cols (dbl_QSprob p, int num, int *dellist),
    dbl_QSdelete_col (dbl_QSprob p, int colindex),
    dbl_QSdelete_setcols (dbl_QSprob p, int *flags),
    dbl_QSdelete_named_column (dbl_QSprob p, const char *colname),
    dbl_QSdelete_named_columns_list (dbl_QSprob p, int num,
        const char **colnames),
    dbl_QSchange_senses (dbl_QSprob p, int num, int *rowlist, char *sense),
    dbl_QSchange_sense (dbl_QSprob p, int rowindex, int sense),
    dbl_QSchange_coef (dbl_QSprob p, int rowindex, int colindex, double coef),
    dbl_QSchange_objcoef (dbl_QSprob p, int indx, double coef),
    dbl_QSchange_rhscoef (dbl_QSprob p, int indx, double coef),
    dbl_QSchange_bounds (dbl_QSprob p, int num, int *collist, char *lu,
        double * bounds),
    dbl_QSchange_bound (dbl_QSprob p, int indx, int lu, double bound),
    dbl_QSload_basis (dbl_QSprob p, dbl_QSbas B),
    dbl_QSread_and_load_basis (dbl_QSprob p, const char *filename),
    dbl_QSload_basis_array (dbl_QSprob p, char *cstat, char *rstat),
    dbl_QSload_basis_and_row_norms_array (dbl_QSprob p, char *cstat,
        char *rstat, double * rownorms),
    dbl_QSget_basis_array (dbl_QSprob p, char *cstat, char *rstat),
    dbl_QSget_basis_and_row_norms_array (dbl_QSprob p, char *cstat,
        char *rstat, double * rownorms),
    dbl_QSget_binv_row (dbl_QSprob p, int indx, double * binvrow),
    dbl_QSget_tableau_row (dbl_QSprob p, int indx, double * tableaurow),
    dbl_QSget_basis_order (dbl_QSprob p, int *basorder),
    dbl_QSget_status (dbl_QSprob p, int *status),
    dbl_QSget_solution (dbl_QSprob p, double * value, double * x, double * pi,
        double * slack, double * rc),
    dbl_QSget_objval (dbl_QSprob p, double * value),
    dbl_QSget_pi_array (dbl_QSprob p, double * pi),
    dbl_QSget_rc_array (dbl_QSprob p, double * rc),
    dbl_QSget_x_array (dbl_QSprob p, double * x),
    dbl_QSget_slack_array (dbl_QSprob p, double * slack),
    dbl_QSget_infeas_array (dbl_QSprob p, double * pi),
    dbl_QSget_colcount (dbl_QSprob p),
    dbl_QSget_rowcount (dbl_QSprob p),
    dbl_QSget_nzcount (dbl_QSprob p),
    dbl_QSget_obj (dbl_QSprob p, double * obj),
    dbl_QSget_rhs (dbl_QSprob p, double * rhs),
    dbl_QSget_rows_list (dbl_QSprob p, int num, int *rowlist, int **rowcnt,
        int **rowbeg, int **rowind, double ** rowval, double ** rhs,
        char **sense, char ***names),
    dbl_QSget_rows (dbl_QSprob p, int **rowcnt, int **rowbeg, int **rowind,
        double ** rowval, double ** rhs, char **sense, char ***names),
    dbl_QSget_columns_list (dbl_QSprob p, int num, int *collist, int **colcnt,
        int **colbeg, int **colind, double ** colval, double ** obj,
        double ** lower, double ** upper, char ***names),
    dbl_QSget_columns (dbl_QSprob p, int **colcnt, int **colbeg, int **colind,
        double ** colval, double ** obj, double ** lower, double ** upper,
        char ***names),
    dbl_QSget_rownames (dbl_QSprob p, char **rownames),
    dbl_QSget_colnames (dbl_QSprob p, char **colnames),
    dbl_QSget_bound (dbl_QSprob p, int colindex, int lu, double * bound),
    dbl_QSget_bounds (dbl_QSprob p, double * lower, double * upper),
    dbl_QSget_intflags (dbl_QSprob p, int *intflags),
    dbl_QSget_intcount (dbl_QSprob p, int *count),
    dbl_QSget_column_index (dbl_QSprob p, const char *name, int *colindex),
    dbl_QSget_row_index (dbl_QSprob p, const char *name, int *rowindex),
    dbl_QSget_named_x (dbl_QSprob p, const char *colname, double * val),
    dbl_QSget_named_rc (dbl_QSprob p, const char *colname, double * val),
    dbl_QSget_named_pi (dbl_QSprob p, const char *rowname, double * val),
    dbl_QSget_named_slack (dbl_QSprob p, const char *rowname, double * val),
    dbl_QScompute_row_norms (dbl_QSprob p),
    dbl_QSwrite_prob (dbl_QSprob p, const char *filename,
        const char *filetype),
    dbl_QSwrite_prob_file (dbl_QSprob p, FILE * file, const char *filetype),
    dbl_QSwrite_basis (dbl_QSprob p, dbl_QSbas B, const char *filename),
    dbl_QStest_row_norms (dbl_QSprob p),
    dbl_QSset_param (dbl_QSprob p, int whichparam, int newvalue),
    dbl_QSset_param_EGlpNum (dbl_QSprob p, int whichparam, double newvalue),
    dbl_QSget_param (dbl_QSprob p, int whichparam, int *value),
    dbl_QSget_param_EGlpNum (dbl_QSprob p, int whichparam, double * value);

    dbl_QSLIB_INTERFACE char *dbl_QSget_probname (dbl_QSprob p);
    dbl_QSLIB_INTERFACE char *dbl_QSget_objname (dbl_QSprob p);
    dbl_QSLIB_INTERFACE char *dbl_QSversion (void);

    dbl_QSLIB_INTERFACE dbl_QSprob dbl_QScreate_prob (const char *name,
        int objsense),
    dbl_QSread_prob (const char *filename, const char *filetype),
    dbl_QSload_prob (const char *probname, int ncols, int nrows, int *cmatcnt,
        int *cmatbeg, int *cmatind, double * cmatval, int objsense,
        double * obj, double * rhs, char *sense, double * lower,
        double * upper, const char **colnames, const char **rownames),
    dbl_QScopy_prob (dbl_QSprob p, const char *newname);

    dbl_QSLIB_INTERFACE dbl_QSbas dbl_QSget_basis (dbl_QSprob p),
    dbl_QSread_basis (dbl_QSprob p, const char *filename);

#ifdef  __cplusplus
}
#endif

/****************************************************************************
 *
 * This is the undocumented part of the QSlib interface 
 *
 ****************************************************************************/
/* 
 * functions to facilitate line by line reading from other sources than 
 * files from within MPS/LP parsers  
 * 
 * functions to facilitate the collection of error information instead of 
 * having the parsers print messages to stderr
 *                              by mps/lp format writers
 * 
 * a problem's reporter is used by the solver code to provide important 
 * dbl_feedback/progress information
 */

#ifdef WIN32
typedef struct dbl_QSLIB_INTERFACE dbl_qsline_reader *dbl_QSline_reader;
typedef struct dbl_QSLIB_INTERFACE dbl_qsformat_error *dbl_QSformat_error;
typedef struct dbl_QSLIB_INTERFACE dbl_qserror_collector *dbl_QSerror_collector;
typedef struct dbl_QSLIB_INTERFACE dbl_qserror_memory *dbl_QSerror_memory;
#else
typedef struct dbl_qsline_reader *dbl_QSline_reader;
typedef struct dbl_qsformat_error *dbl_QSformat_error;
typedef struct dbl_qserror_collector *dbl_QSerror_collector;
typedef struct dbl_qserror_memory *dbl_QSerror_memory;
#endif

#ifdef  __cplusplus
extern "C"
{
#endif
    dbl_QSLIB_INTERFACE const char *dbl_QSformat_error_type_string (int tp);

    dbl_QSLIB_INTERFACE int dbl_QSerror_get_type (dbl_QSformat_error error);
    dbl_QSLIB_INTERFACE const char *dbl_QSerror_get_desc (dbl_QSformat_error error);
    dbl_QSLIB_INTERFACE int dbl_QSerror_get_line_number (dbl_QSformat_error error);
    dbl_QSLIB_INTERFACE int dbl_QSerror_get_pos (dbl_QSformat_error error);
    dbl_QSLIB_INTERFACE const char *dbl_QSerror_get_line (dbl_QSformat_error error);
    dbl_QSLIB_INTERFACE void dbl_QSerror_print (FILE * f,
                                                                            dbl_QSformat_error error);

    dbl_QSLIB_INTERFACE dbl_QSerror_collector dbl_QSerror_collector_new (void *fct,
        void *dest);
    dbl_QSLIB_INTERFACE
    dbl_QSerror_collector dbl_QSerror_memory_collector_new (dbl_QSerror_memory mem);
    dbl_QSLIB_INTERFACE void dbl_QSerror_collector_free (dbl_QSerror_collector c);

/****************************************************************************
 * line reader 
 */
    dbl_QSLIB_INTERFACE dbl_QSline_reader dbl_QSline_reader_new (void *fct,
         void *data_src);

    dbl_QSLIB_INTERFACE void dbl_QSline_reader_free (dbl_QSline_reader reader);

    dbl_QSLIB_INTERFACE void dbl_QSline_reader_set_error_collector (dbl_QSline_reader reader,
        dbl_QSerror_collector collector);

    dbl_QSLIB_INTERFACE char *dbl_QSline_reader_get (dbl_QSline_reader reader,
        char *s, int size);

    dbl_QSLIB_INTERFACE dbl_QSprob dbl_QSget_prob (dbl_QSline_reader reader,
        const char *probname, const char *filetype);


    dbl_QSLIB_INTERFACE dbl_QSerror_memory dbl_QSerror_memory_create (int takeErrorLines);
    dbl_QSLIB_INTERFACE void dbl_QSerror_memory_free (dbl_QSerror_memory mem);

    dbl_QSLIB_INTERFACE int dbl_QSerror_memory_get_nof (dbl_QSerror_memory mem,
        int error_type);
    dbl_QSLIB_INTERFACE int dbl_QSerror_memory_get_nerrors (dbl_QSerror_memory mem);

    dbl_QSLIB_INTERFACE dbl_QSformat_error
        dbl_QSerror_memory_get_last_error (dbl_QSerror_memory mem);
    dbl_QSLIB_INTERFACE dbl_QSformat_error
        dbl_QSerror_memory_get_prev_error (dbl_QSformat_error e);

    dbl_QSLIB_INTERFACE void dbl_QSset_reporter (dbl_QSprob prob,
        int iterskip, void *fct, void *dest);

    dbl_QSLIB_INTERFACE int dbl_QSreport_prob (dbl_QSprob p,
        const char *filetype, dbl_QSerror_collector c);

#ifdef  __cplusplus
}
#endif
#endif  /* dbl___QS_QSOPT_H */

/****************************************************************************/
/****************************************************************************/

#ifndef mpq___QS_QSOPT_H
#define mpq___QS_QSOPT_H


#ifdef WIN32

#ifdef QSLIB_EXPORTS
#define mpq_QSLIB_INTERFACE __declspec(dllexport)
#else
#define mpq_QSLIB_INTERFACE __declspec(dllimport)
#endif

#else
#define mpq_QSLIB_INTERFACE extern
#endif

#ifdef WIN32
typedef struct mpq_QSLIB_INTERFACE mpq_qsdata *mpq_QSprob;
typedef struct mpq_QSLIB_INTERFACE qsbasis *mpq_QSbas;
#else
typedef struct mpq_qsdata *mpq_QSprob;
typedef struct qsbasis *mpq_QSbas;
#endif

/****************************************************************************/
/*                                                                          */
/*                      QSopt Library Functions                             */
/*                                                                          */
/****************************************************************************/
#ifdef  __cplusplus
extern "C"
{
#endif

#ifdef WIN32
/* 
 *  in WINDOWS we make 
 *     mpq_solver_main/mpq_reader_main part of DLL
 */
    mpq_QSLIB_INTERFACE int mpq_solver_main (int argc, char **argv);
    mpq_QSLIB_INTERFACE int mpq_reader_main (int argc, char **argv);
#endif

    mpq_QSLIB_INTERFACE void mpq_QSfree (void *ptr),
    mpq_QSfree_prob (mpq_QSprob p),
    mpq_QSfree_basis (mpq_QSbas B),
    mpq_QSset_precision (const unsigned prec),
       /* set the mpq_precision for floats to the given number of bits */

    mpq_QSstart (void),  
       /* when we use non native numbers, we need to make some */
       /* initializations before operating with the library */

    mpq_QSend (void);
       /* to free any internal static data needed by the variable */
       /* mpq_precision numbers */

    mpq_QSLIB_INTERFACE int mpq_QSopt_primal (mpq_QSprob p, int *status),
    mpq_QSopt_dual (mpq_QSprob p, int *status),
    mpq_QSopt_pivotin_col (mpq_QSprob p, int ccnt, int *clist),
    mpq_QSopt_pivotin_row (mpq_QSprob p, int rcnt, int *rlist),
    mpq_QSopt_strongbranch (mpq_QSprob p, int ncand, int *candidatelist,
        mpq_t * xlist, mpq_t * down_vals, mpq_t * up_vals,
        int iterations, mpq_t objbound),
    mpq_QSchange_objsense (mpq_QSprob p, int newsense),
    mpq_QSget_objsense (mpq_QSprob p, int *newsense),
    mpq_QSnew_col (mpq_QSprob p, mpq_t obj, mpq_t lower, mpq_t upper,
        const char *name),
    mpq_QSadd_cols (mpq_QSprob p, int num, int *cmatcnt, int *cmatbeg,
        int *cmatind, mpq_t * cmatval, mpq_t * obj, mpq_t * lower,
        mpq_t * upper, const char **names),
    mpq_QSadd_col (mpq_QSprob p, int cnt, int *cmatind, mpq_t * cmatval,
        mpq_t obj, mpq_t lower, mpq_t upper, const char *name),
    mpq_QSnew_row (mpq_QSprob p, mpq_t rhs, int sense, const char *name),
    mpq_QSadd_rows (mpq_QSprob p, int num, int *rmatcnt, int *rmatbeg,
        int *rmatind, mpq_t * rmatval, mpq_t * rhs, char *sense,
        const char **names),
    mpq_QSadd_row (mpq_QSprob p, int cnt, int *rmatind, mpq_t * rmatval,
         mpq_t * rhs, int sense, const char *name),
    mpq_QSdelete_rows (mpq_QSprob p, int num, int *dellist),
    mpq_QSdelete_row (mpq_QSprob p, int rowindex),
    mpq_QSdelete_setrows (mpq_QSprob p, int *flags),
    mpq_QSdelete_named_row (mpq_QSprob p, const char *rowname),
    mpq_QSdelete_named_rows_list (mpq_QSprob p, int num,
        const char **rownames),
    mpq_QSdelete_cols (mpq_QSprob p, int num, int *dellist),
    mpq_QSdelete_col (mpq_QSprob p, int colindex),
    mpq_QSdelete_setcols (mpq_QSprob p, int *flags),
    mpq_QSdelete_named_column (mpq_QSprob p, const char *colname),
    mpq_QSdelete_named_columns_list (mpq_QSprob p, int num,
        const char **colnames),
    mpq_QSchange_senses (mpq_QSprob p, int num, int *rowlist, char *sense),
    mpq_QSchange_sense (mpq_QSprob p, int rowindex, int sense),
    mpq_QSchange_coef (mpq_QSprob p, int rowindex, int colindex, mpq_t coef),
    mpq_QSchange_objcoef (mpq_QSprob p, int indx, mpq_t coef),
    mpq_QSchange_rhscoef (mpq_QSprob p, int indx, mpq_t coef),
    mpq_QSchange_bounds (mpq_QSprob p, int num, int *collist, char *lu,
        mpq_t * bounds),
    mpq_QSchange_bound (mpq_QSprob p, int indx, int lu, mpq_t bound),
    mpq_QSload_basis (mpq_QSprob p, mpq_QSbas B),
    mpq_QSread_and_load_basis (mpq_QSprob p, const char *filename),
    mpq_QSload_basis_array (mpq_QSprob p, char *cstat, char *rstat),
    mpq_QSload_basis_and_row_norms_array (mpq_QSprob p, char *cstat,
        char *rstat, mpq_t * rownorms),
    mpq_QSget_basis_array (mpq_QSprob p, char *cstat, char *rstat),
    mpq_QSget_basis_and_row_norms_array (mpq_QSprob p, char *cstat,
        char *rstat, mpq_t * rownorms),
    mpq_QSget_binv_row (mpq_QSprob p, int indx, mpq_t * binvrow),
    mpq_QSget_tableau_row (mpq_QSprob p, int indx, mpq_t * tableaurow),
    mpq_QSget_basis_order (mpq_QSprob p, int *basorder),
    mpq_QSget_status (mpq_QSprob p, int *status),
    mpq_QSget_solution (mpq_QSprob p, mpq_t * value, mpq_t * x,
        mpq_t * pi, mpq_t * slack, mpq_t * rc),
    mpq_QSget_objval (mpq_QSprob p, mpq_t * value),
    mpq_QSget_pi_array (mpq_QSprob p, mpq_t * pi),
    mpq_QSget_rc_array (mpq_QSprob p, mpq_t * rc),
    mpq_QSget_x_array (mpq_QSprob p, mpq_t * x),
    mpq_QSget_slack_array (mpq_QSprob p, mpq_t * slack),
    mpq_QSget_infeas_array (mpq_QSprob p, mpq_t * pi),
    mpq_QSget_colcount (mpq_QSprob p),
    mpq_QSget_rowcount (mpq_QSprob p),
    mpq_QSget_nzcount (mpq_QSprob p),
    mpq_QSget_obj (mpq_QSprob p, mpq_t * obj),
    mpq_QSget_rhs (mpq_QSprob p, mpq_t * rhs),
    mpq_QSget_rows_list (mpq_QSprob p, int num, int *rowlist, int **rowcnt,
        int **rowbeg, int **rowind, mpq_t ** rowval, mpq_t ** rhs,
        char **sense, char ***names),
    mpq_QSget_rows (mpq_QSprob p, int **rowcnt, int **rowbeg, int **rowind,
        mpq_t ** rowval, mpq_t ** rhs, char **sense, char ***names),
    mpq_QSget_columns_list (mpq_QSprob p, int num, int *collist, int **colcnt,
        int **colbeg, int **colind, mpq_t ** colval, mpq_t ** obj,
        mpq_t ** lower, mpq_t ** upper, char ***names),
    mpq_QSget_columns (mpq_QSprob p, int **colcnt, int **colbeg, int **colind,
        mpq_t ** colval, mpq_t ** obj, mpq_t ** lower, mpq_t ** upper,
        char ***names),
    mpq_QSget_rownames (mpq_QSprob p, char **rownames),
    mpq_QSget_colnames (mpq_QSprob p, char **colnames),
    mpq_QSget_bound (mpq_QSprob p, int colindex, int lu, mpq_t * bound),
    mpq_QSget_bounds (mpq_QSprob p, mpq_t * lower, mpq_t * upper),
    mpq_QSget_intflags (mpq_QSprob p, int *intflags),
    mpq_QSget_intcount (mpq_QSprob p, int *count),
    mpq_QSget_column_index (mpq_QSprob p, const char *name, int *colindex),
    mpq_QSget_row_index (mpq_QSprob p, const char *name, int *rowindex),
    mpq_QSget_named_x (mpq_QSprob p, const char *colname, mpq_t * val),
    mpq_QSget_named_rc (mpq_QSprob p, const char *colname, mpq_t * val),
    mpq_QSget_named_pi (mpq_QSprob p, const char *rowname, mpq_t * val),
    mpq_QSget_named_slack (mpq_QSprob p, const char *rowname, mpq_t * val),
    mpq_QScompute_row_norms (mpq_QSprob p),
    mpq_QSwrite_prob (mpq_QSprob p, const char *filename,
        const char *filetype),
    mpq_QSwrite_prob_file (mpq_QSprob p, FILE * file, const char *filetype),
    mpq_QSwrite_basis (mpq_QSprob p, mpq_QSbas B, const char *filename),
    mpq_QStest_row_norms (mpq_QSprob p),
    mpq_QSset_param (mpq_QSprob p, int whichparam, int newvalue),
    mpq_QSset_param_EGlpNum (mpq_QSprob p, int whichparam, mpq_t newvalue),
    mpq_QSget_param (mpq_QSprob p, int whichparam, int *value),
    mpq_QSget_param_EGlpNum (mpq_QSprob p, int whichparam, mpq_t * value);

    mpq_QSLIB_INTERFACE char *mpq_QSget_probname (mpq_QSprob p);
    mpq_QSLIB_INTERFACE char *mpq_QSget_objname (mpq_QSprob p);
    mpq_QSLIB_INTERFACE char *mpq_QSversion (void);

    mpq_QSLIB_INTERFACE mpq_QSprob mpq_QScreate_prob (const char *name,
        int objsense),
    mpq_QSread_prob (const char *filename, const char *filetype),
    mpq_QSload_prob (const char *probname, int ncols, int nrows, int *cmatcnt,
        int *cmatbeg, int *cmatind, mpq_t * cmatval, int objsense, mpq_t * obj,
        mpq_t * rhs, char *sense, mpq_t * lower, mpq_t * upper,
        const char **colnames, const char **rownames),
    mpq_QScopy_prob (mpq_QSprob p, const char *newname);

    mpq_QSLIB_INTERFACE mpq_QSbas mpq_QSget_basis (mpq_QSprob p),
    mpq_QSread_basis (mpq_QSprob p, const char *filename);

#ifdef  __cplusplus
}
#endif

/****************************************************************************
 *
 * This is the undocumented part of the QSlib interface 
 *
 ****************************************************************************/
/* 
 * functions to facilitate line by line reading from other sources than 
 * files from within MPS/LP parsers  
 * 
 * functions to facilitate the collection of error information instead of 
 * having the parsers print messages to stderr by mps/lp format writers
 * 
 * a problem's reporter is used by the solver code to provide important 
 * mpq_feedback/progress information
 */

#ifdef WIN32
typedef struct mpq_QSLIB_INTERFACE mpq_qsline_reader *mpq_QSline_reader;
typedef struct mpq_QSLIB_INTERFACE mpq_qsformat_error *mpq_QSformat_error;
typedef struct mpq_QSLIB_INTERFACE mpq_qserror_collector *mpq_QSerror_collector;
typedef struct mpq_QSLIB_INTERFACE mpq_qserror_memory *mpq_QSerror_memory;
#else
typedef struct mpq_qsline_reader *mpq_QSline_reader;
typedef struct mpq_qsformat_error *mpq_QSformat_error;
typedef struct mpq_qserror_collector *mpq_QSerror_collector;
typedef struct mpq_qserror_memory *mpq_QSerror_memory;
#endif

#ifdef  __cplusplus
extern "C"
{
#endif
    mpq_QSLIB_INTERFACE const char *mpq_QSformat_error_type_string (int tp);

    mpq_QSLIB_INTERFACE int mpq_QSerror_get_type (mpq_QSformat_error error);
    mpq_QSLIB_INTERFACE const char *mpq_QSerror_get_desc (mpq_QSformat_error error);
    mpq_QSLIB_INTERFACE int mpq_QSerror_get_line_number (mpq_QSformat_error error);
    mpq_QSLIB_INTERFACE int mpq_QSerror_get_pos (mpq_QSformat_error error);
    mpq_QSLIB_INTERFACE const char *mpq_QSerror_get_line (mpq_QSformat_error error);
    mpq_QSLIB_INTERFACE void mpq_QSerror_print (FILE * f,
        mpq_QSformat_error error);

    mpq_QSLIB_INTERFACE mpq_QSerror_collector mpq_QSerror_collector_new (void *fct,
         void *dest);
    mpq_QSLIB_INTERFACE
    mpq_QSerror_collector mpq_QSerror_memory_collector_new (mpq_QSerror_memory mem);
    mpq_QSLIB_INTERFACE void mpq_QSerror_collector_free (mpq_QSerror_collector c);

    mpq_QSLIB_INTERFACE mpq_QSline_reader mpq_QSline_reader_new (void *fct,
         void *data_src);

    mpq_QSLIB_INTERFACE void mpq_QSline_reader_free (mpq_QSline_reader reader);

    mpq_QSLIB_INTERFACE void mpq_QSline_reader_set_error_collector (mpq_QSline_reader reader,
        mpq_QSerror_collector collector);

    mpq_QSLIB_INTERFACE char *mpq_QSline_reader_get (mpq_QSline_reader reader,
        char *s, int size);

    mpq_QSLIB_INTERFACE mpq_QSprob mpq_QSget_prob (mpq_QSline_reader reader,
        const char *probname, const char *filetype);


    mpq_QSLIB_INTERFACE mpq_QSerror_memory mpq_QSerror_memory_create (int takeErrorLines);
    mpq_QSLIB_INTERFACE void mpq_QSerror_memory_free (mpq_QSerror_memory mem);

    mpq_QSLIB_INTERFACE int mpq_QSerror_memory_get_nof (mpq_QSerror_memory mem,
                                                                                            int error_type);
    mpq_QSLIB_INTERFACE int mpq_QSerror_memory_get_nerrors (mpq_QSerror_memory mem);

    mpq_QSLIB_INTERFACE mpq_QSformat_error
        mpq_QSerror_memory_get_last_error (mpq_QSerror_memory mem);
    mpq_QSLIB_INTERFACE mpq_QSformat_error
        mpq_QSerror_memory_get_prev_error (mpq_QSformat_error e);

    mpq_QSLIB_INTERFACE void mpq_QSset_reporter (mpq_QSprob prob,
        int iterskip, void *fct, void *dest);

    mpq_QSLIB_INTERFACE int mpq_QSreport_prob (mpq_QSprob p,
        const char *filetype, mpq_QSerror_collector c);

#ifdef  __cplusplus
}
#endif
#endif           /* mpq___QS_QSOPT_H */

/****************************************************************************/
/****************************************************************************/

#ifndef mpf___QS_QSOPT_H
#define mpf___QS_QSOPT_H


#ifdef WIN32

#ifdef QSLIB_EXPORTS
#define mpf_QSLIB_INTERFACE __declspec(dllexport)
#else
#define mpf_QSLIB_INTERFACE __declspec(dllimport)
#endif

#else
#define mpf_QSLIB_INTERFACE extern
#endif

#ifdef WIN32
typedef struct mpf_QSLIB_INTERFACE mpf_qsdata *mpf_QSprob;
typedef struct mpf_QSLIB_INTERFACE qsbasis *mpf_QSbas;
#else
typedef struct mpf_qsdata *mpf_QSprob;
typedef struct qsbasis *mpf_QSbas;
#endif

/****************************************************************************/
/*                                                                          */
/*                      QSopt Library Functions                             */
/*                                                                          */
/****************************************************************************/
#ifdef  __cplusplus
extern "C"
{
#endif

#ifdef WIN32
/* 
 *  in WINDOWS we make 
 *     mpf_solver_main/mpf_reader_main part of DLL
 */
    mpf_QSLIB_INTERFACE int mpf_solver_main (int argc, char **argv);
    mpf_QSLIB_INTERFACE int mpf_reader_main (int argc, char **argv);
#endif

    mpf_QSLIB_INTERFACE void mpf_QSfree (void *ptr),
    mpf_QSfree_prob (mpf_QSprob p),
    mpf_QSfree_basis (mpf_QSbas B),
    mpf_QSset_precision (const unsigned prec),

    mpf_QSstart (void), 
    mpf_QSend (void); 

    mpf_QSLIB_INTERFACE int mpf_QSopt_primal (mpf_QSprob p, int *status),
    mpf_QSopt_dual (mpf_QSprob p, int *status),
    mpf_QSopt_pivotin_col (mpf_QSprob p, int ccnt, int *clist),
    mpf_QSopt_pivotin_row (mpf_QSprob p, int rcnt, int *rlist),
    mpf_QSopt_strongbranch (mpf_QSprob p, int ncand, int *candidatelist,
        mpf_t * xlist, mpf_t * down_vals, mpf_t * up_vals, int iterations,
        mpf_t objbound),
    mpf_QSchange_objsense (mpf_QSprob p, int newsense),
    mpf_QSget_objsense (mpf_QSprob p, int *newsense),
    mpf_QSnew_col (mpf_QSprob p, mpf_t obj, mpf_t lower, mpf_t upper,
        const char *name),
    mpf_QSadd_cols (mpf_QSprob p, int num, int *cmatcnt, int *cmatbeg,
        int *cmatind, mpf_t * cmatval, mpf_t * obj, mpf_t * lower,
        mpf_t * upper, const char **names),
    mpf_QSadd_col (mpf_QSprob p, int cnt, int *cmatind, mpf_t * cmatval,
        mpf_t obj, mpf_t lower, mpf_t upper, const char *name),
    mpf_QSnew_row (mpf_QSprob p, mpf_t rhs, int sense, const char *name),
    mpf_QSadd_rows (mpf_QSprob p, int num, int *rmatcnt, int *rmatbeg,
        int *rmatind, mpf_t * rmatval, mpf_t * rhs, char *sense,
        const char **names),
    mpf_QSadd_row (mpf_QSprob p, int cnt, int *rmatind, mpf_t * rmatval,
        mpf_t * rhs, int sense, const char *name),
    mpf_QSdelete_rows (mpf_QSprob p, int num, int *dellist),
    mpf_QSdelete_row (mpf_QSprob p, int rowindex),
    mpf_QSdelete_setrows (mpf_QSprob p, int *flags),
    mpf_QSdelete_named_row (mpf_QSprob p, const char *rowname),
    mpf_QSdelete_named_rows_list (mpf_QSprob p, int num,
        const char **rownames),
    mpf_QSdelete_cols (mpf_QSprob p, int num, int *dellist),
    mpf_QSdelete_col (mpf_QSprob p, int colindex),
    mpf_QSdelete_setcols (mpf_QSprob p, int *flags),
    mpf_QSdelete_named_column (mpf_QSprob p, const char *colname),
    mpf_QSdelete_named_columns_list (mpf_QSprob p, int num,
        const char **colnames),
    mpf_QSchange_senses (mpf_QSprob p, int num, int *rowlist, char *sense),
    mpf_QSchange_sense (mpf_QSprob p, int rowindex, int sense),
    mpf_QSchange_coef (mpf_QSprob p, int rowindex, int colindex, mpf_t coef),
    mpf_QSchange_objcoef (mpf_QSprob p, int indx, mpf_t coef),
    mpf_QSchange_rhscoef (mpf_QSprob p, int indx, mpf_t coef),
    mpf_QSchange_bounds (mpf_QSprob p, int num, int *collist, char *lu,
        mpf_t * bounds),
    mpf_QSchange_bound (mpf_QSprob p, int indx, int lu, mpf_t bound),
    mpf_QSload_basis (mpf_QSprob p, mpf_QSbas B),
    mpf_QSread_and_load_basis (mpf_QSprob p, const char *filename),
    mpf_QSload_basis_array (mpf_QSprob p, char *cstat, char *rstat),
    mpf_QSload_basis_and_row_norms_array (mpf_QSprob p, char *cstat,
        char *rstat, mpf_t * rownorms),
    mpf_QSget_basis_array (mpf_QSprob p, char *cstat, char *rstat),
    mpf_QSget_basis_and_row_norms_array (mpf_QSprob p, char *cstat,
        char *rstat, mpf_t * rownorms),
    mpf_QSget_binv_row (mpf_QSprob p, int indx, mpf_t * binvrow),
    mpf_QSget_tableau_row (mpf_QSprob p, int indx, mpf_t * tableaurow),
    mpf_QSget_basis_order (mpf_QSprob p, int *basorder),
    mpf_QSget_status (mpf_QSprob p, int *status),
    mpf_QSget_solution (mpf_QSprob p, mpf_t * value, mpf_t * x, mpf_t * pi,
        mpf_t * slack, mpf_t * rc),
    mpf_QSget_objval (mpf_QSprob p, mpf_t * value),
    mpf_QSget_pi_array (mpf_QSprob p, mpf_t * pi),
    mpf_QSget_rc_array (mpf_QSprob p, mpf_t * rc),
    mpf_QSget_x_array (mpf_QSprob p, mpf_t * x),
    mpf_QSget_slack_array (mpf_QSprob p, mpf_t * slack),
    mpf_QSget_infeas_array (mpf_QSprob p, mpf_t * pi),
    mpf_QSget_colcount (mpf_QSprob p),
    mpf_QSget_rowcount (mpf_QSprob p),
    mpf_QSget_nzcount (mpf_QSprob p),
    mpf_QSget_obj (mpf_QSprob p, mpf_t * obj),
    mpf_QSget_rhs (mpf_QSprob p, mpf_t * rhs),
    mpf_QSget_rows_list (mpf_QSprob p, int num, int *rowlist, int **rowcnt,
        int **rowbeg, int **rowind, mpf_t ** rowval, mpf_t ** rhs,
        char **sense, char ***names),
    mpf_QSget_rows (mpf_QSprob p, int **rowcnt, int **rowbeg, int **rowind,
        mpf_t ** rowval, mpf_t ** rhs, char **sense, char ***names),
    mpf_QSget_columns_list (mpf_QSprob p, int num, int *collist, int **colcnt,
        int **colbeg, int **colind, mpf_t ** colval, mpf_t ** obj,
        mpf_t ** lower, mpf_t ** upper, char ***names),
    mpf_QSget_columns (mpf_QSprob p, int **colcnt, int **colbeg, int **colind,
        mpf_t ** colval, mpf_t ** obj, mpf_t ** lower, mpf_t ** upper,
        char ***names),
    mpf_QSget_rownames (mpf_QSprob p, char **rownames),
    mpf_QSget_colnames (mpf_QSprob p, char **colnames),
    mpf_QSget_bound (mpf_QSprob p, int colindex, int lu, mpf_t * bound),
    mpf_QSget_bounds (mpf_QSprob p, mpf_t * lower, mpf_t * upper),
    mpf_QSget_intflags (mpf_QSprob p, int *intflags),
    mpf_QSget_intcount (mpf_QSprob p, int *count),
    mpf_QSget_column_index (mpf_QSprob p, const char *name, int *colindex),
    mpf_QSget_row_index (mpf_QSprob p, const char *name, int *rowindex),
    mpf_QSget_named_x (mpf_QSprob p, const char *colname, mpf_t * val),
    mpf_QSget_named_rc (mpf_QSprob p, const char *colname, mpf_t * val),
    mpf_QSget_named_pi (mpf_QSprob p, const char *rowname, mpf_t * val),
    mpf_QSget_named_slack (mpf_QSprob p, const char *rowname, mpf_t * val),
    mpf_QScompute_row_norms (mpf_QSprob p),
    mpf_QSwrite_prob (mpf_QSprob p, const char *filename,
        const char *filetype),
    mpf_QSwrite_prob_file (mpf_QSprob p, FILE * file, const char *filetype),
    mpf_QSwrite_basis (mpf_QSprob p, mpf_QSbas B, const char *filename),
    mpf_QStest_row_norms (mpf_QSprob p),
    mpf_QSset_param (mpf_QSprob p, int whichparam, int newvalue),
    mpf_QSset_param_EGlpNum (mpf_QSprob p, int whichparam, mpf_t newvalue),
    mpf_QSget_param (mpf_QSprob p, int whichparam, int *value),
    mpf_QSget_param_EGlpNum (mpf_QSprob p, int whichparam, mpf_t * value);

    mpf_QSLIB_INTERFACE char *mpf_QSget_probname (mpf_QSprob p);
    mpf_QSLIB_INTERFACE char *mpf_QSget_objname (mpf_QSprob p);
    mpf_QSLIB_INTERFACE char *mpf_QSversion (void);

    mpf_QSLIB_INTERFACE mpf_QSprob mpf_QScreate_prob (const char *name,
        int objsense),
    mpf_QSread_prob (const char *filename, const char *filetype),
    mpf_QSload_prob (const char *probname, int ncols, int nrows, int *cmatcnt,
        int *cmatbeg, int *cmatind, mpf_t * cmatval, int objsense,
        mpf_t * obj, mpf_t * rhs, char *sense, mpf_t * lower, mpf_t * upper,
        const char **colnames, const char **rownames),
    mpf_QScopy_prob (mpf_QSprob p, const char *newname);

    mpf_QSLIB_INTERFACE mpf_QSbas mpf_QSget_basis (mpf_QSprob p),
    mpf_QSread_basis (mpf_QSprob p, const char *filename);

#ifdef  __cplusplus
}
#endif

/****************************************************************************
 *
 * This is the undocumented part of the QSlib interface 
 *
 ****************************************************************************/

#ifdef WIN32
typedef struct mpf_QSLIB_INTERFACE mpf_qsline_reader *mpf_QSline_reader;
typedef struct mpf_QSLIB_INTERFACE mpf_qsformat_error *mpf_QSformat_error;
typedef struct mpf_QSLIB_INTERFACE mpf_qserror_collector *mpf_QSerror_collector;
typedef struct mpf_QSLIB_INTERFACE mpf_qserror_memory *mpf_QSerror_memory;
#else
typedef struct mpf_qsline_reader *mpf_QSline_reader;
typedef struct mpf_qsformat_error *mpf_QSformat_error;
typedef struct mpf_qserror_collector *mpf_QSerror_collector;
typedef struct mpf_qserror_memory *mpf_QSerror_memory;
#endif

#ifdef  __cplusplus
extern "C"
{
#endif
    mpf_QSLIB_INTERFACE const char *mpf_QSformat_error_type_string (int tp);

    mpf_QSLIB_INTERFACE int mpf_QSerror_get_type (mpf_QSformat_error error);
    mpf_QSLIB_INTERFACE const char *mpf_QSerror_get_desc (mpf_QSformat_error error);
    mpf_QSLIB_INTERFACE int mpf_QSerror_get_line_number (mpf_QSformat_error error);
    mpf_QSLIB_INTERFACE int mpf_QSerror_get_pos (mpf_QSformat_error error);
    mpf_QSLIB_INTERFACE const char *mpf_QSerror_get_line (mpf_QSformat_error error);
    mpf_QSLIB_INTERFACE void mpf_QSerror_print (FILE * f,
        mpf_QSformat_error error);

    mpf_QSLIB_INTERFACE mpf_QSerror_collector mpf_QSerror_collector_new (void *fct,
         void *dest);
    mpf_QSLIB_INTERFACE
    mpf_QSerror_collector mpf_QSerror_memory_collector_new (mpf_QSerror_memory mem);
    mpf_QSLIB_INTERFACE void mpf_QSerror_collector_free (mpf_QSerror_collector c);

    mpf_QSLIB_INTERFACE mpf_QSline_reader mpf_QSline_reader_new (void *fct,
         void *data_src);

    mpf_QSLIB_INTERFACE void mpf_QSline_reader_free (mpf_QSline_reader reader);

    mpf_QSLIB_INTERFACE void mpf_QSline_reader_set_error_collector (mpf_QSline_reader reader,
        mpf_QSerror_collector collector);

    mpf_QSLIB_INTERFACE char *mpf_QSline_reader_get (mpf_QSline_reader reader,
         char *s, int size);

    mpf_QSLIB_INTERFACE mpf_QSprob mpf_QSget_prob (mpf_QSline_reader reader,
        const char *probname, const char *filetype);

    mpf_QSLIB_INTERFACE mpf_QSerror_memory mpf_QSerror_memory_create (int takeErrorLines);
    mpf_QSLIB_INTERFACE void mpf_QSerror_memory_free (mpf_QSerror_memory mem);

    mpf_QSLIB_INTERFACE int mpf_QSerror_memory_get_nof (mpf_QSerror_memory mem,
        int error_type);
    mpf_QSLIB_INTERFACE int mpf_QSerror_memory_get_nerrors (mpf_QSerror_memory mem);

    mpf_QSLIB_INTERFACE mpf_QSformat_error
        mpf_QSerror_memory_get_last_error (mpf_QSerror_memory mem);
    mpf_QSLIB_INTERFACE mpf_QSformat_error
        mpf_QSerror_memory_get_prev_error (mpf_QSformat_error e);

    mpf_QSLIB_INTERFACE void mpf_QSset_reporter (mpf_QSprob prob,
        int iterskip, void *fct, void *dest);

    mpf_QSLIB_INTERFACE int mpf_QSreport_prob (mpf_QSprob p,
        const char *filetype, mpf_QSerror_collector c);

#ifdef  __cplusplus
}
#endif
#endif  /* mpf___QS_QSOPT_H */

/*****************************************************************************/
/*****************************************************************************/

#ifndef __EXACT_H__
#define __EXACT_H__

/* ========================================================================= */
/* Here we define an interface to solve LP's (QSexact_solver) */
/* ========================================================================= */

/* ========================================================================= */
/* If enabled, save the last problem proved to be optimal, and its solution. */
#define QSEXACT_SAVE_OPTIMAL 0

/* ========================================================================= */
/** If enabled, save the intermediate problems created by the functions      */
/*  QScopy_prob_mpq_dbl and QScopy_prob_mpq_mpf                              */
#define QSEXACT_SAVE_INT 0

/* ========================================================================= */
/* Copy an exact problem (mpq_QSdata) to a regular double version of the     */
/* problem (dbl_QSdata)                                                      */
dbl_QSdata *QScopy_prob_mpq_dbl (mpq_QSdata * p, const char *newname);

/* ========================================================================= */
/* Copy an exact problem (mpq_QSdata) to a regular double version of the     */
/*  problem (dbl_QSdata)                                                     */
mpf_QSdata *QScopy_prob_mpq_mpf (mpq_QSdata * p, const char *newname);

/* ========================================================================= */
/* Print into a file the optimal solution.
 * param p original problem.
 * param out_f file where to write the solution.
 * return zero on success, non-zero otherwise.
 * */
int QSexact_print_sol (mpq_QSdata * p, FILE * out_f);


/* ========================================================================= */
/* Set the number of bits to use with mpf_t type numbers and change all
 * internal constants as needed. */
#define QSexact_set_precision(precision) mpf_QSset_precision(precision)

/* ========================================================================= */
/* Given an mpq_QSdata problem, solve it exactly.
 * param x if not null, we store here the primal solution to the 
 * problem (if it exist).
 * param y if not null, we store here the dual solution to the
 * problem, 
 * param p_mpq problem to solve exactly.
 * param status pointer to the integer where we will return the status
 * of the problem, either optimal, infeasible, or unbounded (we could also 
 * return time out).
 * param simplexalgo whether to use primal or dual simplex while solving
 * to optimality the problem.
 * param basis if not null, use the given basis to start the
 * iteration of simplex, and store here the optimal basis (if found).
 * return zero on success, non-zero otherwise. */

int QSexact_solver (mpq_QSdata * p_mpq, mpq_t * const x, mpq_t * const y,
        QSbasis * const basis, int simplexalgo, int *status);

#endif  /* __EXACT_H__ */
#endif  /* __QSopt_ex__  */
