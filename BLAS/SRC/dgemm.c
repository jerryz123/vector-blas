/* dgemm.f -- translated by f2c (version 20160102).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"
#include "vec-utils.h"
/* > \brief \b DGEMM */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC) */

/*       .. Scalar Arguments .. */
/*       DOUBLE PRECISION ALPHA,BETA */
/*       INTEGER K,LDA,LDB,LDC,M,N */
/*       CHARACTER TRANSA,TRANSB */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION A(LDA,*),B(LDB,*),C(LDC,*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DGEMM  performs one of the matrix-matrix operations */
/* > */
/* >    C := alpha*op( A )*op( B ) + beta*C, */
/* > */
/* > where  op( X ) is one of */
/* > */
/* >    op( X ) = X   or   op( X ) = X**T, */
/* > */
/* > alpha and beta are scalars, and A, B and C are matrices, with op( A ) */
/* > an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] TRANSA */
/* > \verbatim */
/* >          TRANSA is CHARACTER*1 */
/* >           On entry, TRANSA specifies the form of op( A ) to be used in */
/* >           the matrix multiplication as follows: */
/* > */
/* >              TRANSA = 'N' or 'n',  op( A ) = A. */
/* > */
/* >              TRANSA = 'T' or 't',  op( A ) = A**T. */
/* > */
/* >              TRANSA = 'C' or 'c',  op( A ) = A**T. */
/* > \endverbatim */
/* > */
/* > \param[in] TRANSB */
/* > \verbatim */
/* >          TRANSB is CHARACTER*1 */
/* >           On entry, TRANSB specifies the form of op( B ) to be used in */
/* >           the matrix multiplication as follows: */
/* > */
/* >              TRANSB = 'N' or 'n',  op( B ) = B. */
/* > */
/* >              TRANSB = 'T' or 't',  op( B ) = B**T. */
/* > */
/* >              TRANSB = 'C' or 'c',  op( B ) = B**T. */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >           On entry,  M  specifies  the number  of rows  of the  matrix */
/* >           op( A )  and of the  matrix  C.  M  must  be at least  zero. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >           On entry,  N  specifies the number  of columns of the matrix */
/* >           op( B ) and the number of columns of the matrix C. N must be */
/* >           at least zero. */
/* > \endverbatim */
/* > */
/* > \param[in] K */
/* > \verbatim */
/* >          K is INTEGER */
/* >           On entry,  K  specifies  the number of columns of the matrix */
/* >           op( A ) and the number of rows of the matrix op( B ). K must */
/* >           be at least  zero. */
/* > \endverbatim */
/* > */
/* > \param[in] ALPHA */
/* > \verbatim */
/* >          ALPHA is DOUBLE PRECISION. */
/* >           On entry, ALPHA specifies the scalar alpha. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension ( LDA, ka ), where ka is */
/* >           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise. */
/* >           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k */
/* >           part of the array  A  must contain the matrix  A,  otherwise */
/* >           the leading  k by m  part of the array  A  must contain  the */
/* >           matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >           On entry, LDA specifies the first dimension of A as declared */
/* >           in the calling (sub) program. When  TRANSA = 'N' or 'n' then */
/* >           LDA must be at least  max( 1, m ), otherwise  LDA must be at */
/* >           least  max( 1, k ). */
/* > \endverbatim */
/* > */
/* > \param[in] B */
/* > \verbatim */
/* >          B is DOUBLE PRECISION array, dimension ( LDB, kb ), where kb is */
/* >           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise. */
/* >           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n */
/* >           part of the array  B  must contain the matrix  B,  otherwise */
/* >           the leading  n by k  part of the array  B  must contain  the */
/* >           matrix B. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >           On entry, LDB specifies the first dimension of B as declared */
/* >           in the calling (sub) program. When  TRANSB = 'N' or 'n' then */
/* >           LDB must be at least  max( 1, k ), otherwise  LDB must be at */
/* >           least  max( 1, n ). */
/* > \endverbatim */
/* > */
/* > \param[in] BETA */
/* > \verbatim */
/* >          BETA is DOUBLE PRECISION. */
/* >           On entry,  BETA  specifies the scalar  beta.  When  BETA  is */
/* >           supplied as zero then C need not be set on input. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* >          C is DOUBLE PRECISION array, dimension ( LDC, N ) */
/* >           Before entry, the leading  m by n  part of the array  C must */
/* >           contain the matrix  C,  except when  beta  is zero, in which */
/* >           case C need not be set on entry. */
/* >           On exit, the array  C  is overwritten by the  m by n  matrix */
/* >           ( alpha*op( A )*op( B ) + beta*C ). */
/* > \endverbatim */
/* > */
/* > \param[in] LDC */
/* > \verbatim */
/* >          LDC is INTEGER */
/* >           On entry, LDC specifies the first dimension of C as declared */
/* >           in  the  calling  (sub)  program.   LDC  must  be  at  least */
/* >           max( 1, m ). */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup double_blas_level3 */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  Level 3 Blas routine. */
/* > */
/* >  -- Written on 8-February-1989. */
/* >     Jack Dongarra, Argonne National Laboratory. */
/* >     Iain Duff, AERE Harwell. */
/* >     Jeremy Du Croz, Numerical Algorithms Group Ltd. */
/* >     Sven Hammarling, Numerical Algorithms Group Ltd. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dgemm_(char *transa, char *transb, integer *m, integer *
	n, integer *k, doublereal *alpha, doublereal *a, integer *lda,
	doublereal *b, integer *ldb, doublereal *beta, doublereal *c__,
	integer *ldc, ftnlen transa_len, ftnlen transb_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2,
	    i__3;

    /* Local variables */
    static integer i__, j, l, info;
    static logical nota, notb;
    static doublereal temp;
    static integer ncola;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer nrowa, nrowb;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);


    /*  -- Reference BLAS level3 routine (version 3.7.0) -- */
    /*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    -- */
    /*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /*     December 2016 */

    /*     .. Scalar Arguments .. */
    /*     .. */
    /*     .. Array Arguments .. */
    /*     .. */

    /*  ===================================================================== */

    /*     .. External Functions .. */
    /*     .. */
    /*     .. External Subroutines .. */
    /*     .. */
    /*     .. Intrinsic Functions .. */
    /*     .. */
    /*     .. Local Scalars .. */
    /*     .. */
    /*     .. Parameters .. */
    /*     .. */

    /*     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not */
    /*     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows */
    /*     and  columns of  A  and the  number of  rows  of  B  respectively. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;

    /* Function Body */
    nota = lsame_(transa, "N", (ftnlen)1, (ftnlen)1);
    notb = lsame_(transb, "N", (ftnlen)1, (ftnlen)1);
    if (nota) {
      nrowa = *m;
      ncola = *k;
    } else {
      nrowa = *k;
      ncola = *m;
    }
    if (notb) {
      nrowb = *k;
    } else {
      nrowb = *n;
    }

    /*     Test the input parameters. */

    info = 0;
    if (! nota && ! lsame_(transa, "C", (ftnlen)1, (ftnlen)1) && ! lsame_(
                                                                          transa, "T", (ftnlen)1, (ftnlen)1)) {
      info = 1;
    } else if (! notb && ! lsame_(transb, "C", (ftnlen)1, (ftnlen)1) && !
               lsame_(transb, "T", (ftnlen)1, (ftnlen)1)) {
      info = 2;
    } else if (*m < 0) {
      info = 3;
    } else if (*n < 0) {
      info = 4;
    } else if (*k < 0) {
      info = 5;
    } else if (*lda < max(1,nrowa)) {
      info = 8;
    } else if (*ldb < max(1,nrowb)) {
      info = 10;
    } else if (*ldc < max(1,*m)) {
      info = 13;
    }
    if (info != 0) {
      xerbla_("DGEMM ", &info, (ftnlen)6);
      return 0;
    }

    if (*m == 0 || *n == 0 || (*alpha == 0. || *k == 0) && *beta == 1.) {
      return 0;
    }

    if (*alpha == 0.) {
      if (*beta == 0.) {
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
          i__2 = *m;
          for (i__ = 1; i__ <= i__2; ++i__) {
            c__[i__ + j * c_dim1] = 0.;
          }
        }
      } else {
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
          i__2 = *m;
          for (i__ = 1; i__ <= i__2; ++i__) {
            c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
          }
        }
      }
      return 0;
    }
    if (notb) {
      if (nota) {
        setvcfg0(VFP64,VFP64,VFP64,VFP64);
        setvcfg2(VFP64,SFP64,SFP64,SFP64);
        setvcfg4(SFP64,SFP64,SFP64,SFP64);
        int vl = 0;
        double* ca;
        double* cb;
        i__ = 0;
        i__1 = *n;
        i__2 = *k;
        i__3 = *m;
        ca = a;
        cb = b;
        asm volatile ("vinsert v5, %0, x0" : : "r" (*alpha));
        asm volatile ("vinsert v6, %0, x0" : : "r" (*beta));
        while (i__3 - i__ > 0) {
          setvl(vl, i__3 - i__);
          j = 1;
          for (j = 1; j <= i__1 - 4; j += 4) {
            asm volatile ("vsne v0, v0, v0");
            asm volatile ("vsne v1, v0, v0");
            asm volatile ("vsne v2, v0, v0");
            asm volatile ("vsne v3, v0, v0");
            ca = a + i__ + 1 + a_dim1;
            cb = b + 1 + j * b_dim1;
            for (l = 1; l <= i__2; ++l) {
              asm volatile ("vld v4, 0(%0)" : : "r" (ca));
              asm volatile ("vinsert v8, %0, x0" : : "r" (cb[0]));
              asm volatile ("vinsert v9, %0, x0" : : "r" (cb[b_dim1]));
              asm volatile ("vinsert v10, %0, x0" : : "r" (cb[b_dim1*2]));
              asm volatile ("vinsert v11, %0, x0" : : "r" (cb[b_dim1*3]));
              asm volatile ("vmadd v0, v4, v8,  v0");
              asm volatile ("vmadd v1, v4, v9,  v1");
              asm volatile ("vmadd v2, v4, v10, v2");
              asm volatile ("vmadd v3, v4, v11, v3");
              ca += a_dim1;
              cb += 1;
            }
            asm volatile ("vmul v0, v0, v5");
            asm volatile ("vmul v1, v1, v5");
            asm volatile ("vmul v2, v2, v5");
            asm volatile ("vmul v3, v3, v5");
            if (*beta != 0.f) {
              asm volatile ("vld v4, 0(%0)" : : "r" (&c__[i__ + 1 + j * c_dim1]));
              asm volatile ("vmadd v0, v4, v6, v0");
              asm volatile ("vld v4, 0(%0)" : : "r" (&c__[i__ + 1 + j * c_dim1 + c_dim1]));
              asm volatile ("vmadd v1, v4, v6, v1");
              asm volatile ("vld v4, 0(%0)" : : "r" (&c__[i__ + 1 + j * c_dim1 + c_dim1 * 2]));
              asm volatile ("vmadd v2, v4, v6, v2");
              asm volatile ("vld v4, 0(%0)" : : "r" (&c__[i__ + 1 + j * c_dim1 + c_dim1 * 3]));
              asm volatile ("vmadd v3, v4, v6, v3");
            }
            asm volatile ("vst v0, 0(%0)" : : "r" (&c__[i__ + 1 + j * c_dim1]));
            asm volatile ("vst v1, 0(%0)" : : "r" (&c__[i__ + 1 + j * c_dim1 + c_dim1]));
            asm volatile ("vst v2, 0(%0)" : : "r" (&c__[i__ + 1 + j * c_dim1 + c_dim1 * 2]));
            asm volatile ("vst v3, 0(%0)" : : "r" (&c__[i__ + 1 + j * c_dim1 + c_dim1 * 3]));
          }
          for (; j <= i__1; ++j) {
            asm volatile ("vsne v0, v0, v0");
            ca = a + i__ + 1 + a_dim1;
            cb = b + 1 + j * b_dim1;
            for (l = 1; l <= i__2; ++l) {
              asm volatile ("vld v4, 0(%0)" : : "r" (ca));
              asm volatile ("vinsert v8, %0, x0" : : "r" (*cb));
              asm volatile ("vmadd v0, v4, v8, v0");
              ca += a_dim1;
              cb += 1;
            }
            asm volatile ("vmul v0, v0, v5");
            if (*beta != 0.f) {
              asm volatile ("vld v1, 0(%0)" : : "r" (&c__[i__ + 1 + j * c_dim1]));
              asm volatile ("vmadd v0, v1, v6, v0");
            }
            asm volatile ("vst v0, 0(%0)" : : "r" (&c__[i__ + 1 + j * c_dim1]));
          }
            i__ += vl;
        }
      } else {
        double* ta = malloc(*k * *m * sizeof(double));

        for (i__ = 1; i__ <= *m; i__++)
          for (l = 1; l <= *k; l++)
            {
              ta[i__ - 1 + (l - 1) * *m] = a[l + i__ * a_dim1];
            }
        dgemm_("N", "N", m, n, k, alpha, ta, m, b + b_offset, ldb, beta, c__ + c_offset, ldc, 1, 1);
        free(ta);
        return 0;
      }
    } else {
      if (nota) {
        setvcfg0(VFP64,VFP64,VFP64,VFP64);
        setvcfg2(VFP64,SFP64,SFP64,SFP64);
        setvcfg4(SFP64,SFP64,SFP64,SFP64);
/*           Form  C := alpha*A*B + beta*C. */
          int vl = 0;
          double* ca;
          double* cb;
          i__ = 0;
          i__1 = *n;
          i__2 = *k;
          i__3 = *m;
          ca = a;
          cb = b;
          asm volatile ("vinsert v5, %0, x0" : : "r" (*alpha));
          asm volatile ("vinsert v6, %0, x0" : : "r" (*beta));
          while (i__3 - i__ > 0) {
            setvl(vl, i__3 - i__);
            j = 1;
             for (j = 1; j <= i__1 - 4; j += 4) {
              asm volatile ("vsne v0, v0, v0");
              asm volatile ("vsne v1, v0, v0");
              asm volatile ("vsne v2, v0, v0");
              asm volatile ("vsne v3, v0, v0");
              ca = a + i__ + 1 + a_dim1;
              cb = b + j + b_dim1;
              for (l = 1; l <= i__2; ++l) {
                asm volatile ("vld v4, 0(%0)" : : "r" (ca));
                asm volatile ("vinsert v8, %0, x0"  : : "r" (cb[0]));
                asm volatile ("vinsert v9, %0, x0"  : : "r" (cb[1]));
                asm volatile ("vinsert v10, %0, x0" : : "r" (cb[2]));
                asm volatile ("vinsert v11, %0, x0" : : "r" (cb[3]));
                asm volatile ("vmadd v0, v4, v8,  v0");
                asm volatile ("vmadd v1, v4, v9,  v1");
                asm volatile ("vmadd v2, v4, v10, v2");
                asm volatile ("vmadd v3, v4, v11, v3");
                ca += a_dim1;
                cb += b_dim1;
              }
              asm volatile ("vmul v0, v0, v5");
              asm volatile ("vmul v1, v1, v5");
              asm volatile ("vmul v2, v2, v5");
              asm volatile ("vmul v3, v3, v5");
              if (*beta != 0.f) {
                asm volatile ("vld v4, 0(%0)" : : "r" (&c__[i__ + 1 + j * c_dim1]));
                asm volatile ("vmadd v0, v4, v6, v0");
                asm volatile ("vld v4, 0(%0)" : : "r" (&c__[i__ + 1 + j * c_dim1 + c_dim1]));
                asm volatile ("vmadd v1, v4, v6, v1");
                asm volatile ("vld v4, 0(%0)" : : "r" (&c__[i__ + 1 + j * c_dim1 + c_dim1 * 2]));
                asm volatile ("vmadd v2, v4, v6, v2");
                asm volatile ("vld v4, 0(%0)" : : "r" (&c__[i__ + 1 + j * c_dim1 + c_dim1 * 3]));
                asm volatile ("vmadd v3, v4, v6, v3");

              }
              asm volatile ("vst v0, 0(%0)" : : "r" (&c__[i__ + 1 + j * c_dim1]));
              asm volatile ("vst v1, 0(%0)" : : "r" (&c__[i__ + 1 + j * c_dim1 + c_dim1]));
              asm volatile ("vst v2, 0(%0)" : : "r" (&c__[i__ + 1 + j * c_dim1 + c_dim1 * 2]));
              asm volatile ("vst v3, 0(%0)" : : "r" (&c__[i__ + 1 + j * c_dim1 + c_dim1 * 3]));
            }
            for (; j <= i__1; ++j) {
              asm volatile ("vsne v0, v0, v0");
              ca = a + i__ + 1 + a_dim1;
              cb = b + j + b_dim1;
              for (l = 1; l <= i__2; ++l) {
                asm volatile ("vld v4, 0(%0)" : : "r" (ca));
                asm volatile ("vinsert v8, %0, x0" : : "r" (cb[0]));
                asm volatile ("vmadd v0, v4, v8, v0");
                ca += a_dim1;
                cb += b_dim1;
              }
              asm volatile ("vmul v0, v0, v5");
              if (*beta != 0.f) {
                asm volatile ("vld v1, 0(%0)" : : "r" (&c__[i__ + 1 + j * c_dim1]));
                asm volatile ("vmadd v0, v1, v6, v0");
              }
              asm volatile ("vst v0, 0(%0)" : : "r" (&c__[i__ + 1 + j * c_dim1]));
            }
            i__ += vl;
          }
      } else {
        double* ta = malloc(*k * *m * sizeof(double));

        for (i__ = 1; i__ <= *m; i__++)
          for (l = 1; l <= *k; l++)
            {
            ta[i__ - 1 + (l - 1) * *m] = a[l + i__ * a_dim1];
            }
        dgemm_("N", "T", m, n, k, alpha, ta, m, b + b_offset, ldb, beta, c__ + c_offset, ldc, 1, 1);
        free(ta);
        return 0;
      }
    }
    return 0;
} /* dgemm_ */
