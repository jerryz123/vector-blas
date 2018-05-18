/* caxpy.f -- translated by f2c (version 20160102).
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

/* > \brief \b CAXPY */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CAXPY(N,CA,CX,INCX,CY,INCY) */

/*       .. Scalar Arguments .. */
/*       COMPLEX CA */
/*       INTEGER INCX,INCY,N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX CX(*),CY(*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >    CAXPY constant times a vector plus a vector. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >         number of elements in input vector(s) */
/* > \endverbatim */
/* > */
/* > \param[in] CA */
/* > \verbatim */
/* >          CA is COMPLEX */
/* >           On entry, CA specifies the scalar alpha. */
/* > \endverbatim */
/* > */
/* > \param[in] CX */
/* > \verbatim */
/* >          CX is COMPLEX array, dimension ( 1 + ( N - 1 )*abs( INCX ) ) */
/* > \endverbatim */
/* > */
/* > \param[in] INCX */
/* > \verbatim */
/* >          INCX is INTEGER */
/* >         storage spacing between elements of CX */
/* > \endverbatim */
/* > */
/* > \param[in,out] CY */
/* > \verbatim */
/* >          CY is COMPLEX array, dimension ( 1 + ( N - 1 )*abs( INCY ) ) */
/* > \endverbatim */
/* > */
/* > \param[in] INCY */
/* > \verbatim */
/* >          INCY is INTEGER */
/* >         storage spacing between elements of CY */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2017 */

/* > \ingroup complex_blas_level1 */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >     jack dongarra, linpack, 3/11/78. */
/* >     modified 12/3/93, array(1) declarations changed to array(*) */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int caxpy_(integer *n, complex *ca, complex *cx, integer *
	incx, complex *cy, integer *incy)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    complex q__1, q__2;

    /* Local variables */
    static integer i__, ix, iy;
    extern doublereal scabs1_(complex *);


/*  -- Reference BLAS level1 routine (version 3.8.0) -- */
/*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2017 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
    /* Parameter adjustments */
    //--cy;
    //--cx;

    /* Function Body */
    if (*n <= 0)
      return 0;
    if (scabs1_(ca) == 0.f)
      return 0;

    setvcfg0(VFP32, // y[] real
             VFP32, // x[] real
             SFP32, // da  real
             VFP32);// y[] imag
    setvcfg2(VFP32, // x[] imag
             SFP32, // da  imag
             SFP32, //
             SFP32);
    int vl = 1;
    asm volatile ("vinsert v2, %0, x0" : : "r" (ca->r));
    asm volatile ("vinsert v5, %0, x0" : : "r" (ca->i));
    i__ = 0;

    if (*incx < 0)
      cx += (-(*n) + 1) * *incx;

    if (*incy < 0)
      cy += (-(*n) + 1) * *incy;

    int sincx = *incx << 3;
    int sincy = *incy << 3;
    int dincx = *incx;
    int dincy = *incy;
    int dn = *n;
    while (i__ < dn)
      {
        setvl(vl, dn - i__);
        asm volatile ("vlds  v0, 0(%0), %1" : : "r" (cx), "r" (sincx));
        asm volatile ("vlds  v1, 0(%0), %1" : : "r" (cy), "r" (sincy));
        asm volatile ("vlds  v4, 4(%0), %1" : : "r" (cx), "r" (sincx));
        asm volatile ("vlds  v3, 4(%0), %1" : : "r" (cy), "r" (sincy));

        asm volatile ("vmadd  v1, v2, v0, v1");
        asm volatile ("vnmsub v1, v5, v4, v1");
        asm volatile ("vmadd  v3, v2, v4, v3");
        asm volatile ("vmadd  v3, v5, v0, v3");
        asm volatile ("vsts  v1, 0(%0), %1" : : "r" (cy), "r" (sincy));
        asm volatile ("vsts  v3, 4(%0), %1" : : "r" (cy), "r" (sincy));

        i__ += vl;
        cx += vl * dincx;
        cy += vl * dincy;
      }

    return 0;
} /* caxpy_ */
