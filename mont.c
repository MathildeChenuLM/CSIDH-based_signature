
#include <assert.h>

#include "mont.h"

void xDBLADD(proj *R, proj *S, proj const *P, proj const *Q, proj const *PQ, proj const *A)
{
    fp a, b, c, d;

    fp_add3(&a, &Q->x, &Q->z);
    fp_sub3(&b, &Q->x, &Q->z);
    fp_add3(&c, &P->x, &P->z);
    fp_sub3(&d, &P->x, &P->z);
    fp_sq2(&R->x, &c);
    fp_sq2(&S->x, &d);
    fp_mul2(&c, &b);
    fp_mul2(&d, &a);
    fp_sub3(&b, &R->x, &S->x);
    fp_add3(&a, &A->z, &A->z); /* multiplication by 2 */
    fp_mul3(&R->z, &a, &S->x);
    fp_add3(&S->x, &A->x, &a);
    fp_add2(&R->z, &R->z); /* multiplication by 2 */
    fp_mul2(&R->x, &R->z);
    fp_mul2(&S->x, &b);
    fp_sub3(&S->z, &c, &d);
    fp_add2(&R->z, &S->x);
    fp_add3(&S->x, &c, &d);
    fp_mul2(&R->z, &b);
    fp_sq2(&d, &S->z);
    fp_sq2(&b, &S->x);
    fp_mul3(&S->x, &PQ->z, &b);
    fp_mul3(&S->z, &PQ->x, &d);
}


void xDBL(proj *Q, proj const *A, proj const *P)
{
    fp a, b, c;
    fp_add3(&a, &P->x, &P->z);
    fp_sq1(&a);
    fp_sub3(&b, &P->x, &P->z);
    fp_sq1(&b);
    fp_sub3(&c, &a, &b);
    fp_add2(&b, &b); fp_add2(&b, &b); /* multiplication by 4 */
    fp_mul2(&b, &A->z);
    fp_mul3(&Q->x, &a, &b);
    fp_add3(&a, &A->z, &A->z); /* multiplication by 2 */
    fp_add2(&a, &A->x);
    fp_mul2(&a, &c);
    fp_add2(&a, &b);
    fp_mul3(&Q->z, &a, &c);
}//4M + 2S + 8A


void xADD(proj *S, proj const *P, proj const *Q, proj const *PQ)
{
    fp a, b, c, d;
    fp_add3(&a, &P->x, &P->z);
    fp_sub3(&b, &P->x, &P->z);
    fp_add3(&c, &Q->x, &Q->z);
    fp_sub3(&d, &Q->x, &Q->z);
    fp_mul2(&a, &d);
    fp_mul2(&b, &c);
    fp_add3(&c, &a, &b);
    fp_sub3(&d, &a, &b);
    fp_sq1(&c);
    fp_sq1(&d);
//     fp_mul3(&S->x, &PQ->z, &c);
    fp_mul3(&a, &PQ->z, &c);
    fp_mul3(&S->z, &PQ->x, &d);
    S->x = a;
}//Cost 4M + 2S + 6A



/* Montgomery ladder. */
/* P must not be the unique point of order 2. */
/* not constant-time! */
void xMUL(proj *Q, proj const *A, proj const *P, u512 const *k)
{
    proj R = *P;
    const proj Pcopy = *P; /* in case Q = P */

    Q->x = fp_1;
    Q->z = fp_0;

    
    
    unsigned long i = 512;
    while (--i && !u512_bit(k, i));

    do {

        bool bit = u512_bit(k, i);

        if (bit) { proj T = *Q; *Q = R; R = T; } /* not constant-time */
        //fp_cswap(&Q->x, &R.x, bit);
        //fp_cswap(&Q->z, &R.z, bit);

        xDBLADD(Q, &R, Q, &R, &Pcopy, A);

        if (bit) { proj T = *Q; *Q = R; R = T; } /* not constant-time */
        //fp_cswap(&Q->x, &R.x, bit);
        //fp_cswap(&Q->z, &R.z, bit);

    } while (i--);
}

/* computes the isogeny with kernel point K of order k */
/* returns the new curve coefficient A and the image of P */
/* (obviously) not constant time in k */
void xISOG(proj *A, proj *P, proj const *K, uint64_t k)
{
    assert (k >= 3);
    assert (k % 2 == 1);

    fp tmp0, tmp1;
    fp T[4] = {K->z, K->x, K->x, K->z};
    proj Q;
// For Evaluation
    fp_mul3(&Q.x,  &P->x, &K->x);
    fp_mul3(&tmp0, &P->z, &K->z);
    fp_sub2(&Q.x,  &tmp0);

    fp_mul3(&Q.z,  &P->x, &K->z);
    fp_mul3(&tmp0, &P->z, &K->x);
    fp_sub2(&Q.z,  &tmp0);
//=======
    proj M[3] = {*K};
    xDBL(&M[1], A, K);

    for (uint64_t i = 1; i < k / 2; ++i) {

        if (i >= 2)
            xADD(&M[i % 3], &M[(i - 1) % 3], K, &M[(i - 2) % 3]);

        fp_mul3(&tmp0, &M[i % 3].x, &T[0]);
        fp_mul3(&tmp1, &M[i % 3].z, &T[1]);
        fp_add3(&T[0], &tmp0, &tmp1);

        fp_mul2(&T[1], &M[i % 3].x);

        fp_mul3(&tmp0, &M[i % 3].z, &T[2]);
        fp_mul3(&tmp1, &M[i % 3].x, &T[3]);
        fp_add3(&T[2], &tmp0, &tmp1);

        fp_mul2(&T[3], &M[i % 3].z);

// //         For evaluation
        fp_mul3(&tmp0, &P->x, &M[i % 3].x);
        fp_mul3(&tmp1, &P->z, &M[i % 3].z);
        fp_sub2(&tmp0, &tmp1);
        fp_mul2(&Q.x,  &tmp0);

        fp_mul3(&tmp0, &P->x, &M[i % 3].z);
        fp_mul3(&tmp1, &P->z, &M[i % 3].x);
        fp_sub2(&tmp0, &tmp1);
        fp_mul2(&Q.z,  &tmp0);
    }

    fp_mul2(&T[0], &T[1]);
    fp_add2(&T[0], &T[0]); /* multiplication by 2 */

    fp_sq1(&T[1]);

    fp_mul2(&T[2], &T[3]);
    fp_add2(&T[2], &T[2]); /* multiplication by 2 */

    fp_sq1(&T[3]);

    /* Ax := T[1] * T[3] * Ax - 3 * Az * (T[1] * T[2] - T[0] * T[3]) */
    fp_mul3(&tmp0, &T[1], &T[2]);
    fp_mul3(&tmp1, &T[0], &T[3]);
    fp_sub2(&tmp0, &tmp1);
    fp_mul2(&tmp0, &A->z);
    fp_add3(&tmp1, &tmp0, &tmp0); fp_add2(&tmp0, &tmp1); /* multiplication by 3 */

    fp_mul3(&tmp1, &T[1], &T[3]);
    fp_mul2(&tmp1, &A->x);

    fp_sub3(&A->x, &tmp1, &tmp0);

    /* Az := Az * T[3]^2 */
    fp_sq1(&T[3]);
    fp_mul2(&A->z, &T[3]);

    /* X := X * Xim^2, Z := Z * Zim^2 */
    fp_sq1(&Q.x);
    fp_sq1(&Q.z);
    fp_mul2(&P->x, &Q.x);
    fp_mul2(&P->z, &Q.z);
}

//simultaneous square-and-multiply, computes x^exp and y^exp 
void exp_by_squaring_(fp* x, fp* y, uint64_t exp)
{
    fp result1, result2;
    fp_set(&result1, 1);
    fp_set(&result2, 1);

    while (exp)
    {
        if (exp & 1){
          fp_mul2(&result1, x);
          fp_mul2(&result2, y);
    }
    
        fp_sq1(x);
    fp_sq1(y);
        exp >>= 1;
    }

    fp_cswap(&result1, x, 1);
    fp_cswap(&result2, y, 1);

}

void xISOG_MR(proj *A, proj *P, proj const *K, uint64_t k)
{
    assert (k >= 3);
    assert (k % 2 == 1);

    fp tmp0, tmp1, tmp2, Psum, Pdif;
    proj Q, Aed, prod;

    fp_add3(&Aed.z, &A->z, &A->z);  //compute twisted Edwards curve coefficients
    fp_add3(&Aed.x, &A->x, &Aed.z);
    fp_sub3(&Aed.z, &A->x, &Aed.z);
   
    fp_add3(&Psum, &P->x, &P->z);   //precomputations
    fp_sub3(&Pdif, &P->x, &P->z);

    fp_sub3(&prod.x, &K->x, &K->z);
    fp_add3(&prod.z, &K->x, &K->z);
    
    fp_mul3(&tmp1, &prod.x, &Psum);
    fp_mul3(&tmp0, &prod.z, &Pdif);
    fp_add3(&Q.x, &tmp0, &tmp1);
    fp_sub3(&Q.z, &tmp0, &tmp1);

    proj M[3] = {*K};
    xDBL(&M[1], A, K);

    for (uint64_t i = 1; i < k / 2; ++i) {

        if (i >= 2)
            xADD(&M[i % 3], &M[(i - 1) % 3], K, &M[(i - 2) % 3]);

    fp_sub3(&tmp1, &M[i % 3].x, &M[i % 3].z);
        fp_add3(&tmp0, &M[i % 3].x, &M[i % 3].z);
    fp_mul2(&prod.x, &tmp1);
        fp_mul2(&prod.z, &tmp0);
        fp_mul2(&tmp1, &Psum);
        fp_mul2(&tmp0, &Pdif);
        fp_add3(&tmp2, &tmp0, &tmp1);
    fp_mul2(&Q.x, &tmp2);
        fp_sub3(&tmp2, &tmp0, &tmp1);
    fp_mul2(&Q.z, &tmp2);

    }


    // point evaluation
    fp_sq1(&Q.x);
    fp_sq1(&Q.z);
    fp_mul2(&P->x, &Q.x);
    fp_mul2(&P->z, &Q.z);

    //compute Aed.x^k, Aed.z^k
    exp_by_squaring_(&Aed.x, &Aed.z, k);

    //compute prod.x^8, prod.z^8
    fp_sq1(&prod.x);
    fp_sq1(&prod.x);
    fp_sq1(&prod.x);
    fp_sq1(&prod.z);
    fp_sq1(&prod.z);
    fp_sq1(&prod.z);

    //compute image curve parameters
    fp_mul2(&Aed.z, &prod.x);
    fp_mul2(&Aed.x, &prod.z);

    //compute Montgomery params
    fp_add3(&A->x, &Aed.x, &Aed.z);
    fp_sub3(&A->z, &Aed.x, &Aed.z);
    fp_add2(&A->x, &A->x);
}

