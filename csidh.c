
#include <string.h>
#include <assert.h>
#include <stdio.h>

#include "csidh.h"
#include "rng.h"

/* specific to p, should perhaps be somewhere else */
const unsigned primes[num_primes] = {
      3,   5,   7,  11,  13,  17,  19,  23,  29,  31,  37,  41,  43,  47,  53,  59,
     61,  67,  71,  73,  79,  83,  89,  97, 101, 103, 107, 109, 113, 127, 131, 137,
    139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227,
    229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313,
    317, 331, 337, 347, 349, 353, 359, 367, 373, 587,
};

const u512 four_sqrt_p = {{
    0x85e2579c786882cf, 0x4e3433657e18da95, 0x850ae5507965a0b3, 0xa15bc4e676475964,
}};
// why do we need square roots ?


const public_key base = {0}; /* A = 0 */


void csidh_private(private_key *priv)
{
    /* Private key generation.
    */
    memset(&priv->e, 0, sizeof(priv->e));
    for (size_t i = 0; i < num_primes; ) {
        int8_t buf[64];
        randombytes(buf, sizeof(buf));
        for (size_t j = 0; j < sizeof(buf); ++j) {
            if (buf[j] <= max_exponent && buf[j] >= -max_exponent) {
                priv->e[i / 2] |= (buf[j] & 0xf) << i % 2 * 4;
                if (++i >= num_primes)
                    break;
            }
        }
    }
}

/* compute [(p+1)/l] P for all l in our list of primes. */
/* divide and conquer is much faster than doing it naively,
 * but uses more memory. */
static void cofactor_multiples(proj *P, const proj *A, size_t lower, size_t upper)
{
    assert(lower < upper);

    if (upper - lower == 1)
        return;

    size_t mid = lower + (upper - lower + 1) / 2;

    u512 cl = u512_1, cu = u512_1;
    for (size_t i = lower; i < mid; ++i)
        u512_mul3_64(&cu, &cu, primes[i]);
    for (size_t i = mid; i < upper; ++i)
        u512_mul3_64(&cl, &cl, primes[i]);

    xMUL(&P[mid], A, &P[lower], &cu);
    xMUL(&P[lower], A, &P[lower], &cl);

    cofactor_multiples(P, A, lower, mid);
    cofactor_multiples(P, A, mid, upper);
}

/* never accepts invalid keys. */
bool validate(public_key const *in)
{
    const proj A = {in->A, fp_1};
    do {
        proj P[num_primes];
        fp_random(&P->x);
        P->z = fp_1;
        /* maximal 2-power in p+1 */
        xDBL(P, &A, P);
        xDBL(P, &A, P);
        cofactor_multiples(P, &A, 0, num_primes);
        u512 order = u512_1;
        for (size_t i = num_primes - 1; i < num_primes; --i) {
            /* we only gain information if [(p+1)/l] P is non-zero */
            if (memcmp(&P[i].z, &fp_0, sizeof(fp))) {
                u512 tmp;
                u512_set(&tmp, primes[i]);
                xMUL(&P[i], &A, &P[i], &tmp);
                if (memcmp(&P[i].z, &fp_0, sizeof(fp)))
                    /* P does not have order dividing p+1. */
                    return false;
                u512_mul3_64(&order, &order, primes[i]);
                if (u512_sub3(&tmp, &four_sqrt_p, &order)) /* returns borrow */
                    /* order > 4 sqrt(p), hence definitely supersingular */
                    return true;
            }
        }
    /* P didn't have big enough order to prove supersingularity. */
    } while (1);
}

/* never accepts invalid keys. */
bool validate_2(public_key const *in)
{
    /* same than above, except that the multiplication is changed.
    */
    const proj A = {in->A, fp_1};
    do {
        proj P[num_primes];
        fp_random(&P->x);
        P->z = fp_1;
        /* maximal 2-power in p+1 */
        xDBL(P, &A, P);
        xDBL(P, &A, P);
        cofactor_multiples(P, &A, 0, num_primes);
        u512 order = u512_1;
        for (size_t i = num_primes - 1; i < num_primes; --i) {
            /* we only gain information if [(p+1)/l] P is non-zero */
            if (memcmp(&P[i].z, &fp_0, sizeof(fp))) {
                u512 tmp;
                u512_set(&tmp, primes[i]);
                xMUL_A24afin(&P[i], &A, &P[i], &tmp); // changed here !!
                if (memcmp(&P[i].z, &fp_0, sizeof(fp)))
                    /* P does not have order dividing p+1. */
                    return false;
                u512_mul3_64(&order, &order, primes[i]);
                if (u512_sub3(&tmp, &four_sqrt_p, &order)) /* returns borrow */
                    /* order > 4 sqrt(p), hence definitely supersingular */
                    return true;
            }
        }
    /* P didn't have big enough order to prove supersingularity. */
    } while (1);
}


/* compute x^3 + Ax^2 + x */
static void montgomery_rhs(fp *rhs, fp const *A, fp const *x)
{
    fp tmp;
    *rhs = *x;
    fp_sq1(rhs);
    fp_mul3(&tmp, A, x);
    fp_add2(rhs, &tmp);
    fp_add2(rhs, &fp_1);
    fp_mul2(rhs, x);
}//Cost 2M + 1S + 2A


/* compute x^3 + Ax^2 + x as x(x(x+A) + 1) */
static void montgomery_rhs_own(fp *rhs, fp const *A, fp const *x)
{
    //Daniel and Eduardo, better formula.
    fp_add3(rhs, x, A);
    fp_mul2(rhs, x);
    fp_add2(rhs, &fp_1);
    fp_mul2(rhs, x);
}//Cost 2M + 2A

// // compute x^3 + (A/C)x^2 + x as Cx(x(Cx+A) + c) 
static void montgomery_rhs_proj(fp *rhs, proj const *A, fp const *x)
{
    //WATCH OUT ! No usar !!!
    fp cx;
    fp_mul3(&cx, &A->z, x);//cx := c*x;
    fp_add3(rhs, &cx, &A->x);// rhs := (cx + A);
    fp_mul2(rhs, x);//rhs := (cx + A)*x;
    fp_add2(rhs, &A->z);//rhs := (cx + A)*x + C;
    fp_mul2(rhs, &cx);//rhs := cx*((cx + A)*x + C);
}//Cost 3M + 2A

/* totally not constant-time. */
void action(public_key *out, public_key const *in, private_key const *priv)
{
    u512 k[2];
    u512_set(&k[0], 4); /* maximal 2-power in p+1 */
    u512_set(&k[1], 4); /* maximal 2-power in p+1 */

    uint8_t e[2][num_primes];
    for (size_t i = 0; i < num_primes; ++i) {
        int8_t t = (int8_t) (priv->e[i / 2] << i % 2 * 4) >> 4;
        //int8_t t = (int8_t) (priv->e[i / 2] << ((i % 2) * 4)) >> 4;
        if (t > 0) {
            e[0][i] = t;
            e[1][i] = 0;
            u512_mul3_64(&k[1], &k[1], primes[i]);
        }
        else if (t < 0) {
            e[1][i] = -t;
            e[0][i] = 0;
            u512_mul3_64(&k[0], &k[0], primes[i]);
        }
        else {
            e[0][i] = 0;
            e[1][i] = 0;
            u512_mul3_64(&k[0], &k[0], primes[i]);
            u512_mul3_64(&k[1], &k[1], primes[i]);
        }
    }
    proj A = {in->A, fp_1};
    bool done[2] = {false, false};
    do {
        assert(!memcmp(&A.z, &fp_1, sizeof(fp))); // memcopy returns 0 if it could copy everything.
        proj P;
        fp_random(&P.x);
        P.z = fp_1;
        fp rhs;
        montgomery_rhs(&rhs, &A.x, &P.x);
        bool sign = !fp_issquare(&rhs);
        if (done[sign])
            continue;
        xMUL(&P, &A, &P, &k[sign]);
        done[sign] = true;
        for (size_t i = 0; i < num_primes; ++i) {
            if (e[sign][i]) {

                u512 cof = u512_1;
                for (size_t j = i + 1; j < num_primes; ++j)
                    if (e[sign][j])
                        u512_mul3_64(&cof, &cof, primes[j]);

                proj K;
                xMUL(&K, &A, &P, &cof);

                if (memcmp(&K.z, &fp_0, sizeof(fp))) {

                    xISOG(&A, &P, &K, primes[i]);

                    if (!--e[sign][i])
                        u512_mul3_64(&k[sign], &k[sign], primes[i]);

                }

            }
            done[sign] &= !e[sign][i];
        }

        fp_inv(&A.z);
        fp_mul2(&A.x, &A.z);
        A.z = fp_1;
    } while (!(done[0] && done[1]));
    out->A = A.x;
}

/* totally not constant-time. */
void action_2(public_key *out, public_key const *in, private_key const *priv)
{
    u512 k[2];
    u512_set(&k[0], 4); /* maximal 2-power in p+1 */
    u512_set(&k[1], 4); /* maximal 2-power in p+1 */
    uint8_t e[2][num_primes];
    for (size_t i = 0; i < num_primes; ++i) {
        int8_t t = (int8_t) (priv->e[i / 2] << i % 2 * 4) >> 4;
        if (t > 0) {
            e[0][i] = t;
            e[1][i] = 0;
            u512_mul3_64(&k[1], &k[1], primes[i]);
        }
        else if (t < 0) {
            e[1][i] = -t;
            e[0][i] = 0;
            u512_mul3_64(&k[0], &k[0], primes[i]);
        }
        else {
            e[0][i] = 0;
            e[1][i] = 0;
            u512_mul3_64(&k[0], &k[0], primes[i]);
            u512_mul3_64(&k[1], &k[1], primes[i]);
        }
    }
    proj A = {in->A, fp_1};
    bool done[2] = {false, false};
    do {
//         assert(!memcmp(&A.z, &fp_1, sizeof(fp)));
        proj P;
        fp_random(&P.x);
        P.z = fp_1;
        fp rhs;
//         montgomery_rhs_proj(&rhs, &A, &P.x);
        montgomery_rhs_own(&rhs, &A.x, &P.x);
//         bool sign = !fp_issquare(&rhs);
        bool sign = isSquare(&rhs);
        if (done[sign])
            continue;
        xMUL_afin(&P, &A, &P, &k[sign]);
        done[sign] = true;
        for (size_t i = 0; i < num_primes; ++i) {
            if (e[sign][i]) {
                u512 cof = u512_1;
                for (size_t j = i + 1; j < num_primes; ++j)
                    if (e[sign][j])
                        u512_mul3_64(&cof, &cof, primes[j]);
                    proj K;
                    xMUL(&K, &A, &P, &cof);
                    if (memcmp(&K.z, &fp_0, sizeof(fp))) {
                        xISOG_2(&A, &P, &K, primes[i]);
                        if (!--e[sign][i])
                            u512_mul3_64(&k[sign], &k[sign], primes[i]);
                }
            }
            done[sign] &= !e[sign][i];
        }
        fp_inv_own(&A.z, &A.z);
//         fp_inv(&A.z);
        fp_mul2(&A.x, &A.z);
        A.z = fp_1;
    } while (!(done[0] && done[1]));
    out->A = A.x;
}

void action_3(public_key *out, public_key const *in, private_key const *priv)
{
    u512 k[2];
    u512_set(&k[0], 4); /* maximal 2-power in p+1 */ // product of l_i with sign -
    u512_set(&k[1], 4); /* maximal 2-power in p+1 */ // product of l_i with sign +
    uint8_t e[2][num_primes];
    for (size_t i = 0; i < num_primes; ++i) {
        int8_t t = (int8_t) (priv->e[i / 2] << i % 2 * 4) >> 4;
        if (t > 0) {
            e[0][i] = t;
            e[1][i] = 0;
            u512_mul3_64(&k[1], &k[1], primes[i]);
        }
        else if (t < 0) {
            e[1][i] = -t;
            e[0][i] = 0;
            u512_mul3_64(&k[0], &k[0], primes[i]);
        }
        else {
            e[0][i] = 0;
            e[1][i] = 0;
            u512_mul3_64(&k[0], &k[0], primes[i]);
            u512_mul3_64(&k[1], &k[1], primes[i]);
        }
    }



    proj A = {in->A, fp_1};  // fp_1 = 1.
    bool done[2] = {false, false};

    // TORSION PART.
    
    fp Tx = {0x980BC417887630B8, 0x3223E5E8BC417B6A, 0x7239597EF4B7DEAC, 0x45AA9B1894E63200, 
    		0x45B6E360BBB5BCFD, 0x2303C882AD3A21D0, 0xF9B62E2CFCEBB476, 0x5C1DCF5FCF76F5E0};

    fp Tz = fp_1;
    proj T = {Tx, Tz};

    proj T_i;
    done[0] = true;
    for(int i=0; i<num_primes; i++){
    	u512 coef;
    	if( e[0][i] > 0 ){

    		// Setting the coef.
    		u512_set(&coef, 4);
    		for(int j=0; j<num_primes; j++){
    			if(i != j){
    				u512_mul3_64(&coef, &coef, primes[j]); // coef = (p+1) / l_i.
    			}
    		}
    			 

    		// Computing l_i torsion point.
    		xMUL(&T_i, &A, &T, &coef);

    		// Computing isogeny.
    		xISOG_2(&A, &T, &T_i, primes[i]);
    		//e[0][i] = e[0][i] - 1;
    		if (!--e[0][i])
                        u512_mul3_64(&k[0], &k[0], primes[i]);

    	}
    	done[0] &= !e[0][i];
    }
    // NORMALISATION OF THE CURVE (X : Z) -> (X/Z : 1).
    fp_inv_own(&A.z, &A.z);
    //fp_inv(&A.z);
    fp_mul2(&A.x, &A.z);
    A.z = fp_1;



    do {

    	// POINT GENERATION
    	//assert(!memcmp(&A.z, &fp_1, sizeof(fp)));
        proj P;
        fp_random(&P.x);
        P.z = fp_1;

        // CHECKING POINT SIGN.
        fp rhs;
        //montgomery_rhs_proj(&rhs, &A, &P.x);
        montgomery_rhs_own(&rhs, &A.x, &P.x);
        //bool sign = !fp_issquare(&rhs);
        bool sign = isSquare(&rhs);


        if (done[sign])
            continue;

        // GENERATING THE POINT Q.
        xMUL_afin(&P, &A, &P, &k[sign]); // (Q =) P = kP.
        done[sign] = true;

        for (size_t i = 0; i < num_primes; ++i) {
            if (e[sign][i]) { // means we have to compute an l_i isogeny with sign "sign".

            	// COMPUTING THE MULTIPLICATIVE COEF.
                u512 cof = u512_1;
                for (size_t j = i + 1; j < num_primes; ++j){
                    if (e[sign][j])
                        u512_mul3_64(&cof, &cof, primes[j]);
                }

                // COMPUTING THE l_i TORSION POINT R.
                proj K;
                xMUL(&K, &A, &P, &cof);

                // COMPUTING THE ISOGENY AND THE IMAGE.
                if (memcmp(&K.z, &fp_0, sizeof(fp))) {
                    xISOG_2(&A, &P, &K, primes[i]);

                    // UPDATING THE COOEFFICIENT k.
                    if (!--e[sign][i])
                        u512_mul3_64(&k[sign], &k[sign], primes[i]);
                }
            }
            done[sign] &= !e[sign][i];
        }

        // NORMALISATION OF THE CURVE (X : Z) -> (X/Z : 1).
        fp_inv_own(&A.z, &A.z);
        //fp_inv(&A.z);
        fp_mul2(&A.x, &A.z);
        A.z = fp_1;

    } while (!(done[0] && done[1]));
    out->A = A.x;
}

/* totally not constant-time. */
void action_4(public_key *out, public_key const *in, private_key const *priv)
{
    u512 k[2];
    u512_set(&k[0], 4); /* maximal 2-power in p+1 */ // product of l_i with sign -
    u512_set(&k[1], 4); /* maximal 2-power in p+1 */ // product of l_i with sign +
    uint8_t e[2][num_primes];
    for (size_t i = 0; i < num_primes; ++i) {
        int8_t t = (int8_t) (priv->e[i / 2] << i % 2 * 4) >> 4;
        if (t > 0) {
            e[0][i] = t;
            e[1][i] = 0;
            u512_mul3_64(&k[1], &k[1], primes[i]);
        }
        else if (t < 0) {
            e[1][i] = -t;
            e[0][i] = 0;
            u512_mul3_64(&k[0], &k[0], primes[i]);
        }
        else {
            e[0][i] = 0;
            e[1][i] = 0;
            u512_mul3_64(&k[0], &k[0], primes[i]);
            u512_mul3_64(&k[1], &k[1], primes[i]);
        }
    }



    proj A = {in->A, fp_1};  // fp_1 = 1.
    bool done[2] = {false, false};

    // TORSION PART.
    
    fp Tx = {0x980BC417887630B8, 0x3223E5E8BC417B6A, 0x7239597EF4B7DEAC, 0x45AA9B1894E63200, 
    		0x45B6E360BBB5BCFD, 0x2303C882AD3A21D0, 0xF9B62E2CFCEBB476, 0x5C1DCF5FCF76F5E0};
	fp Tz = fp_1;
    proj T = {Tx, Tz};

    // GENERATING THE POINT Q.
    xMUL_afin(&T, &A, &T, &k[0]); // (Q =) P = kP.
    done[0] = true;

        for (size_t i = 0; i < num_primes; ++i) {
            if (e[0][i]) { // means we have to compute an l_i isogeny with sign "sign".

            	// COMPUTING THE MULTIPLICATIVE COEF.
                u512 cof = u512_1;
                for (size_t j = i + 1; j < num_primes; ++j){
                    if (e[0][j])
                        u512_mul3_64(&cof, &cof, primes[j]);
                }

            // COMPUTING THE l_i TORSION POINT R.
            proj T_i;
            xMUL(&T_i, &A, &T, &cof);

            // COMPUTING THE ISOGENY AND THE IMAGE.
        	xISOG_2(&A, &T, &T_i, primes[i]);

            // UPDATING THE COOEFFICIENT k.
            if (!--e[0][i])
                u512_mul3_64(&k[0], &k[0], primes[i]);
                
            }
            done[0] &= !e[0][i];
        }
    // NORMALISATION OF THE CURVE (X : Z) -> (X/Z : 1).
    fp_inv_own(&A.z, &A.z);
    //fp_inv(&A.z);
    fp_mul2(&A.x, &A.z);
    A.z = fp_1;



    do {

    	// POINT GENERATION
    	//assert(!memcmp(&A.z, &fp_1, sizeof(fp)));
        proj P;
        fp_random(&P.x);
        P.z = fp_1;

        // CHECKING POINT SIGN.
        fp rhs;
        //montgomery_rhs_proj(&rhs, &A, &P.x);
        montgomery_rhs_own(&rhs, &A.x, &P.x);
        //bool sign = !fp_issquare(&rhs);
        bool sign = isSquare(&rhs);


        if (done[sign])
            continue;

        // GENERATING THE POINT Q.
        xMUL_afin(&P, &A, &P, &k[sign]); // (Q =) P = kP.
        done[sign] = true;

        for (size_t i = 0; i < num_primes; ++i) {
            if (e[sign][i]) { // means we have to compute an l_i isogeny with sign "sign".

            	// COMPUTING THE MULTIPLICATIVE COEF.
                u512 cof = u512_1;
                for (size_t j = i + 1; j < num_primes; ++j){
                    if (e[sign][j])
                        u512_mul3_64(&cof, &cof, primes[j]);
                }

                // COMPUTING THE l_i TORSION POINT R.
                proj K;
                xMUL(&K, &A, &P, &cof);

                // COMPUTING THE ISOGENY AND THE IMAGE.
                if (memcmp(&K.z, &fp_0, sizeof(fp))) {
                    xISOG_2(&A, &P, &K, primes[i]);

                    // UPDATING THE COOEFFICIENT k.
                    if (!--e[sign][i])
                        u512_mul3_64(&k[sign], &k[sign], primes[i]);
                }
            }
            done[sign] &= !e[sign][i];
        }

        // NORMALISATION OF THE CURVE (X : Z) -> (X/Z : 1).
        fp_inv_own(&A.z, &A.z);
        //fp_inv(&A.z);
        fp_mul2(&A.x, &A.z);
        A.z = fp_1;

    } while (!(done[0] && done[1]));
    out->A = A.x;
}


/* totally not constant-time. */
void action_MR(public_key *out, public_key const *in, private_key const *priv)
{
    u512 k[2];
    u512_set(&k[0], 4); /* maximal 2-power in p+1 */
    u512_set(&k[1], 4); /* maximal 2-power in p+1 */

    uint8_t e[2][num_primes];

    for (size_t i = 0; i < num_primes; ++i) {

        int8_t t = (int8_t) (priv->e[i / 2] << i % 2 * 4) >> 4;

        if (t > 0) {
            e[0][i] = t;
            e[1][i] = 0;
            u512_mul3_64(&k[1], &k[1], primes[i]);
        }
        else if (t < 0) {
            e[1][i] = -t;
            e[0][i] = 0;
            u512_mul3_64(&k[0], &k[0], primes[i]);
        }
        else {
            e[0][i] = 0;
            e[1][i] = 0;
            u512_mul3_64(&k[0], &k[0], primes[i]);
            u512_mul3_64(&k[1], &k[1], primes[i]);
        }
    }

    proj A = {in->A, fp_1};

    bool done[2] = {false, false};

    do {

        assert(!memcmp(&A.z, &fp_1, sizeof(fp)));

        proj P;
        fp_random(&P.x);
        P.z = fp_1;

        fp rhs;
        montgomery_rhs(&rhs, &A.x, &P.x);
        bool sign = !fp_issquare(&rhs);

        if (done[sign])
            continue;

        xMUL(&P, &A, &P, &k[sign]);

        done[sign] = true;

        for (size_t i = num_primes-1; i < num_primes; --i) {  //changed loop direction

            if (e[sign][i]) {

                u512 cof = u512_1;
                for (size_t j = i - 1; j < num_primes; --j)   //changed loop direction
                    if (e[sign][j])
                        u512_mul3_64(&cof, &cof, primes[j]);

                proj K;
                xMUL(&K, &A, &P, &cof);

                if (memcmp(&K.z, &fp_0, sizeof(fp))) {

                    xISOG_MR(&A, &P, &K, primes[i]);

                    if (!--e[sign][i])
                        u512_mul3_64(&k[sign], &k[sign], primes[i]);

                }

            }

            done[sign] &= !e[sign][i];
        }

        fp_inv(&A.z);
        fp_mul2(&A.x, &A.z);
        A.z = fp_1;

    } while (!(done[0] && done[1]));

    out->A = A.x;
}

void action_MR_torsion(public_key *out, public_key const *in, private_key const *priv)
{
    u512 k[2];
    u512_set(&k[0], 4); /* maximal 2-power in p+1 */
    u512_set(&k[1], 4); /* maximal 2-power in p+1 */

    uint8_t e[2][num_primes];

    for (size_t i = 0; i < num_primes; ++i) {

        int8_t t = (int8_t) (priv->e[i / 2] << i % 2 * 4) >> 4;

        if (t > 0) {
            e[0][i] = t;
            e[1][i] = 0;
            u512_mul3_64(&k[1], &k[1], primes[i]);
        }
        else if (t < 0) {
            e[1][i] = -t;
            e[0][i] = 0;
            u512_mul3_64(&k[0], &k[0], primes[i]);
        }
        else {
            e[0][i] = 0;
            e[1][i] = 0;
            u512_mul3_64(&k[0], &k[0], primes[i]);
            u512_mul3_64(&k[1], &k[1], primes[i]);
        }
    }

    proj A = {in->A, fp_1};

    bool done[2] = {false, false};

    fp Tx = {0x980BC417887630B8, 0x3223E5E8BC417B6A, 0x7239597EF4B7DEAC, 0x45AA9B1894E63200, 
            0x45B6E360BBB5BCFD, 0x2303C882AD3A21D0, 0xF9B62E2CFCEBB476, 0x5C1DCF5FCF76F5E0};
    fp Tz = fp_1;
    proj T = {Tx, Tz};

    // GENERATING THE POINT Q.
    xMUL_afin(&T, &A, &T, &k[0]); // (Q =) P = kP.
    done[0] = true;

        for (size_t i = num_primes-1; i < num_primes; --i) {  //changed loop direction

            if (e[0][i]) {

                u512 cof = u512_1;
                for (size_t j = i - 1; j < num_primes; --j)   //changed loop direction
                    if (e[0][j])
                        u512_mul3_64(&cof, &cof, primes[j]);

                proj T_i;
                xMUL(&T_i, &A, &T, &cof);

                if (memcmp(&T_i.z, &fp_0, sizeof(fp))) {

                    xISOG_MR(&A, &T, &T_i, primes[i]);

                    if (!--e[0][i])
                        u512_mul3_64(&k[0], &k[0], primes[i]);

                }

            }

            done[0] &= !e[0][i];
        }
    // NORMALISATION OF THE CURVE (X : Z) -> (X/Z : 1).
    fp_inv_own(&A.z, &A.z);
    //fp_inv(&A.z);
    fp_mul2(&A.x, &A.z);
    A.z = fp_1;

    do {

        assert(!memcmp(&A.z, &fp_1, sizeof(fp)));

        proj P;
        fp_random(&P.x);
        P.z = fp_1;

        fp rhs;
        montgomery_rhs(&rhs, &A.x, &P.x);
        bool sign = !fp_issquare(&rhs);

        if (done[sign])
            continue;

        xMUL(&P, &A, &P, &k[sign]);

        done[sign] = true;

        for (size_t i = num_primes-1; i < num_primes; --i) {  //changed loop direction

            if (e[sign][i]) {

                u512 cof = u512_1;
                for (size_t j = i - 1; j < num_primes; --j)   //changed loop direction
                    if (e[sign][j])
                        u512_mul3_64(&cof, &cof, primes[j]);

                proj K;
                xMUL(&K, &A, &P, &cof);

                if (memcmp(&K.z, &fp_0, sizeof(fp))) {

                    xISOG_MR(&A, &P, &K, primes[i]);

                    if (!--e[sign][i])
                        u512_mul3_64(&k[sign], &k[sign], primes[i]);

                }

            }

            done[sign] &= !e[sign][i];
        }

        fp_inv(&A.z);
        fp_mul2(&A.x, &A.z);
        A.z = fp_1;

    } while (!(done[0] && done[1]));

    out->A = A.x;
}


/* includes public-key validation. */
bool csidh(public_key *out, public_key const *in, private_key const *priv)
{
    if (!validate(in)) {
        fp_random(&out->A);
        return false;
    }
    action(out, in, priv);
    return true;
}

/* includes public-key validation. */
bool csidh_2(public_key *out, public_key const *in, private_key const *priv)
{
    if (!validate_2(in)) {
        fp_random(&out->A);
        return false;
    }
    action_2(out, in, priv);
    return true;
}

bool csidh_3(public_key *out, public_key const *in, private_key const *priv)
{
    if (!validate_2(in)) {
        fp_random(&out->A);
        printf("FALSE !\n");
        return false;
    }
    action_3(out, in, priv);
    return true;
}

bool csidh_4(public_key *out, public_key const *in, private_key const *priv)
{
    if (!validate_2(in)) {
        fp_random(&out->A);
        printf("FALSE !\n");
        return false;
    }
    action_4(out, in, priv);
    return true;
}

bool csidh_MR(public_key *out, public_key const *in, private_key const *priv)
{
    if (!validate_2(in)) {
        fp_random(&out->A);
        printf("FALSE !\n");
        return false;
    }
    action_MR(out, in, priv);
    return true;
}
bool csidh_MR_torsion(public_key *out, public_key const *in, private_key const *priv)
{
    if (!validate_2(in)) {
        fp_random(&out->A);
        printf("FALSE !\n");
        return false;
    }
    action_MR_torsion(out, in, priv);
    return true;
}