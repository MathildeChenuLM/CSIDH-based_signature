#ifndef MONT_H
#define MONT_H

#include "u512.h"
#include "fp.h"

/* P^1 over fp. */
typedef struct proj {
    fp x;
    fp z;
} proj;

void xDBL(proj *Q, proj const *A, proj const *P);
//void xDBL24(proj *Q, proj const *A, proj const *P);
void xADD(proj *S, proj const *P, proj const *Q, proj const *PQ);
//void xADD_zafin(proj *S, proj const *P, proj const *Q, proj const *PQ);
void xDBLADD(proj *R, proj *S, proj const *P, proj const *Q, proj const *PQ, proj const *A);
void xMUL(proj *Q, proj const *A, proj const *P, u512 const *k);
//void xDBLADD_zafin(proj *R, proj *S, proj const *P, proj const *Q, proj const *PQ, proj const *A);
// void xMUL_zafin(proj *Q, proj const *A, proj const *P, u512 const *k);
void xISOG(proj *A, proj *P, proj const *K, uint64_t k); // P = phi(P) where phi is the isogeny of kernel <K>
void exp_by_squaring_(fp* x, fp* y, uint64_t exp);
void xISOG_MR(proj *A, proj *P, proj const *K, uint64_t k); // P = phi(P) where phi is the isogeny of kernel <K>
// void xISOG_2(proj *A, proj *P, proj const *K, uint64_t k); // same with X-wing :)
#endif
