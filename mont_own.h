#ifndef MONTOWN_H
#define MONTOWN_H

#include "u512.h"
#include "fp.h"
#include "mont.h"

/* P^1 over fp. */
// typedef struct proj {
//     fp x;
//     fp z;
// } proj;

void xDBL_own(proj *Q, proj const *A, proj const *P);
void xDBL24(proj *Q, proj const *A, proj const *P);
// void xADD_own(proj *S, proj const *P, proj const *Q, proj const *PQ);
void xADD_zafin(proj *S, proj const *P, proj const *Q, proj const *PQ);
// void xMUL(proj *Q, proj const *A, proj const *P, u512 const *k);
void xDBLADD_zafin(proj *R, proj *S, proj const *P, proj const *Q, proj const *PQ, proj const *A);
void xDBLADD_afin(proj *R, proj *S, proj const *P, proj const *Q, proj const *PQ, proj const *A);
void xDBLADD_A24afin(proj *R, proj *S, proj const *P, proj const *Q, proj const *PQ, proj const *A);

void xMUL_zafin(proj *Q, proj const *A, proj const *P, u512 const *k);
void xMUL_afin(proj *Q, proj const *A, proj const *P, u512 const *k);
void xMUL_A24afin(proj *Q, proj const *A, proj const *P, u512 const *k);


// void xISOG(proj *A, proj *P, proj const *K, uint64_t k);
void xISOG_2(proj *A, proj *P, proj const *K, uint64_t k);
#endif
