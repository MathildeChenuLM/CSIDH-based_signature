#ifndef CSIDH_H
#define CSIDH_H

#include "u512.h"
#include "fp.h"
#include "mont.h"
#include "mont_own.h"

/* specific to p, should perhaps be somewhere else */
#define num_primes 74
#define max_exponent 5 /* (2*5+1)^74 is roughly 2^256 */


typedef struct private_key {
    int8_t e[(num_primes + 1) / 2]; /* packed int4_t */
} private_key;

typedef struct public_key {
    fp A; /* Montgomery coefficient: represents y^2 = x^3 + Ax^2 + x */
} public_key;

extern const public_key base;

void csidh_private(private_key *priv);
bool csidh(public_key *out, public_key const *in, private_key const *priv);
bool csidh_2(public_key *out, public_key const *in, private_key const *priv);
bool csidh_3(public_key *out, public_key const *in, private_key const *priv);
bool csidh_4(public_key *out, public_key const *in, private_key const *priv);
bool csidh_MR(public_key *out, public_key const *in, private_key const *priv);
bool csidh_MR_torsion(public_key *out, public_key const *in, private_key const *priv);

bool validate(public_key const *in);
void action(public_key *out, public_key const *in, private_key const *priv);
void action_2(public_key *out, public_key const *in, private_key const *priv);
void action_3(public_key *out, public_key const *in, private_key const *priv);
void action_4(public_key *out, public_key const *in, private_key const *priv);
void action_MR(public_key *out, public_key const *in, private_key const *priv);
void action_MR_torsion(public_key *out, public_key const *in, private_key const *priv);


#endif
