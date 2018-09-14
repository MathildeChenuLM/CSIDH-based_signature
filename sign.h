//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
//            Signatures based on CSIDH (with Unruh transform)
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
// Written by Mathilde, started on 24/08/18.

#ifndef SIGN_H
#define SIGN_H

#include "u512.h"
#include "fp.h"
#include "mont.h"
#include "mont_own.h"
#include "csidh.h"

#define nb_of_comm 128
#define half_nb_of_primes 37

typedef struct signature {
    unsigned char *hash;
    int8_t resp_sent[nb_of_comm][half_nb_of_primes];
    unsigned char Gresp_sent[nb_of_comm][half_nb_of_primes];
} signature;

void signature_init(signature *sig);
    /* Initialise a signature
    */

void signature_clear(signature *sig);
    /* Free the space allocated to a signature.
    */

void hash_G(int8_t in_int8[half_nb_of_primes], unsigned char *out);
 	/* Transforms a private key into a string, and returns its hash.
 	The G hash function is designed specificaly for answers in Unruh transform,
 	since the size of the input is the size of the output. 
 	*/

void hash_H(unsigned char *out, uint64_t comm[nb_of_comm][8], 
	unsigned char Gresp0[nb_of_comm][half_nb_of_primes], 
	unsigned char Gresp1[nb_of_comm][half_nb_of_primes]);
 	/* Returns the hash of the concatenation of all the elements listed above.
 	H(comm, Gresp0, Gresp1). 
 	*/

void sign_original( signature *sig, private_key *priv, public_key *base );
    /* Signs an empty message according to the private key of the signer.
    */

void verify_original( signature *sig, public_key *pub_a, public_key *base);
    /* Verifies that the signature is valid.
    Should ne preceeded by public key validation.
    */

void sign_xwing( signature *sig, private_key *priv, public_key *base );
    /* Signs an empty message according to the private key of the signer.
    */

void verify_xwing( signature *sig, public_key *pub_a, public_key *base);
    /* Verifies that the signature is valid.
    Should ne preceeded by public key validation.
    */

void sign_torsion( signature *sig, private_key *priv, public_key *base );
    /* Signs an empty message according to the private key of the signer.
    */

void verify_torsion( signature *sig, public_key *pub_a, public_key *base);
    /* Verifies that the signature is valid.
    Should ne preceeded by public key validation.
    */

void sign_MR( signature *sig, private_key *priv, public_key *base );
    /* Signs an empty message according to the private key of the signer.
    */

void verify_MR( signature *sig, public_key *pub_a, public_key *base);
    /* Verifies that the signature is valid.
    Should ne preceeded by public key validation.
    */

void sign_MR_torsion( signature *sig, private_key *priv, public_key *base );
    /* Signs an empty message according to the private key of the signer.
    */

void verify_MR_torsion( signature *sig, public_key *pub_a, public_key *base);
    /* Verifies that the signature is valid.
    Should ne preceeded by public key validation.
    */

#endif
