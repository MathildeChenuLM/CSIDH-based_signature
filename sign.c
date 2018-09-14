//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
//            Signatures based on CSIDH (with Unruh transform)
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
// Written by Mathilde, started on 24/08/18. 

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "compact_shake.c"
#include "csidh.h"
#include "tools_sign.h"
#include "sign.h"
#include "rng.h"

void signature_init(signature *sig){
    /* Initialise a signature.
    */
    sig->hash = (unsigned char *) malloc(32);
}

void signature_clear(signature *sig){
    /* Free the space allocated to a signature.
    */
    free(sig->hash);
}

void hash_G(int8_t in_int8[half_nb_of_primes], unsigned char *out){ 
 	/* Transforms a private key into a string, and returns its hash.
 	The G hash function is designed specificaly for answers in Unruh transform,
 	since the size of the input is the size of the output. 
 	*/
	unsigned char in_char[half_nb_of_primes];

	for(int i = 0; i<half_nb_of_primes; i++){
		in_char[i] = (unsigned char) in_int8[i];
	}
	// Keccak with rate 1088, capacity 512, inLen and outLen 37.
	Keccak(1088, in_char, half_nb_of_primes, 0x1F, out, half_nb_of_primes); 

    return;
}

void hash_H(unsigned char *out, uint64_t comm[nb_of_comm][8], 
	unsigned char Gresp0[nb_of_comm][half_nb_of_primes], 
	unsigned char Gresp1[nb_of_comm][half_nb_of_primes]){ 
 	/* Returns the hash of the concatenation of all the elements listed above.
 	H(comm, Gresp0, Gresp1). 
 	*/

 	int size = (64 +2*half_nb_of_primes)*nb_of_comm;
	unsigned char buf[size];

	//comm_i
	unsigned char tmp[64];
	for(int i = 0; i<nb_of_comm; i++){
		from_tab_int64_to_tab_char(tmp, comm[i]);
		concat_char(size, buf, 64, tmp, i*64);
	}

	//G(resp0)	
	for(int i = 0; i<nb_of_comm; i++){
		concat_char(size, buf, half_nb_of_primes, Gresp0[i], (64*nb_of_comm + half_nb_of_primes*i) ) ;
	}

	//G(resp1)
	for(int i = 0; i<nb_of_comm; i++){
		concat_char(size, buf, half_nb_of_primes, Gresp1[i], (64 + half_nb_of_primes)*nb_of_comm + half_nb_of_primes*i);
	}

	// Keccak with rate 1088, capacity 512, inLen and outLen 32.
	Keccak(1088, buf, size, 0x1F, out, 32); 

    return;
} 

void sign_original( signature *sig, private_key *priv, public_key *base ){
    /* Signs an empty message according to the private key of the signer.
    */

    uint64_t comm[nb_of_comm][8]; // public_keys

    int8_t resp0[nb_of_comm][half_nb_of_primes]; // private keys 
    int8_t resp1[nb_of_comm][half_nb_of_primes]; // private keys

    unsigned char Gresp0[nb_of_comm][half_nb_of_primes];
    unsigned char Gresp1[nb_of_comm][half_nb_of_primes];

    private_key resp0_i;
    public_key res_action;
    for(int i=0; i<nb_of_comm; i++){

        ideal_gen(&resp0_i);
        copy_tab_int8(half_nb_of_primes, resp0[i], resp0_i.e);
        substract_tab(half_nb_of_primes, resp1[i], resp0[i], priv->e );
        hash_G(resp0[i], Gresp0[i]);
        hash_G(resp1[i], Gresp1[i]);

        action(&res_action, base, &resp0_i);

        for(int j=0; j<8; j++){
            comm[i][j] = (((res_action.A).x).c)[j];
        }
    }   

    hash_H(sig->hash, comm, Gresp0, Gresp1);
    
    for(int i=0; i<nb_of_comm; i++){

        if( bit_at_i(sig->hash, i) == 0 ){
            // the challenge is 0;
            copy_tab_int8(half_nb_of_primes, (sig->resp_sent)[i], resp0[i]);
            copy_tab_uchar(half_nb_of_primes, (sig->Gresp_sent)[i], Gresp1[i]); 
        }
        else{
            // the challenge is 1;
            copy_tab_int8(half_nb_of_primes, (sig->resp_sent)[i], resp1[i]);
            copy_tab_uchar(half_nb_of_primes, (sig->Gresp_sent)[i], Gresp0[i]); 
        }
    }

    return;
}



void sign_xwing( signature *sig, private_key *priv, public_key *base ){
    /* Signs an empty message according to the private key of the signer.
    */

    uint64_t comm[nb_of_comm][8]; // public_keys

    int8_t resp0[nb_of_comm][half_nb_of_primes]; // private keys 
    int8_t resp1[nb_of_comm][half_nb_of_primes]; // private keys

    unsigned char Gresp0[nb_of_comm][half_nb_of_primes];
    unsigned char Gresp1[nb_of_comm][half_nb_of_primes];

    private_key resp0_i;
    public_key res_action;
    for(int i=0; i<nb_of_comm; i++){

    	ideal_gen(&resp0_i);
    	copy_tab_int8(half_nb_of_primes, resp0[i], resp0_i.e);
    	substract_tab(half_nb_of_primes, resp1[i], resp0[i], priv->e );
    	hash_G(resp0[i], Gresp0[i]);
    	hash_G(resp1[i], Gresp1[i]);

    	action_2(&res_action, base, &resp0_i);

    	for(int j=0; j<8; j++){
    		comm[i][j] = (((res_action.A).x).c)[j];
    	}
    }   

    hash_H(sig->hash, comm, Gresp0, Gresp1);
    
    for(int i=0; i<nb_of_comm; i++){

    	if( bit_at_i(sig->hash, i) == 0 ){
    		// the challenge is 0;
    		copy_tab_int8(half_nb_of_primes, (sig->resp_sent)[i], resp0[i]);
    		copy_tab_uchar(half_nb_of_primes, (sig->Gresp_sent)[i], Gresp1[i]); 
    	}
    	else{
    		// the challenge is 1;
    		copy_tab_int8(half_nb_of_primes, (sig->resp_sent)[i], resp1[i]);
    		copy_tab_uchar(half_nb_of_primes, (sig->Gresp_sent)[i], Gresp0[i]); 
    	}
    }

    return;
}

void sign_torsion( signature *sig, private_key *priv, public_key *base ){
    /* Signs an empty message according to the private key of the signer.
    */

    uint64_t comm[nb_of_comm][8]; // public_keys

    int8_t resp0[nb_of_comm][half_nb_of_primes]; // private keys 
    int8_t resp1[nb_of_comm][half_nb_of_primes]; // private keys

    unsigned char Gresp0[nb_of_comm][half_nb_of_primes];
    unsigned char Gresp1[nb_of_comm][half_nb_of_primes];

    private_key resp0_i;
    public_key res_action;
    for(int i=0; i<nb_of_comm; i++){

        ideal_gen(&resp0_i);
        copy_tab_int8(half_nb_of_primes, resp0[i], resp0_i.e);
        substract_tab(half_nb_of_primes, resp1[i], resp0[i], priv->e );
        hash_G(resp0[i], Gresp0[i]);
        hash_G(resp1[i], Gresp1[i]);

        action_4(&res_action, base, &resp0_i);

        for(int j=0; j<8; j++){
            comm[i][j] = (((res_action.A).x).c)[j];
        }
    }   

    hash_H(sig->hash, comm, Gresp0, Gresp1);
    
    for(int i=0; i<nb_of_comm; i++){

        if( bit_at_i(sig->hash, i) == 0 ){
            // the challenge is 0;
            copy_tab_int8(half_nb_of_primes, (sig->resp_sent)[i], resp0[i]);
            copy_tab_uchar(half_nb_of_primes, (sig->Gresp_sent)[i], Gresp1[i]); 
        }
        else{
            // the challenge is 1;
            copy_tab_int8(half_nb_of_primes, (sig->resp_sent)[i], resp1[i]);
            copy_tab_uchar(half_nb_of_primes, (sig->Gresp_sent)[i], Gresp0[i]); 
        }
    }

    return;
}

void sign_MR( signature *sig, private_key *priv, public_key *base ){
    /* Signs an empty message according to the private key of the signer.
    */

    uint64_t comm[nb_of_comm][8]; // public_keys

    int8_t resp0[nb_of_comm][half_nb_of_primes]; // private keys 
    int8_t resp1[nb_of_comm][half_nb_of_primes]; // private keys

    unsigned char Gresp0[nb_of_comm][half_nb_of_primes];
    unsigned char Gresp1[nb_of_comm][half_nb_of_primes];

    private_key resp0_i;
    public_key res_action;
    for(int i=0; i<nb_of_comm; i++){

        ideal_gen(&resp0_i);
        copy_tab_int8(half_nb_of_primes, resp0[i], resp0_i.e);
        substract_tab(half_nb_of_primes, resp1[i], resp0[i], priv->e );
        hash_G(resp0[i], Gresp0[i]);
        hash_G(resp1[i], Gresp1[i]);

        action_MR(&res_action, base, &resp0_i);

        for(int j=0; j<8; j++){
            comm[i][j] = (((res_action.A).x).c)[j];
        }
    }   

    hash_H(sig->hash, comm, Gresp0, Gresp1);
    
    for(int i=0; i<nb_of_comm; i++){

        if( bit_at_i(sig->hash, i) == 0 ){
            // the challenge is 0;
            copy_tab_int8(half_nb_of_primes, (sig->resp_sent)[i], resp0[i]);
            copy_tab_uchar(half_nb_of_primes, (sig->Gresp_sent)[i], Gresp1[i]); 
        }
        else{
            // the challenge is 1;
            copy_tab_int8(half_nb_of_primes, (sig->resp_sent)[i], resp1[i]);
            copy_tab_uchar(half_nb_of_primes, (sig->Gresp_sent)[i], Gresp0[i]); 
        }
    }

    return;
}

void sign_MR_torsion( signature *sig, private_key *priv, public_key *base ){
    /* Signs an empty message according to the private key of the signer.
    */

    uint64_t comm[nb_of_comm][8]; // public_keys

    int8_t resp0[nb_of_comm][half_nb_of_primes]; // private keys 
    int8_t resp1[nb_of_comm][half_nb_of_primes]; // private keys

    unsigned char Gresp0[nb_of_comm][half_nb_of_primes];
    unsigned char Gresp1[nb_of_comm][half_nb_of_primes];

    private_key resp0_i;
    public_key res_action;
    for(int i=0; i<nb_of_comm; i++){

        ideal_gen(&resp0_i);
        copy_tab_int8(half_nb_of_primes, resp0[i], resp0_i.e);
        substract_tab(half_nb_of_primes, resp1[i], resp0[i], priv->e );
        hash_G(resp0[i], Gresp0[i]);
        hash_G(resp1[i], Gresp1[i]);

        action_MR_torsion(&res_action, base, &resp0_i);

        for(int j=0; j<8; j++){
            comm[i][j] = (((res_action.A).x).c)[j];
        }
    }   

    hash_H(sig->hash, comm, Gresp0, Gresp1);
    
    for(int i=0; i<nb_of_comm; i++){

        if( bit_at_i(sig->hash, i) == 0 ){
            // the challenge is 0;
            copy_tab_int8(half_nb_of_primes, (sig->resp_sent)[i], resp0[i]);
            copy_tab_uchar(half_nb_of_primes, (sig->Gresp_sent)[i], Gresp1[i]); 
        }
        else{
            // the challenge is 1;
            copy_tab_int8(half_nb_of_primes, (sig->resp_sent)[i], resp1[i]);
            copy_tab_uchar(half_nb_of_primes, (sig->Gresp_sent)[i], Gresp0[i]); 
        }
    }

    return;
}

void verify_original( signature *sig, public_key *pub_a, public_key *base){
    /* Verifies that the signature is valid.
    Should ne preceeded by public key validation.
    */

    uint64_t comm[nb_of_comm][8]; // public_keys

    unsigned char Gresp0[nb_of_comm][half_nb_of_primes];
    unsigned char Gresp1[nb_of_comm][half_nb_of_primes];

    unsigned char *digest;
    digest = (unsigned char *) malloc(32);

    private_key resp_i;
    public_key comm_i;

    for(int i=0; i<nb_of_comm; i++){

        //transforming sig.resp_sent[i] into a private key. 
        for(int j=0; j<half_nb_of_primes; j++){
             (resp_i.e)[j] = ((sig->resp_sent)[i])[j];
        }

        if( bit_at_i(sig->hash, i) == 0 ){
            action(&comm_i, base, &resp_i);
        }
        else{
            action(&comm_i, pub_a, &resp_i);    
        }

        //transforming the public key into int8_t[nb_of_comm][half_nb_of_primes]
        for(int j=0; j<8; j++){
            comm[i][j] = (((comm_i.A).x).c)[j];
        }
    }

    // GRESP PART (OK)
    for(int i=0; i<nb_of_comm; i++){
        if( bit_at_i(sig->hash, i) == 0 ){
            hash_G( (sig->resp_sent)[i], Gresp0[i] );
            copy_tab_uchar(half_nb_of_primes, Gresp1[i], (sig->Gresp_sent)[i]);
        }
        else{
            copy_tab_uchar(half_nb_of_primes, Gresp0[i], (sig->Gresp_sent)[i]);
            hash_G( (sig->resp_sent)[i], Gresp1[i] );       
        }
    }

    hash_H(digest, comm, Gresp0, Gresp1);

    printf("Test ! Compare the two hash : \n");
    print_hash(digest, 32); // 32 char = 32*8 = 256 bits.
    printf("\n");
    print_hash(sig->hash, 32);
    printf("\n");

    free(digest);
    return;
}

void verify_xwing( signature *sig, public_key *pub_a, public_key *base){
    /* Verifies that the signature is valid.
    Should ne preceeded by public key validation.
    */

	uint64_t comm[nb_of_comm][8]; // public_keys

    unsigned char Gresp0[nb_of_comm][half_nb_of_primes];
    unsigned char Gresp1[nb_of_comm][half_nb_of_primes];

    unsigned char *digest;
    digest = (unsigned char *) malloc(32);

    private_key resp_i;
    public_key comm_i;

    for(int i=0; i<nb_of_comm; i++){

        //transforming sig.resp_sent[i] into a private key. 
        for(int j=0; j<half_nb_of_primes; j++){
             (resp_i.e)[j] = ((sig->resp_sent)[i])[j];
        }

        if( bit_at_i(sig->hash, i) == 0 ){
            action_2(&comm_i, base, &resp_i);
        }
        else{
            action_2(&comm_i, pub_a, &resp_i);    
        }

        //transforming the public key into int8_t[nb_of_comm][half_nb_of_primes]
        for(int j=0; j<8; j++){
            comm[i][j] = (((comm_i.A).x).c)[j];
        }
    }

    // GRESP PART (OK)
    for(int i=0; i<nb_of_comm; i++){
        if( bit_at_i(sig->hash, i) == 0 ){
            hash_G( (sig->resp_sent)[i], Gresp0[i] );
            copy_tab_uchar(half_nb_of_primes, Gresp1[i], (sig->Gresp_sent)[i]);
        }
        else{
            copy_tab_uchar(half_nb_of_primes, Gresp0[i], (sig->Gresp_sent)[i]);
            hash_G( (sig->resp_sent)[i], Gresp1[i] );       
        }
    }

    hash_H(digest, comm, Gresp0, Gresp1);

    printf("Test ! Compare the two hash : \n");
    print_hash(digest, 32); // 32 char = 32*8 = 256 bits.
    printf("\n");
    print_hash(sig->hash, 32);
    printf("\n");

    free(digest);
    return;
}

void verify_torsion( signature *sig, public_key *pub_a, public_key *base){
    /* Verifies that the signature is valid.
    Should ne preceeded by public key validation.
    */

    uint64_t comm[nb_of_comm][8]; // public_keys

    unsigned char Gresp0[nb_of_comm][half_nb_of_primes];
    unsigned char Gresp1[nb_of_comm][half_nb_of_primes];

    unsigned char *digest;
    digest = (unsigned char *) malloc(32);

    private_key resp_i;
    public_key comm_i;

    for(int i=0; i<nb_of_comm; i++){

        //transforming sig.resp_sent[i] into a private key. 
        for(int j=0; j<half_nb_of_primes; j++){
             (resp_i.e)[j] = ((sig->resp_sent)[i])[j];
        }

        if( bit_at_i(sig->hash, i) == 0 ){
            action_4(&comm_i, base, &resp_i);
        }
        else{
            action_2(&comm_i, pub_a, &resp_i);    
        }

        //transforming the public key into int8_t[nb_of_comm][half_nb_of_primes]
        for(int j=0; j<8; j++){
            comm[i][j] = (((comm_i.A).x).c)[j];
        }
    }

    // GRESP PART (OK)
    for(int i=0; i<nb_of_comm; i++){
        if( bit_at_i(sig->hash, i) == 0 ){
            hash_G( (sig->resp_sent)[i], Gresp0[i] );
            copy_tab_uchar(half_nb_of_primes, Gresp1[i], (sig->Gresp_sent)[i]);
        }
        else{
            copy_tab_uchar(half_nb_of_primes, Gresp0[i], (sig->Gresp_sent)[i]);
            hash_G( (sig->resp_sent)[i], Gresp1[i] );       
        }
    }

    hash_H(digest, comm, Gresp0, Gresp1);

    printf("Test ! Compare the two hash : \n");
    print_hash(digest, 32); // 32 char = 32*8 = 256 bits.
    printf("\n");
    print_hash(sig->hash, 32);
    printf("\n");

    free(digest);
    return;
}

void verify_MR( signature *sig, public_key *pub_a, public_key *base){
    /* Verifies that the signature is valid.
    Should ne preceeded by public key validation.
    */

    uint64_t comm[nb_of_comm][8]; // public_keys

    unsigned char Gresp0[nb_of_comm][half_nb_of_primes];
    unsigned char Gresp1[nb_of_comm][half_nb_of_primes];

    unsigned char *digest;
    digest = (unsigned char *) malloc(32);

    private_key resp_i;
    public_key comm_i;

    for(int i=0; i<nb_of_comm; i++){

        //transforming sig.resp_sent[i] into a private key. 
        for(int j=0; j<half_nb_of_primes; j++){
             (resp_i.e)[j] = ((sig->resp_sent)[i])[j];
        }

        if( bit_at_i(sig->hash, i) == 0 ){
            action_MR(&comm_i, base, &resp_i);
        }
        else{
            action_MR(&comm_i, pub_a, &resp_i);    
        }

        //transforming the public key into int8_t[nb_of_comm][half_nb_of_primes]
        for(int j=0; j<8; j++){
            comm[i][j] = (((comm_i.A).x).c)[j];
        }
    }

    // GRESP PART (OK)
    for(int i=0; i<nb_of_comm; i++){
        if( bit_at_i(sig->hash, i) == 0 ){
            hash_G( (sig->resp_sent)[i], Gresp0[i] );
            copy_tab_uchar(half_nb_of_primes, Gresp1[i], (sig->Gresp_sent)[i]);
        }
        else{
            copy_tab_uchar(half_nb_of_primes, Gresp0[i], (sig->Gresp_sent)[i]);
            hash_G( (sig->resp_sent)[i], Gresp1[i] );       
        }
    }

    hash_H(digest, comm, Gresp0, Gresp1);

    printf("Test ! Compare the two hash : \n");
    print_hash(digest, 32); // 32 char = 32*8 = 256 bits.
    printf("\n");
    print_hash(sig->hash, 32);
    printf("\n");

    free(digest);
    return;
}

void verify_MR_torsion( signature *sig, public_key *pub_a, public_key *base){
    /* Verifies that the signature is valid.
    Should ne preceeded by public key validation.
    */

    uint64_t comm[nb_of_comm][8]; // public_keys

    unsigned char Gresp0[nb_of_comm][half_nb_of_primes];
    unsigned char Gresp1[nb_of_comm][half_nb_of_primes];

    unsigned char *digest;
    digest = (unsigned char *) malloc(32);

    private_key resp_i;
    public_key comm_i;

    for(int i=0; i<nb_of_comm; i++){

        //transforming sig.resp_sent[i] into a private key. 
        for(int j=0; j<half_nb_of_primes; j++){
             (resp_i.e)[j] = ((sig->resp_sent)[i])[j];
        }

        if( bit_at_i(sig->hash, i) == 0 ){
            action_MR_torsion(&comm_i, base, &resp_i);
        }
        else{
            action_MR(&comm_i, pub_a, &resp_i);    
        }

        //transforming the public key into int8_t[nb_of_comm][half_nb_of_primes]
        for(int j=0; j<8; j++){
            comm[i][j] = (((comm_i.A).x).c)[j];
        }
    }

    // GRESP PART (OK)
    for(int i=0; i<nb_of_comm; i++){
        if( bit_at_i(sig->hash, i) == 0 ){
            hash_G( (sig->resp_sent)[i], Gresp0[i] );
            copy_tab_uchar(half_nb_of_primes, Gresp1[i], (sig->Gresp_sent)[i]);
        }
        else{
            copy_tab_uchar(half_nb_of_primes, Gresp0[i], (sig->Gresp_sent)[i]);
            hash_G( (sig->resp_sent)[i], Gresp1[i] );       
        }
    }

    hash_H(digest, comm, Gresp0, Gresp1);

    printf("Test ! Compare the two hash : \n");
    print_hash(digest, 32); // 32 char = 32*8 = 256 bits.
    printf("\n");
    print_hash(sig->hash, 32);
    printf("\n");

    free(digest);
    return;
}



