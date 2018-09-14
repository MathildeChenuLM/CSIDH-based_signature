//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
//        Tools for signatures based on CSIDH (with Unruh transform)
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
// Written by Mathilde, started on 28/08/18.

#include <stdint.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "csidh.h"
#include "tools_sign.h"
#include "sign.h"
#include "rng.h"

bool bit_at_i(unsigned char *digest, int i){
	/* Returns the ith bit of an unsigned char.
	Used on a digest to find the value of the committment. 
	*/
	return ( ( digest[i/8] >> (i%8) ) & 1 );
}

void copy_tab_int8(int size, int8_t dest[size], int8_t src[size]){
    /* Copies an array of int8_t into another array of int8_t.
    Used to copy resp_i.
    */

	for(int i=0; i<size; i++){
		dest[i] = src[i];
	}
}

void copy_tab_uchar(int size, unsigned char dest[size], unsigned char src[size]){
    /* Copies an array of unsigned char into another array of unsigned char.
    Used to copy Gresp_i.
    */
	for(int i=0; i<size; i++){
		dest[i] = src[i];
	}
}

void copy_tab_int64(int size, int64_t dest[size], uint64_t src[size]){
    /* Copies an array of int64 into another array of int_64.
    Used to copy comm_i.
    */
    for(int i=0; i<size; i++){
        dest[i] = src[i];
    }
}


void from_tab_int64_to_tab_char(unsigned char p_char[64], uint64_t p[8]){
	/* Converts an array int64_t[8] into a char[64].
	*/
	for(int i = 0; i < 8; i++){
		for(int j = 0; j < 8; j++){
			p_char[j + i*8] = (char) ( p[i] >> ((7-j)*8) );
		}
	}
}

void concat_char(int out_size, unsigned char out[out_size], 
	int in_size, unsigned char in[in_size], int start){
    /* Concatenates in in the end of out, 
    where in and out are two arrays of chars of respective size out_size and in_size.
    */
	for(int i = 0; i<in_size; i++){
		out[i+start] = in[i];
	}
}

void print_char_tab(int size, unsigned char c[size]){
	/* Prints a string, ie a char[size] with representation in hexa.
	*/
	for(int i = 0; i < size; i++){
		printf("%02x", c[i]);
	}
	printf("\n");
}

void print_int8_tab(int size, int8_t tab[size]){
	/* Prints an array of int8_t, ie a int8_t[size] with representation in hexa.
	*/
	for(int i = 0; i < size; i++){
		printf("%02x", tab[i]);
	}
	printf("\n");
}

void print_uint64_tab(int size, uint64_t tab[size]){
	/* Prints a string, ie a char[size] with representation in hexa.
	*/
	for(int i = 0; i < size; i++){
		printf("%02lx", tab[i]);
	}
	printf("\n");
}

void print_hash(unsigned char *out, int outLen){
	/* Prints the content of a digest.
	(Redondant with print_char, but more explicit name).
    outLen should be nb_of_comm / 8;
	*/

	for(int i = 0; i < outLen; i++){
		printf("%02x ", out[i]); //02 so that 6 would appear as 06.
	}
	printf("\n");
}

void ideal_gen(private_key *priv){
    /* Random ideal generation. 
    Used to generate a random committment, which exponents are in |[-2 ; 2]|.
    Function written by csidh code author, modified by Mathilde.
    WATCH OUT : Keep the exponents in this bound, 
        otherwise the answer of challenge 1 woudn't fit in int4 ! (5+2 = 7). 
    */
    memset(&priv->e, 0, sizeof(priv->e));
    for (size_t i = 0; i < num_primes; ) {
        int8_t buf[64];
        randombytes(buf, sizeof(buf));
        for (size_t j = 0; j < sizeof(buf); ++j) {
            if (buf[j] <= 2 && buf[j] >= -2) {
                priv->e[i / 2] |= (buf[j] & 0xf) << i % 2 * 4;
                if (++i >= num_primes)
                    break;
            }
        }
    }
}

int8_t substract(int8_t e1, int8_t e2){
	/* Returns the substraction of the exponents packed in e1 by those in e2.
	Used to compute i1 // i2, where i1 and i2 are two ideals whose exponents are e1 and e2.
	WATCH OUT ! If the exponents of i1 are between -5 and 5, those of i2 should be 
		between -2 and 2, in order to avoid overflow.
	*/
	int8_t e3 = 0; // initialisation.

	for(int i = 0; i < 2; i++){

		int8_t t1 = (int8_t) ( e1 << ( (i%2)*4 ) ) >> 4; // unpacks the ith exponent.
		int8_t t2 = (int8_t) ( e2 << ( (i%2)*4 ) ) >> 4; // idem.
		int8_t t3 = ((t1<<4) - (t2<<4)) >> 4; // Necessary to shift to keep the right sign.

		// packing again.
		if(! i%2 ){
			e3 = t3 << 4;
		}
		else{
			e3 = e3 ^ (t3 & 0xf);
		}

	}

	return(e3);
}

void substract_tab(int N, int8_t e3[N], int8_t e1[N], int8_t e2[N]){
	/* Returns the substraction of the exponents packed in e1 by those in e2.
	Exponents are really private key, ie ideal factorised representation.
	Used to compute i1 // i2, where i1 and i2 are two ideals whose exponents are e1 and e2.
	WATCH OUT ! If the exponents of i1 are between -5 and 5, those of i2 should be 
		between -2 and 2, in order to avoid overflow.
	*/
	for(int i = 0; i < N; i++){
		e3[i] = substract(e1[i], e2[i]);
	}
}