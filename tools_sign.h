//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
//        Tools for signatures based on CSIDH (with Unruh transform)
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
// Written by Mathilde, started on 28/08/18.

#ifndef TOOLS_SIGN_H
#define TOOLS_SIGN_H

#include "u512.h"
#include "fp.h"
#include "mont.h"
#include "mont_own.h"
#include "csidh.h"

bool bit_at_i(unsigned char *digest, int i);
	/* Returns the ith bit of an unsigned char.
	Used on a digest to find the value of the committment. 
	*/

void copy_tab_int8(int size, int8_t dest[size], int8_t src[size]);
	/* Copies an array of int8_t into another array of int8_t.
	Used to copy resp_i.
    */

void copy_tab_uchar(int size, unsigned char dest[size], unsigned char src[size]);
    /* Copies an array of unsigned char into another array of unsigned char.
    Used to copy Gresp_i.
    */

void copy_tab_int64(int size, int64_t dest[size], uint64_t src[size]);
    /* Copies an array of int64 into another array of int_64.
    Used to copy comm_i.
    */


void from_tab_int64_to_tab_char(unsigned char p_char[64], uint64_t p[8]);
	/* Converts an array int64_t[8] into a char[64].
	*/

void concat_char(int out_size, unsigned char out[out_size], 
	int in_size, unsigned char in[in_size], int start);
    /* Concatenates in in the end of out, 
    where in and out are two arrays of chars of respective size out_size and in_size.
    */

void print_char_tab(int size, unsigned char c[size]);
	/* Prints a string, ie a char[size] with representation in hexa.
	*/

void print_int8_tab(int size, int8_t tab[size]);
	/* Prints an array of int8_t, ie a int8_t[size] with representation in hexa.
	*/

void print_uint64_tab(int size, uint64_t tab[size]);
	/* Prints a string, ie a char[size] with representation in hexa.
	*/

void print_hash(unsigned char *out, int outLen);
	/* Prints the content of a digest.
	(Redondant with print_char, but more explicit name).
    outLen should be nb_of_comm / 8;
	*/

void ideal_gen(private_key *priv);
    /* Random ideal generation.
    Used to generate a random committment, which exponents are in |[-2 ; 2]|.
    WATCH OUT : Keep the exponents in this bound, 
        otherwise the answer of challenge 1 woudn't fit in int4 ! (5+2 = 7). 
    */

int8_t substract(int8_t e1, int8_t e2);
	/* Returns the substraction of the exponents packed in e1 by those in e2.
	Used to compute i1 // i2, where i1 and i2 are two ideals whose exponents are e1 and e2.
	WATCH OUT ! If the exponents of i1 are between -5 and 5, those of i2 should be 
		between -2 and 2, in order to avoid overflow.
	*/

void substract_tab(int N, int8_t e3[N], int8_t e1[N], int8_t e2[N]);
	/* Returns the substraction of the exponents packed in e1 by those in e2.
	Exponents are really private key, ie ideal factorised representation.
	Used to compute i1 // i2, where i1 and i2 are two ideals whose exponents are e1 and e2.
	WATCH OUT ! If the exponents of i1 are between -5 and 5, those of i2 should be 
		between -2 and 2, in order to avoid overflow.
	*/







#endif