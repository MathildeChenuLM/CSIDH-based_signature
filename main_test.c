
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <assert.h>

#include "u512.h"
#include "fp.h"
#include "mont.h"
#include "csidh.h"

#include <inttypes.h>

void u512_print(u512 const *x)
{
    for (size_t i = 63; i < 64; --i)
        printf("%02hhx", i[(unsigned char *) x->c]);
}

void fp_print(fp const *x)
{
    u512 y;
    fp_dec(&y, x);
    u512_print(&y);
}

static __inline__ uint64_t rdtsc(void)
{
    uint32_t hi, lo;
    __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
    return lo | (uint64_t) hi << 32;
}




int main()
{
//     clock_t t0, t1;
/*
    private_key priv_alice, priv_bob;
    public_key pub_alice, pub_bob;
    public_key shared_alice, shared_bob;*/
    int i;
//     int  j, Maxit = 100000;
// 
//     uint64_t c0, c1, costs[74];
//     proj P, Q,A;
//     const unsigned primes[num_primes] = {
//         3,   5,   7,  11,  13,  17,  19,  23,  29,  31,  37,  41,  43,  47,  53,  59,
//         61,  67,  71,  73,  79,  83,  89,  97, 101, 103, 107, 109, 113, 127, 131, 137,
//         139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227,
//         229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313,
//         317, 331, 337, 347, 349, 353, 359, 367, 373, 587,
//     };
//     
//     printf("\n");
    
//         proj K = { 0x8DF9C614576A396F, 0x3846CE6378B0BD97, 0x5DBCF39E66DF7EAE, 0x2289DA15AC75159B, 0xBD70753692CC160D, 0xB51960B3480DF8C8, 0x74768F44CC94D919, 
//             0x12DB57974B7A2DA4 , 
//             0xC8FC8DF598726F0A, 0x7B1BC81750A6AF95, 0x5D319E67C1E961B4, 0xB0AA7275301955F1, 
//             0x4A080672D9BA6C64, 0x97A5EF8A246EE77B, 0x6EA9E5D4383676A,    0x3496E2E117E0EC80};
        
//         proj A = {fp_0, fp_1}, A1 = {fp_0, fp_1};
//         proj A = {fp_0, fp_1};
//         proj A11 = {
//             0xE36BF4E92FE62B22, 0x37CD2014C2AB0BE0, 0x2A8FA4653A33837, 0x67A30F5DFFB40DA6, 0xEEFC14EC86062C73, 0x79F651EA05F0BF53, 0xB41A266BE4F4FFBF, 0x43DABFC6B2B76C63,
//             0x277DE6220AEAFA9E, 0x6B701D13B3AE55F5, 0xB7A21A3520FFA3BB, 0x15C5A1355F52F935, 0xE44F5A05D620B67D, 0x89D2A780AD403E21, 0x547FBA8A568AB41F, 0x6258172DB4393CEC
//         };
//         proj A24;
//         
//         proj P = {0x86821C366134D51B, 0xA8260866B54ED976, 0xE1BC77EC488D69C3, 0xFCADC8DC2D5B5786, 0xB6DDA8B6B7E292ED, 0x56F17E1CB52BF478, 0x30B715071C664DE1, 0x41E3F175A2E4B381 ,
//             0xC8FC8DF598726F0A, 0x7B1BC81750A6AF95, 0x5D319E67C1E961B4, 0xB0AA7275301955F1, 
//             0x4A080672D9BA6C64, 0x97A5EF8A246EE77B, 0x6EA9E5D4383676A,    0x3496E2E117E0EC80 };
//         proj Q = {0x4F247A466C87CE78, 0xFC575385CCB78D77, 0x44E41498431A696A, 0xEDEE55F1FB6F18B0, 0x70C5570D5E46AA8C, 0xF26C4DB628287B23, 0x4441AC61E8566A6, 0x4E8A99B4E1301565 ,
//             0xC8FC8DF598726F0A, 0x7B1BC81750A6AF95, 0x5D319E67C1E961B4, 0xB0AA7275301955F1, 
//             0x4A080672D9BA6C64, 0x97A5EF8A246EE77B, 0x6EA9E5D4383676A,    0x3496E2E117E0EC80 };   
//         
//         proj PQ = {0x7A5FF408B72034F8, 0x335105E276EB5E70, 0xE8840E0D5F3F2EC2, 0xF1A5D459CE0B65A8, 0x152901705C269E62, 0x822E25F88788DE73, 0xCBD9015C7AD2FDAA, 0xC13A0CE238FA8AB,
//             0xC8FC8DF598726F0A, 0x7B1BC81750A6AF95, 0x5D319E67C1E961B4, 0xB0AA7275301955F1, 
//             0x4A080672D9BA6C64, 0x97A5EF8A246EE77B, 0x6EA9E5D4383676A,    0x3496E2E117E0EC80 };     
//         proj R, S;
//         
//         fp_add3(&A24.z, &A11.z, &A11.z); //2C
//         fp_add3(&A24.x, &A11.x, &A24.z); //A + 2C
//         fp_add3(&A24.z, &A24.z, &A24.z); // 4C
//         
//         xDBLADD_zafin(&R, &S, &P, &Q, &PQ, &A24);
//         printf("R2 := Fp!0x");
//         for(i = 7; i>= 0; i--)
//             printf("%.16jx", R.x.x.c[i]);
//         printf("/ 0x");
//         for(i = 7; i>= 0; i--)
//             printf("%.16jx", R.z.x.c[i]);
//         printf(";\n");
//         
//         printf("S2 := Fp!0x");
//         for(i = 7; i>= 0; i--)
//             printf("%.16jx", S.x.x.c[i]);
//         printf("/ 0x");
//         for(i = 7; i>= 0; i--)
//             printf("%.16jx", S.z.x.c[i]);
//         printf(";\n");
//         
// //         xISOG(&A, &P, &K, 7);
//         xDBLADD(&R, &S, &P, &Q, &PQ, &A11);
//         printf("R1 := Fp!0x");
//         for(i = 7; i>= 0; i--)
//             printf("%.16jx", R.x.x.c[i]);
//         printf("/ 0x");
//         for(i = 7; i>= 0; i--)
//             printf("%.16jx", R.z.x.c[i]);
//         printf(";\n");
//         
//         printf("S1 := Fp!0x");
//         for(i = 7; i>= 0; i--)
//             printf("%.16jx", S.x.x.c[i]);
//         printf("/ 0x");
//         for(i = 7; i>= 0; i--)
//             printf("%.16jx", S.z.x.c[i]);
//         printf(";\n");
//         
//         
        fp a = {0x277DE6220AEAFA9E, 0x6B701D13B3AE55F5, 0xB7A21A3520FFA3BB, 0x15C5A1355F52F935, 0xE44F5A05D620B67D, 0x89D2A780AD403E21, 0x547FBA8A568AB41F, 0x6258172DB4393CEC};
        fp b;
        fp d = {0xDD69B39EF73BB83F, 0xC5A0EE2117ED83D9, 0x47EEF75E0EF4D6F8, 0x6515B6D692BC2FFD, 0x487DC5C137B75BA2, 0xBC43C11206A21EB6, 0xA71142E6DFD886D9, 0x30E776F5B5E498C4};
        fp_sqr_own(&b, &a);
        
        printf("b := Fp!0x");
        for(i = 7; i>= 0; i--)
            printf("%.16jx", b.x.c[i]);
        printf("/ Rm;\n");
        
        if (isSquare(&d)) 
            printf("\npositivo\n");
        else
            printf("Negativo\n");
//         for(i = 7; i>= 0; i--)
//             printf("%.16jx", S.z.x.c[i]);
//         printf(";\n");
        
        
//         printf("/Rm; \nEn := EllipticCurve(X^3 + (An/Cn)*X^2 + X);\n");
        
   /*     
        xISOG_2(&A1, &P, &K, 7);
        printf("An := Fp!0x");
        for(i = 7; i>= 0; i--)
            printf("%.16jx", A.x.x.c[i]);
        printf("/ Rm; \nCn := Fp!0x");
        for(i = 7; i>= 0; i--)
            printf("%.16jx", A.z.x.c[i]);
        printf("/Rm; \nEn := EllipticCurve(X^3 + (An/Cn)*X^2 + X);\n");*/
    
    
   /* 
    for(i = 0; i < 74;i++)
        costs[i] = 0;
    for(j = 0; j < Maxit ; j++)
            for(i = 0 ; i < 74 ; i++){
                c0 = rdtsc();
                xISOG(&A, &P, &Q, primes[i]);
                c1 = rdtsc();
                costs[i] += c1 - c0;
            }
    printf("Costs_xIsog= [\n");
    for(i = 0; i < 73;i++)
        printf("%"PRIu64",",(uint64_t)costs[i]/Maxit);
    printf("%"PRIu64"]\n",(uint64_t)costs[73]/Maxit);
    
    for(i = 0; i < 74;i++)
        costs[i] = 0;
    
    for(j = 0; j < Maxit ; j++)
        for(i = 0 ; i < 74 ; i++){
            c0 = rdtsc();
            xISOG_2(&A, &P, &Q, primes[i]);
            c1 = rdtsc();
            costs[i] += c1 - c0;
        }
        printf("\n\nCosts_xIsog_2 = [");
        for(i = 0; i < 73;i++)
            printf("%"PRIu64",",(uint64_t)costs[i]/Maxit);
        printf("%"PRIu64"]",(uint64_t)costs[73]/Maxit);
        */
        printf("\n\n");
}

