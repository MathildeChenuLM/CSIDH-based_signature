#include "mont_own.h"


void xADD_zafin(proj *S, proj const *P, proj const *Q, proj const *PQ)
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
    fp_mul3(&S->z, &PQ->x, &d);
    S->x = c;
}//Cost 3M + 2S + 6A

void xDBL24(proj *Q, proj const *A, proj const *P)//xDBL:=function(X,Z,A24,C24)
 {  
    fp t0, t1;
//     proj A24;
//     fp_add3(&A24.z, &A->z, &A->z); //2C
//     fp_add3(&A24.x, &A->x, &A24.z); //A + 2C
//     fp_add3(&A24.z, &A24.z, &A24.z); // 4C
    
    fp_sub3(&t0,&P->x, &P->z);//         t0 := X-Z;           
    fp_add3(&t1,&P->x, &P->z);//         t1 := X+Z; 
    fp_sq1(&t0);//         t0 := t0^2; 
    fp_sq1(&t1);//         t1 := t1^2; 
    fp_mul3(&Q->z, &A->z, &t0);//         Z2 := C24*t0; 
    fp_mul3(&Q->x, &t1, &Q->z);//         X2 := Z2*t1;
    fp_sub2(&t1, &t0);//         t1 := t1-t0; 
    fp_mul3(&t0, &A->x, &t1);//         t0 := A24*t1;
    fp_add2(&Q->z, &t0);//         Z2 := Z2+t0;
    fp_mul2(&Q->z, &t1);//         Z2:= Z2*t1;
}//Cost 4M + 2S + 4A

void xDBL_own(proj *Q, proj const *A, proj const *P)//xDBL:=function(X,Z,A24,C24)
 {  
    fp t0, t1;
    proj A24;
    fp_add3(&A24.z, &A->z, &A->z); //2C
    fp_add3(&A24.x, &A->x, &A24.z); //A + 2C
    fp_add3(&A24.z, &A24.z, &A24.z); // 4C
    
    fp_sub3(&t0,&P->x, &P->z);//         t0 := X-Z;           
    fp_add3(&t1,&P->x, &P->z);//         t1 := X+Z; 
    fp_sq1(&t0);//         t0 := t0^2; 
    fp_sq1(&t1);//         t1 := t1^2; 
    fp_mul3(&Q->z, &A24.z, &t0);//         Z2 := C24*t0; 
    fp_mul3(&Q->x, &t1, &Q->z);//         X2 := Z2*t1;
    fp_sub2(&t1, &t0);//         t1 := t1-t0; 
    fp_mul3(&t0, &A24.x, &t1);//         t0 := A24*t1;
    fp_add2(&Q->z, &t0);//         Z2 := Z2+t0;
    fp_mul2(&Q->z, &t1);//         Z2:= Z2*t1;
}//Cost 4M + 2S + 7A


void xDBL_A24afin(proj *Q, proj const *A, proj const *P)//xDBL:=function(X,Z,A24,C24)
 {  
    fp t0, t1;
    
    fp_sub3(&t0,&P->x, &P->z);//         t0 := X-Z;           
    fp_add3(&t1,&P->x, &P->z);//         t1 := X+Z; 
    fp_sq1(&t0);//         t0 := t0^2; 
    fp_sq1(&t1);//         t1 := t1^2; 
    //C24 = 4
    //fp_mul3(&Q->z, &A24.z, &t0);//         Z2 := C24*t0; 
    fp_add3(&Q->z, &t0, &t0);//2*t0
    fp_add2(&Q->z, &Q->z);//4*t0
    //====
    fp_mul3(&Q->x, &t1, &Q->z);//         X2 := Z2*t1;
    fp_sub2(&t1, &t0);//         t1 := t1-t0; 
    // A24 := A + 2
    fp_mul3(&t0, &A->x, &t1);//         t0 := A24*t1;
    fp_add2(&Q->z, &t0);//         Z2 := Z2+t0;
    fp_mul2(&Q->z, &t1);//         Z2:= Z2*t1;
}//Cost 3M + 2S + 6A


// void xDBLADD_zafin(proj *R, proj *S, proj const *P, proj const *Q, proj const *PQ, proj const *A)
// {
//     fp a, b, c, d;
// 
//     fp_add3(&a, &Q->x, &Q->z);
//     fp_sub3(&b, &Q->x, &Q->z);
//     fp_add3(&c, &P->x, &P->z);
//     fp_sub3(&d, &P->x, &P->z);
//     fp_sq2(&R->x, &c);
//     fp_sq2(&S->x, &d);
//     fp_mul2(&c, &b);
//     fp_mul2(&d, &a);
//     fp_sub3(&b, &R->x, &S->x);
//     fp_add3(&a, &A->z, &A->z); /* multiplication by 2 */
//     fp_mul3(&R->z, &a, &S->x);
//     fp_add3(&S->x, &A->x, &a);
//     fp_add2(&R->z, &R->z); /* multiplication by 2 */
//     fp_mul2(&R->x, &R->z);
//     fp_mul2(&S->x, &b);
//     fp_sub3(&S->z, &c, &d);
//     fp_add2(&R->z, &S->x);
//     fp_add3(&S->x, &c, &d);
//     fp_mul2(&R->z, &b);
//     fp_sq2(&d, &S->z);
//     fp_sq1(&S->x);
//     //fp_mul3(&S->x, &PQ->z, &b);
//     fp_mul3(&S->z, &PQ->x, &d);
// }//7M + 4S + 11A;


// void xDBLADD_zafin(proj *R, proj *S, proj const *P, proj const *Q, proj const *PQ, proj const *A)
// {
//     fp  t0, t1, t2, t3;
//     
//     fp_sub3(&t2, &P->x, &P->z);     //    t2 := XP-ZP;      
//     fp_add3(&t3, &P->x, &P->z);     //    t3 := XP+ZP; //   
//     fp_sq2(&t0, &t2);               //    t0 := t2^2; //    
//     fp_sq2(&t1, &t3);               //    t1 := t3^2; //    
//     fp_mul3(&R->z, &A->z, &t0);     //    Z2 := C24*t0; //  
//     fp_mul3(&R->x, &R->z, &t1);     //    X2 := Z2*t1;//    
//     fp_sub2(&t1, &t0);              //    t1 := t1-t0; //   
//     fp_mul3(&t0, &A->x, &t1);       //    t0 := A24*t1;//   
//     fp_add2(&R->z, &t0);            //    Z2 := Z2+t0;//    
//     fp_mul2(&R->z, &t1);            //    Z2:= Z2*t1;    // 
//     // //xADD
//     fp_sub3(&S->x, &Q->x, &Q->z );  //    XQP:= XQ-ZQ;//    
//     fp_add3(&S->z, &Q->x, &Q->z );  //    ZQP:= XQ+ZQ;//    
//     fp_mul3(&t0, &S->x, &t3);       //    t0:= XQP*t3;//    
//     fp_mul3(&t1, &S->z, &t2);       //    t1:= ZQP*t2;//    
//     fp_sub3(&S->z, &t0, &t1);       //    ZQP:=t0-t1;//     
//     fp_add3(&S->x, &t0, &t1);       //    XQP:=t0+t1;//     
//     fp_sq1(&S->z);                  //    ZQP:=ZQP^2;//     
//     fp_sq1(&S->x);                  //    XQP:=XP^2;
//     fp_mul2(&S->x, &PQ->z);         //    XQP := XQP * zPQ;
//     fp_mul2(&S->z, &PQ->x);         //    ZQP:= xPQ*ZQP;
//                                     
// }   //cost 7M + 4S + 8A;

void xDBLADD_zafin(proj *R, proj *S, proj const *P, proj const *Q, proj const *PQ, proj const *A)
{
    xADD_zafin(S, P, Q, PQ);//3M + 2S + 6A
    xDBL24(R, A, P);//4M + 2S + 4A   
}//Cost 7M + 4S + 10A

void xDBLADD_A24afin(proj *R, proj *S, proj const *P, proj const *Q, proj const *PQ, proj const *A)
{
    xADD(S, P, Q, PQ);//3M + 2S + 6A
    xDBL_A24afin(R, A, P);//4M + 2S + 4A   
}//Cost 7M + 4S + 10A


void xDBLADD_afin(proj *R, proj *S, proj const *P, proj const *Q, proj const *PQ, proj const *A)
{
    xADD_zafin(S, P, Q, PQ);//3M + 2S + 6A
    xDBL_A24afin(R, A, P);//3M + 2S + 6A   
}//Cost 6M + 4S + 12A



void xMUL_zafin(proj *Q, proj const *A, proj const *P, u512 const *k)
{
    proj R = *P;
    const proj Pcopy = *P; /* in case Q = P */
    int i;
    proj A24;
    fp_add3(&A24.z, &A->z, &A->z); //2C
    fp_add3(&A24.x, &A->x, &A24.z); //A + 2C
    fp_add3(&A24.z, &A24.z, &A24.z); // 4C
    Q->x = fp_1;
    Q->z = fp_0;
    i = 512;
    while (--i && !u512_bit(k, i));
    do {
        bool bit = u512_bit(k, i);
        if (bit) { proj T = *Q; *Q = R; R = T; } /* not constant-time */
        //fp_cswap(&Q->x, &R.x, bit);
        //fp_cswap(&Q->z, &R.z, bit);
        xDBLADD_zafin(Q, &R, Q, &R, &Pcopy, &A24);
        if (bit) { proj T = *Q; *Q = R; R = T; } /* not constant-time */
        //fp_cswap(&Q->x, &R.x, bit);
        //fp_cswap(&Q->z, &R.z, bit);
    } while (i--);

}


void xMUL_afin(proj *Q, proj const *A, proj const *P, u512 const *k)
{
    proj R = *P;
    const proj Pcopy = *P; /* in case Q = P */
    int i;
    proj A24;
    fp_add3(&A24.z, &A->z, &A->z); //2
    fp_add3(&A24.x, &A->x, &A24.z); //A + 2
    A24.z = fp_1;
    Q->x = fp_1;
    Q->z = fp_0;
    i = 512;
    while (--i && !u512_bit(k, i));
    do {
        bool bit = u512_bit(k, i);
        if (bit) { proj T = *Q; *Q = R; R = T; } /* not constant-time */
        //fp_cswap(&Q->x, &R.x, bit);
        //fp_cswap(&Q->z, &R.z, bit);
        xDBLADD_afin(Q, &R, Q, &R, &Pcopy, &A24);
        if (bit) { proj T = *Q; *Q = R; R = T; } /* not constant-time */
        //fp_cswap(&Q->x, &R.x, bit);
        //fp_cswap(&Q->z, &R.z, bit);
    } while (i--);

}

void xMUL_A24afin(proj *Q, proj const *A, proj const *P, u512 const *k)
{
    proj R = *P;
    const proj Pcopy = *P; /* in case Q = P */
    int i;
    proj A24;
    fp_add3(&A24.z, &A->z, &A->z); //2
    fp_add3(&A24.x, &A->x, &A24.z); //A + 2
    //A24.z = fp_1;
    Q->x = fp_1;
    Q->z = fp_0;
    i = 512;
    while (--i && !u512_bit(k, i));
    do {
        bool bit = u512_bit(k, i);
        if (bit) { proj T = *Q; *Q = R; R = T; } /* not constant-time */
        //fp_cswap(&Q->x, &R.x, bit);
        //fp_cswap(&Q->z, &R.z, bit);
        xDBLADD_A24afin(Q, &R, Q, &R, &Pcopy, &A24);
        if (bit) { proj T = *Q; *Q = R; R = T; } /* not constant-time */
        //fp_cswap(&Q->x, &R.x, bit);
        //fp_cswap(&Q->z, &R.z, bit);
    } while (i--);

}



void cross_2(proj const X1, proj const X2, fp *s, fp *tx, fp *tz)
{
    fp t0;
    fp_add3(&t0, &X1.x, &X1.z);// t0 := x1 + z1;
    fp_add3(s, &X2.x, &X2.z);// s := x2 + z2;
    fp_mul2(s, &t0);// s := s * t0;
    fp_mul3(tx, &X1.x, &X2.x);// tx := x1*x2;
    fp_mul3(tz, &X1.z, &X2.z);// tz := z1 * z2;
    fp_sub2(s, tx);// s := s - tx;
    fp_sub2(s, tz);// s := s - tz;
}//cost: 3M + 4A

void x_wing(fp *sx1, fp *tx1, fp *sz1, fp *tz1, fp sx2, fp tx2, fp sz2, fp tz2)
{
    fp t0;
    //     sx := sx1 * tz2 + sx2 * tz1; //2M + 1A
    fp_mul2(sx1,&tz2 );
    fp_mul3(&t0, &sx2, tz1);
    fp_add2(sx1, &t0);
    //     sz := sz1 * tx2 + sz2 * tx1; //2M + 1A
    fp_mul2(sz1,&tx2 );
    fp_mul3(&t0, &sz2, tx1);
    fp_add2(sz1, &t0);
    
    fp_mul2(tx1, &tx2);//     tx := tx1 * tx2;//M
    fp_mul2(tz1, &tz2);//     tz := tz1 * tz2;//M
}//cost: 6M + 2A

void part_ev(fp Pplus, fp Pminus, proj P, proj *Q)
{
    fp Xplus, Xminus, t0;
    fp_add3(&Xplus, &P.x, &P.z);//     Xplus := X1 + Z1;
    fp_sub3(&Xminus, &P.x, &P.z);//     Xminus := X1 - Z1;
    //     
    fp_mul2(&Xplus, &Pminus);//     Xplus := Pminus * Xplus; //(X - Z)(Xi + Zi);
    fp_mul2(&Xminus, &Pplus);//     Xminus := Pplus * Xminus; //(X + Z)(Xi - Zi);
    fp_add3(&t0, &Xplus, &Xminus);//     xQ := Xplus + Xminus;
    fp_mul2(&Q->x, &t0);//     xQ := xQ * XQ;
    fp_sub3(&t0, &Xplus, &Xminus);//     zQ := Xplus - Xminus;
    fp_mul2(&Q->z, &t0);//     zQ := zQ * ZQ;
}//4M + 4A

void xISOG_2(proj *A, proj *P, proj const *K, uint64_t d)
{
// generic_d_isog_cross2_AC := function(X, Z, A, C, d , XP, ZP)
//     //Require   K := (X:Z) which is the generator of the kernel
//     //          Curve projective constant A := (Ax:Az)
//     //          Isogeny degree d := 2*l + 1
//     //          Point P to be evaluate
//     //Ensure    Isogenous Curve projective constant a':=(An : Cn)
//     //          Image of P := Q = (xQ : zQ)
        int l, i, mask0;
        fp Pplus, Pminus;
        fp sx1, sx2, sz1, tx1, tx2, tz1, tz2;
        proj Q, P2;
        proj P1 = {K->x, K->z};
        Q.x = fp_1;
        Q.z = fp_1;
        if(d != 3){
            
                        l = d >> 1; //     l := d div 2;
                        mask0 = l & 0x01;//     flag := l mod 2;

                        fp_add3(&Pplus, &P->x, &P->z);//     Pplus := XP + ZP;
                        fp_sub3(&Pminus, &P->x, &P->z);//     Pminus := XP - ZP;
                        //For evaluation
                        part_ev(Pplus, Pminus, P1, &Q);//     xQ, zQ := part_ev(Pplus, Pminus, X1, Z1, 1, 1);                        

                        xDBL(&P2, A, &P1);//     X2, Z2 := xDBL(X1, Z1, A24p, C24);
                //     //For evaluation
                        part_ev(Pplus, Pminus, P2, &Q);//     xQ, zQ := part_ev(Pplus, Pminus, X2, Z2, xQ, zQ);
                        
                        cross_2( P1, P2, &sx1, &tx1,  &tz1);//     sx1, tx1, tz1 := Cross_2(X1, Z1, X2, Z2);
                        sz1 = sx1;//     sz1 := sx1;
                        if (l >= 4)    //     if l ge 4 then
                        {
                            for(i = 3; i <= (l - mask0); i +=2)//         for i := 3 to (l-flag) by 2 do
                            {
                                //             //Kernel Points
                                //             //[i]P
                                xADD(&P1, K, &P2, &P1);//             X1, Z1 := xADD_proj(X, Z, X2, Z2, X1, Z1);
                                //             //[i+1]P
                                xADD(&P2, K, &P1, &P2);//             X2, Z2 := xADD_proj(X, Z, X1, Z1, X2, Z2);
                                //             //Evaluation of P
                                part_ev(Pplus, Pminus, P1, &Q);//             xQ, zQ := part_ev(Pplus, Pminus, X1, Z1, xQ, zQ);
                                part_ev(Pplus, Pminus, P2, &Q);//             xQ, zQ := part_ev(Pplus, Pminus, X2, Z2, xQ, zQ);
                                //             //New child
                                cross_2( P1, P2, &sx2, &tx2, &tz2);//             sx2, tx2, tz2 := Cross_2(X1, Z1, X2, Z2);
                                //             //Merging results
                                x_wing(&sx1, &tx1, &sz1, &tz1, sx2, tx2, sx2, tz2);//             sx1, tx1, sz1, tz1 := X_Wing(sx1, tx1, sz1, tz1, sx2, tx2, sz2, tz2);
                            }//         end for;
                        }//     end if;
                // 
                        if (mask0 == 1)//     if flag eq 1 then
                        {
                            xADD(&P1, K, &P2, &P1);//         X1, Z1 := xADD_proj(X, Z, X2, Z2, X1, Z1);
                            //         //For evaluation
                            part_ev(Pplus, Pminus, P1, &Q);//         xQ, zQ := part_ev(Pplus, Pminus, X1, Z1, xQ, zQ);
                            //         // End evaluation
                            x_wing(&sx1, &tx1, &sz1, &tz1, P1.x, P1.x, P1.z, P1.z);//         sx1, tx1, sz1, tz1 := X_Wing(sx1, tx1, sz1, tz1, X1, X1, Z1, Z1);
                        }//     end if;
                        
                        fp_mul2(&sx1, &tx1);//sx := sx * tx;
                        fp_mul2(&sz1, &tz1);//sz := sz * tz;
                        fp_sub2(&sz1, &sx1);//sz := sz - sx; //(sz*tz - sx*tx)
                        fp_add3(&sx1, &sz1, &sz1); //sx1 := sz + sz;// 2*C*(sz*tz - sx*tx)
                        fp_add3(&sz1, &sx1, &sx1);//4*(sz*tz - sx*tx)
                        fp_add2(&sz1, &sx1);//6*(sz*tz - sx*tx)                        
                        fp_mul2(&sz1, &A->z);//sz1 := sz1 * C;//6*C*(sz*tz - sx*tx)
            
                        fp_mul2(&A->x, &tx1);//         An := A * tx; //M
                        fp_mul2(&A->x, &tz1);//         An := An * tz; //M
                        fp_add2(&A->x, &sz1);//         An := An + delta; //A
                        fp_mul2(&A->x, &tx1);//         An := An * tx; //M
                        
                        fp_sq2(&tx1, &tz1);//         t0 := tz^2;//S
                        fp_mul2(&A->z, &tx1);//         Cn := C * t0;//M
                        fp_mul2(&A->z, &tz1);//         Cn := Cn * tz;//M
                        
//                         //Final step for evaluation
                        fp_sq1(&Q.x);//     xQ := xQ^2;
                        fp_sq1(&Q.z);//     zQ := zQ^2;
                        fp_mul2(&P->x, &Q.x);//     xQ := xQ * XP;
                        fp_mul2(&P->z, &Q.z);//     zQ := zQ * ZP;
        }
        else
        {
                    //For d = 3
                    //get_3_isog CLN version 3S + 3M 
                    fp_sq2(&sx1, &K->x);//         t0:=X3^2;
                    fp_add3(&sz1, &sx1, &sx1);//         t1:=t0+t0;
                    fp_add3(&sx1, &sx1, &sz1);//         t0:=t0+t1;
                    fp_sq2(&sz1, &K->z);//         t1:=Z3^2;
                    fp_sq2(&A->x, &sz1);//         A:=t1^2;
                    fp_add2(&sz1, &sz1);//         t1:=t1+t1;
                    fp_add3(&A->z, &sz1, &sz1);//         C:=t1+t1;
                    fp_sub3(&sz1, &sx1, &sz1);//         t1:=t0-t1;
                    fp_mul2(&sz1, &sx1);//         t1:=t1*t0;
                    fp_sub2(&A->x, &sz1);//         A:=A-t1;  
                    fp_sub2(&A->x, &sz1);//         A:=A-t1;
                    fp_sub2(&A->x, &sz1);//         A:=A-t1;
                    fp_mul3(&sz1, &K->x, &K->z);//         t1:=X3*Z3;
                    fp_mul2(&A->z, &sz1 );//         C:=C*t1;
                    
                    //eval_3_isog CLN Version 6M + 2S
                    fp_mul3(&sx1, &K->x, &P->x);//         t0:=X3*X;       
                    fp_mul3(&sz1, &K->z, &P->x);//         t1:=Z3*X;       
                    fp_mul3(&tx1, &K->z, &P->z);//         t2:=Z3*Z;
                    fp_sub2(&sx1, &tx1);//         t0:=t0-t2;              
                    fp_mul3(&tx1, &P->z, &K->x);//         t2:=Z*X3;               
                    fp_sub2(&sz1, &tx1);//         t1:=t1-t2;
                    fp_sq1(&sx1);//         t0:=t0^2;       
                    fp_sq1(&sz1);//         t1:=t1^2;       
                    fp_mul2(&P->x, &sx1);//         X:=X*t0;        
                    fp_mul2(&P->z, &sz1);//         Z:=Z*t1;   
                    
        }

}
