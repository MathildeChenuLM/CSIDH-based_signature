#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <assert.h>

#include "u512.h"
#include "fp.h"
#include "mont.h"
#include "mont_own.h"
#include "csidh.h"
#include "sign.h"


int main()
{	
	

	clock_t t0, t1;

	signature sig;
	signature_init(&sig);

	public_key base = {0};

	private_key priv_a;
	csidh_private(&priv_a);

	public_key pub_a;
	action(&pub_a, &base, &priv_a);

	printf("\n========================================\n");
    printf("\tOriginal");
    printf("\n========================================\n");

	printf("\nSigning... \n");
	t0 = clock();
	sign_original( &sig, &priv_a, &base);
	t1 = clock();
	printf("Signature done in %lf s.\n",   (double)(t1 - t0) / CLOCKS_PER_SEC );

	printf("\nVerifying... \n");
	t0 = clock();
	verify_original(&sig, &pub_a, &base);
	t1 = clock();
	printf("Verification done in %lf s.\n",   (double)(t1 - t0) / CLOCKS_PER_SEC );

	printf("Done !\n");


	printf("\n========================================\n");
    printf("\tX - Wing");
    printf("\n========================================\n");


	printf("\nSigning... \n");
	t0 = clock();
	sign_xwing( &sig, &priv_a, &base);
	t1 = clock();
	printf("Signature done in %lf s.\n",   (double)(t1 - t0) / CLOCKS_PER_SEC );

	printf("\nVerifying... \n");
	t0 = clock();
	verify_xwing(&sig, &pub_a, &base);
	t1 = clock();
	printf("Verification done in %lf s.\n",   (double)(t1 - t0) / CLOCKS_PER_SEC );

	printf("Done !\n");

	printf("\n========================================\n");
    printf("\tX-wing and torsion points");
    printf("\n========================================\n");

	printf("\nSigning... \n");
	t0 = clock();
	sign_torsion( &sig, &priv_a, &base);
	t1 = clock();
	printf("Signature done in %lf s.\n",   (double)(t1 - t0) / CLOCKS_PER_SEC );

	printf("\nVerifying... \n");
	t0 = clock();
	verify_torsion(&sig, &pub_a, &base);
	t1 = clock();
	printf("Verification done in %lf s.\n",   (double)(t1 - t0) / CLOCKS_PER_SEC );

	printf("Done !\n");

	printf("\n========================================\n");
    printf("\tMeyer-Reith");
    printf("\n========================================\n");

	printf("\nSigning... \n");
	t0 = clock();
	sign_MR( &sig, &priv_a, &base);
	t1 = clock();
	printf("Signature done in %lf s.\n",   (double)(t1 - t0) / CLOCKS_PER_SEC );

	printf("\nVerifying... \n");
	t0 = clock();
	verify_MR(&sig, &pub_a, &base);
	t1 = clock();
	printf("Verification done in %lf s.\n",   (double)(t1 - t0) / CLOCKS_PER_SEC );

	printf("Done !\n");

	printf("\n========================================\n");
    printf("\tMeyer-Reith - torsion");
    printf("\n========================================\n");

	printf("\nSigning... \n");
	t0 = clock();
	sign_MR_torsion( &sig, &priv_a, &base);
	t1 = clock();
	printf("Signature done in %lf s.\n",   (double)(t1 - t0) / CLOCKS_PER_SEC );

	printf("\nVerifying... \n");
	t0 = clock();
	verify_MR_torsion(&sig, &pub_a, &base);
	t1 = clock();
	printf("Verification done in %lf s.\n",   (double)(t1 - t0) / CLOCKS_PER_SEC );

	printf("Done !\n");

	signature_clear(&sig);
	return 0;
}
