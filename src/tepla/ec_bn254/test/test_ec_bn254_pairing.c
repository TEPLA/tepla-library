#include <assert.h>

#include "../ec_bn254_lcl.h"

#include "rdtsc.h"

#define N 100

//============================================
//  Element と Point の表示用
//============================================
void print_element(const Element x)
{
	int len = element_get_str_length(x);

	char *s = (char*)malloc(sizeof(char)*len);

	element_get_str(s, x);

	printf("element: %s\n", s);

	free(s);
}

void print_point(const EC_POINT x)
{
	int len = point_get_str_length(x);

	char *s = (char*)malloc(sizeof(char)*len);

	point_get_str(s, x);

	printf("point: %s\n", s);
}

//============================================
//  有限体の確認用
//============================================
void print_field_feature(const Field f)
{
	fprintf(stdout, "---\n");
	fprintf(stdout, "Field Type: %s\n", field_get_name(f));
	gmp_fprintf(stdout, "   order of field: %Zx^%d\n", *field_get_char(f), field_get_degree(f));
	fprintf(stdout, "---\n");
}

//============================================
//  楕円曲線の確認用
//============================================
void print_ec_feature(const EC_GROUP ec)
{
	fprintf(stdout, "---\n");
	fprintf(stdout, "Elliptic Curve Type: %s\n", curve_get_name(ec));
	fprintf(stdout, "   Y^2 = X^3 + aX + b\n");
	fprintf(stdout, "   a: ");
	print_element(ec->a);
	fprintf(stdout, "   b: ");
	print_element(ec->b);
	gmp_fprintf(stdout, "   field order: %Zx^%d\n", *field_get_char(ec->field), field_get_degree(ec->field));
	fprintf(stdout, "   generator of curve: ");
	print_point(ec->generator);
	gmp_fprintf(stdout, "   order: %Zx\n", ec->order);
	gmp_fprintf(stdout, "   trace: %Zx\n", ec->trace);
	gmp_fprintf(stdout, "   cofactor: %Zx\n", ec->cofactor);
	fprintf(stdout, "---\n");
}

//============================================
//  ペアリングの確認用
//============================================
void test_feature(const EC_PAIRING p)
{
	fprintf(stdout, "---\n");
	fprintf(stdout, "Pairing Type: ECBN254\n");
	fprintf(stdout, "Group1: E1[r](Fp)\n");
	print_ec_feature(p->g1);
	fprintf(stdout, "Group2: Etw[r](Fp2)\n");
	print_ec_feature(p->g2);
	fprintf(stdout, "Group3: Fp12\n");
	print_field_feature(p->g3);
	fprintf(stdout, "---\n");
}

void test_pairing(const EC_PAIRING p)
{
	int i;
	unsigned long long int t1, t2;
	EC_POINT P, Q, R, S;
	Element d, e, f, g, h;

	gmp_randstate_t state;
	mpz_t a, b, c, order;

	//-------------------
	//  init
	//-------------------
	point_init(P, p->g1);
	point_init(R, p->g1);

	point_init(Q, p->g2);
	point_init(S, p->g2);

	element_init(d, p->g3);
	element_init(e, p->g3);
	element_init(f, p->g3);
	element_init(g, p->g3);
	element_init(h, p->g3);

	gmp_randinit_default(state);
	gmp_randseed_ui(state, (int)time(NULL));

	mpz_init(a);
	mpz_init(b);
	mpz_init(c);
	mpz_init(order);

	mpz_set(order, *pairing_get_order(p));

	for(i=0;i<10;i++)
	{
		//-------------------
		//  random
		//-------------------
		mpz_urandomm(a, state, order);
		mpz_urandomm(b, state, order);

		mpz_mul(c, a, b);
		mpz_mod(c, c, order);

		point_random(P);
		point_random(Q);

		point_mul(R, a, P);      // R = aP
		point_mul(S, b, Q);      // S = bQ

		pairing_map(g, P, Q, p);  // g = e(Q, P)

		element_pow(d, g, order);

		assert( element_is_one(d) );

		pairing_map(f, R, Q, p);  // f = e( Q, aP)
		pairing_map(h, P, S, p);  // h = e(bQ,  P)
		pairing_map(e, R, S, p);  // e = e(bQ, aP)

		element_pow(d, g, a);	 // d = e(Q, P)^a

		assert( element_cmp(d, f) == 0 );

		element_pow(d, g, b);    // d = e(Q, P)^b

		assert( element_cmp(d, h) == 0 );

		element_pow(d, g, c);    // d = e(Q, P)^c

		assert( element_cmp(d, e) == 0 );
	}

	t1 = rdtsc();
	for(i=0;i<N;i++){ pairing_map(g, P, Q, p); }
	t2 = rdtsc();

	printf("optimal ate pairing: %.5lf [clock]\n", (double)(t2-t1)/N);

	t1 = clock();
	for(i=0;i<N;i++){ p->pairing(g, Q, P, p); }
	t2 = clock();

	printf("optimal ate pairing: %.5lf [sec]\n", (double)(t2-t1)/(CLOCKS_PER_SEC)/N);

	//-------------------
	//  clear
	//-------------------
	point_clear(P);
	point_clear(R);
	point_clear(Q);
	point_clear(S);

	element_clear(d);
	element_clear(e);
	element_clear(f);
	element_clear(g);
	element_clear(h);

	mpz_clear(a);
	mpz_clear(b);
	mpz_clear(c);
	mpz_clear(order);

	gmp_randclear(state);
}

//============================================
//  メインプログラム
//============================================
int main(void)
{
	EC_PAIRING p;

	pairing_init(p, "ECBN254");

	test_feature(p);

	test_pairing(p);

	pairing_clear(p);

	fprintf(stderr,"ok\n");

	return 0;
}
