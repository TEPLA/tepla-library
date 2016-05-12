#include <assert.h>

#include "rdtsc.h"

#include "../ec_bn254_lcl.h"

#define N 1000
#define M 100

#define t 80            // security Level for SHA
#define MAP_STR "LCIS"  // string for map to point

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

	free(s);
}

//============================================
//  楕円曲線の確認用
//============================================
void test_feature(const EC_GROUP ec)
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
//  楕円曲線の演算テスト
//============================================
void test_arithmetic_operation(const EC_GROUP ec)
{
	int i;
	unsigned long long int t1, t2;
	EC_POINT a, b, c, d;

	Element dx, dy;

	mpz_t scalar;

	//-------------------
	//  init
	//-------------------
	point_init(a, ec);
	point_init(b, ec);
	point_init(c, ec);
	point_init(d, ec);

	//-------------------
	//  random
	//-------------------
	point_random(a);

	assert( point_is_on_curve(a) );

	t1 = rdtsc();
	for(i=0;i<M;i++){ point_random(a); }
	t2 = rdtsc();

	printf("point random: %.2lf [clock]\n", (double)(t2-t1)/M);

	//-------------------
	//  add/dob
	//-------------------
	point_add(b, b, a);
	point_add(c, b, c);

	assert( point_cmp(b, a) == 0 );
	assert( point_cmp(c, b) == 0 );

	point_set_infinity(d);

	point_dob(d, d);

	assert( point_is_infinity(d) );

	point_set_str(a, "[6D2E4115FA177379A504A0EE4EF53767DE51C6364AAB69D4064529EC1FD047A 635B2C858AA4F4A3DB8AA17A588B037CAFFD36678F76E3F3369DFC90C6878C7,193E877C82EFCA81EC2815906630B837BBC6976CC8A7958E6A40D1B190FF2E5F E8A77E88AFCEE9F806DC15BF50EADD138320F1A5A87E78DDE86FA7A867300D]");
	point_set_str(b, "[1A8F5DAB09EE4290F95FE4C824C153E355D55B6CF94B998C6203FEC3D81377CF 15A19F2704C4BDBAAE39A5E26772A3E4E7EC7A9E205651F8822298766DE044FF,1C566EB3917F06B05E0A786BD8030CAFCCDB62864DD0E2A22A9B6817B310FD53 6A0927BB33EB263F45CAB921A20E67A1BD8A791D6EB0415AC92C9B1F74D16D1]");
	point_set_str(d, "[143D414F99AA18C844B331064C9DD66363EBA3D852250CBCF8C9D4B33E0C4C1C 225865D85EC7A34647CA55E026BD1FA201E0C4E8C66F7A43E69AF708F410A0FF,FACA1388C034CF614A72E06EE60DEDC4880CDBD368E5BEC2795130B266FFB9E 1681217E50705AB9A21FEB62E0BF9A5657EB27C3AED3323FE9C57058358735A9]");

	point_add(c, a, b);

	assert( point_cmp(c, d) == 0 );

	t1 = rdtsc();
	for(i=0;i<N;i++){ point_add(c, a, b); }
	t2 = rdtsc();

	printf("point add: %.2lf [clock]\n", (double)(t2-t1)/N);

	element_init(dx, ec->field);
	element_init(dy, ec->field);

	element_set_str(dx, "33F550F9A63EF53C786BF7BDFDAB1538CD76A3FCED3C9DBC3307DD4F354775A C814AE99C91C71845F0B51E4349520908E48C70181313D70C05F6ED24EC1F36");
	element_set_str(dy, "3766ED0DD7C988DB76770081A298DAA924D0E3279726F9B5504129AFA3E57B9 520CE8A563F88AF882AB99086BFDBBDCEF9DE65879AB234DFF5AAFD5BEE7E4F");

	point_set_xy(d, dx, dy);

	point_dob(c, a);

	assert( point_cmp(c, d) == 0 );

	t1 = rdtsc();
	for(i=0;i<N;i++){ point_dob(c, a); }
	t2 = rdtsc();

	printf("point dob: %.2lf [clock]\n", (double)(t2-t1)/N);

	element_clear(dx);
	element_clear(dy);

	//------------------
	// neg/sub
	//------------------
	point_neg(c, a);
	point_add(d, a, c);
	point_sub(b, a, a);

	assert( point_is_infinity(d) );
	assert( point_is_infinity(b) );

	//-------------------
	//  mul
	//-------------------
	mpz_init(scalar);

	mpz_set(scalar, ec->order);

	for(i=0;i<100;i++)
	{
		point_random(a);

		point_mul(b, scalar, a);

		ec_bn254_fp2_mul_end(c, scalar, a);

		assert( point_is_infinity(b) );
		assert( point_cmp(c, b) == 0 );
	}

	t1 = rdtsc();
	for(i=0;i<M;i++){ point_mul(b, scalar, a); }
	t2 = rdtsc();

	printf("point mul with endomorphism: %.2lf [clock]\n", (double)(t2-t1)/M);

	t1 = rdtsc();
	for(i=0;i<M;i++){ ec_bn254_fp2_mul(b, scalar, a); }
	t2 = rdtsc();

	printf("point mul with binary method: %.2lf [clock]\n", (double)(t2-t1)/M);

	mpz_clear(scalar);

	//-------------------
	//  clear
	//-------------------
	point_clear(a);
	point_clear(b);
	point_clear(c);
	point_clear(d);
}

//============================================
//  MAP to POINT テスト
//============================================
void test_map_to_point(const EC_GROUP ec)
{
	int i;
	unsigned long long int t1, t2;
	EC_POINT P, Q;

	point_init(P, ec);
	point_init(Q, ec);

	point_map_to_point(P, MAP_STR, sizeof(MAP_STR), t);

	assert( point_is_on_curve(P) );

	point_mul(Q, ec->order, P);

	assert( point_is_infinity(Q) );

	point_map_to_point(Q, MAP_STR, sizeof(MAP_STR), t);

	assert( point_cmp(Q, P) == 0 );

	t1 = rdtsc();
	for(i=0;i<M;i++){ point_map_to_point(P, MAP_STR, sizeof(MAP_STR), t); }
	t2 = rdtsc();

	printf("point map to point in 128 security: %.2lf [clock]\n", (double)(t2-t1)/M);

	point_clear(P);
	point_clear(Q);
}

//============================================
//  楕円曲線の入出力テスト
//============================================
void test_io(const EC_GROUP ec)
{
	int i;
	unsigned long long int t1, t2;
	EC_POINT P, Q, R;

	size_t osize;
	unsigned char os[130];

	char str[262];

	point_init(P, ec);
	point_init(Q, ec);
	point_init(R, ec);

	//---------------------
	//  octet string
	//---------------------
	point_set_infinity(R);

	point_to_oct(os, &osize, R);
	point_from_oct(Q, os, osize);

	assert( point_is_infinity(Q) );

	for(i=0;i<100;i++)
	{
		point_random(P);

		point_to_oct(os, &osize, P);
		point_from_oct(Q, os, osize);

		assert( point_cmp(P, Q) == 0 );
	}

	t1 = rdtsc();
	for(i=0;i<N;i++){ point_to_oct(os, &osize, P); }
	t2 = rdtsc();

	printf("point to octet string: %.2lf [clock]\n", (double)(t2-t1)/N);

	t1 = rdtsc();
	for(i=0;i<N;i++){ point_from_oct(Q, os, osize); }
	t2 = rdtsc();

	printf("point from octet string: %.2lf [clock]\n", (double)(t2-t1)/N);

	//---------------------
	//  string
	//---------------------
	point_set_infinity(R);

	point_get_str(str, R);
	point_set_str(Q, str);

	assert( point_is_infinity(Q) );

	for(i=0;i<100;i++)
	{
		point_get_str(str, P);
		point_set_str(Q, str);

		assert( point_cmp(P, Q) == 0 );
	}

	t1 = rdtsc();
	for(i=0;i<N;i++){ point_get_str(str, P); }
	t2 = rdtsc();

	printf("point get string: %.2lf [clock]\n", (double)(t2-t1)/N);

	t1 = rdtsc();
	for(i=0;i<N;i++){ point_set_str(Q, str); }
	t2 = rdtsc();

	printf("point set string: %.2lf [clock]\n", (double)(t2-t1)/N);

	point_clear(P);
	point_clear(Q);
	point_clear(R);
}

//============================================
//  メインプログラム
//============================================
int main(void)
{
	EC_GROUP ec;

	curve_init(ec, "ec_bn254_tw");

	test_feature(ec);

	test_arithmetic_operation(ec);

	test_map_to_point(ec);

	test_io(ec);

	curve_clear(ec);

	fprintf(stderr,"ok\n");

	return 0;
}
