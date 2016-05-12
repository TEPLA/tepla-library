#include <assert.h>

#include "rdtsc.h"

#include "../ec_bn254_lcl.h"

#define N 10000

//============================================
//   Feature of Field
//============================================
void test_feature(Field f)
{
    fprintf(stdout, "---\n");
    fprintf(stdout, "FieldType: %s\n", field_get_name(f));
    gmp_fprintf(stdout, "   order of field: %Zx\n", f->order);
    fprintf(stdout, "---\n");
}

//============================================
//   test for arithmetic operations
//============================================
void test_arithmetic_operation_beuchat(Field f)
{
    int i;
    unsigned long long int t1, t2;
    Element a, b, c, d;

    char loop[] = "1000";

    mpz_t exp;

    //--------------------
    //  init
    //--------------------
    element_init(a, f);
    element_init(b, f);
    element_init(c, f);
    element_init(d, f);

    //--------------------
    //  add
    //--------------------
    element_set_str(a, "214C703B4E25BB6E6729B1F388989CF38F33FB3685956D7F39DB604FCC4FB79C");
    element_set_str(b, "1AC62B775FF7059B6E5A03A65097C6C73BD85A016CD6452D3E99C936216BBEEA");
    element_set_str(d, "18A1A0AE10DBB14B870D9B01524B61794E0A15F7DA6BB2AAF9F4C985EDBB7685");

    element_add(c, a, b);

    assert(element_cmp(c, d) == 0);

    t1 = rdtsc();
    for (i = 0; i < N; i++) {
        element_add(c, a, b);
    }
    t2 = rdtsc();

    printf("element add: %.2lf [clock]\n", (double)(t2 - t1) / N);

    //--------------------
    //  sub
    //--------------------
    element_set(d, c);
    element_sub(c, c, d);

    assert(element_is_zero(c));

    //--------------------
    //  mul
    //--------------------
    element_mul(c, a, b);
    element_set_str(d, "EFAE50A68D3DC9F7A6D75A1D0CE556C074FCD360B4B58BD2AA489722F7378AD");

    assert(element_cmp(c, d) == 0);

    t1 = rdtsc();
    for (i = 0; i < N; i++) {
        element_mul(c, a, b);
    }
    t2 = rdtsc();

    printf("element mul: %.2lf [clock]\n", (double)(t2 - t1) / N);

    //--------------------
    //  random
    //--------------------
    element_random(a);
    element_random(b);

    //--------------------
    //  inv
    //--------------------
    element_mul(c, a, b);
    element_inv(b, b);
    element_mul(c, c, b);
    element_inv(d, a);
    element_mul(d, a, d);

    assert(element_cmp(c, a) == 0);
    assert(element_is_one(d));

    t1 = rdtsc();
    for (i = 0; i < N; i++) {
        element_inv(b, a);
    }
    t2 = rdtsc();

    printf("element inv: %.2lf [clock]\n", (double)(t2 - t1) / N);

    //--------------------
    //  pow
    //--------------------

    mpz_init_set_str(exp, loop, 10);

    element_set_one(b);

    for (i = 0; i < atoi(loop); i++) {
        element_mul(b, b, a);
    }

    element_pow(c, a, exp);

    assert(element_cmp(b, c) == 0);

    mpz_set(exp, f->order);

    for (i = 0; i < 100; i++)
    {
        element_random(a);
        element_pow(b, a, exp);

        assert(element_cmp(a, b) == 0);
    }

    t1 = rdtsc();
    for (i = 0; i < N; i++) {
        element_pow(b, a, exp);
    }
    t2 = rdtsc();

    printf("element pow with order: %.2lf [clock]\n", (double)(t2 - t1) / N);

    mpz_clear(exp);

    //--------------------
    //  clear
    //--------------------
    element_clear(a);
    element_clear(b);
    element_clear(c);
    element_clear(d);
}

void test_arithmetic_operation_aranha(Field f)
{
    int i;
    unsigned long long int t1, t2;
    Element a, b, c, d;

    char loop[] = "1000";

    mpz_t exp;

    //--------------------
    //  init
    //--------------------
    element_init(a, f);
    element_init(b, f);
    element_init(c, f);
    element_init(d, f);

    //--------------------
    //  add
    //--------------------
    element_set_str(a, "214C703B4E25BB6E6729B1F388989CF38F33FB3685956D7F39DB604FCC4FB79C");
    element_set_str(b, "1AC62B775FF7059B6E5A03A65097C6C73BD85A016CD6452D3E99C936216BBEEA");
    element_set_str(d, "3C129BB2AE1CC109D583B599D93063BACB0C5537F26BB2AC78752985EDBB7686");

    element_add(c, a, b);

    assert(element_cmp(c, d) == 0);

    t1 = rdtsc();
    for (i = 0; i < N; i++) {
        element_add(c, a, b);
    }
    t2 = rdtsc();

    printf("element add: %.2lf [clock]\n", (double)(t2 - t1) / N);

    //--------------------
    //  sub
    //--------------------
    element_set(d, c);
    element_sub(c, c, d);

    assert(element_is_zero(c));

    //--------------------
    //  mul
    //--------------------
    element_mul(c, a, b);
    element_set_str(d, "C100D8F6E60A9C9CEA9C26A22FDC75486BCEBC00E12C7B67307DB79358EC3AA");
    assert(element_cmp(c, d) == 0);

    t1 = rdtsc();
    for (i = 0; i < N; i++) {
        element_mul(c, a, b);
    }
    t2 = rdtsc();

    printf("element mul: %.2lf [clock]\n", (double)(t2 - t1) / N);

    //--------------------
    //  random
    //--------------------
    element_random(a);
    element_random(b);

    //--------------------
    //  inv
    //--------------------
    element_mul(c, a, b);
    element_inv(b, b);
    element_mul(c, c, b);
    element_inv(d, a);
    element_mul(d, a, d);

    assert(element_cmp(c, a) == 0);
    assert(element_is_one(d));

    t1 = rdtsc();
    for (i = 0; i < N; i++) {
        element_inv(b, a);
    }
    t2 = rdtsc();

    printf("element inv: %.2lf [clock]\n", (double)(t2 - t1) / N);

    //--------------------
    //  pow
    //--------------------

    mpz_init_set_str(exp, loop, 10);

    element_set_one(b);

    for (i = 0; i < atoi(loop); i++) {
        element_mul(b, b, a);
    }

    element_pow(c, a, exp);

    assert(element_cmp(b, c) == 0);

    mpz_set(exp, f->order);

    for (i = 0; i < 100; i++)
    {
        element_random(a);
        element_pow(b, a, exp);

        assert(element_cmp(a, b) == 0);
    }

    t1 = rdtsc();
    for (i = 0; i < N; i++) {
        element_pow(b, a, exp);
    }
    t2 = rdtsc();

    printf("element pow with order: %.2lf [clock]\n", (double)(t2 - t1) / N);

    mpz_clear(exp);

    //--------------------
    //  clear
    //--------------------
    element_clear(a);
    element_clear(b);
    element_clear(c);
    element_clear(d);
}
//============================================
//   test for sqrt
//============================================
void test_sqrt(Field f)
{
    int i;
    unsigned long long int t1, t2;
    Element a, b, c, d;

    element_init(a, f);
    element_init(b, f);
    element_init(c, f);
    element_init(d, f);

    for (i = 0; i < 1000; i++)
    {
        element_random(a);
        element_sqr(b, a);

        assert(element_is_sqr(b));

        element_sqrt(c, b);
        element_sqr(d, c);

        assert(element_cmp(d, b) == 0);
    }

    t1 = rdtsc();
    for (i = 0; i < N; i++) {
        element_is_sqr(b);
    }
    t2 = rdtsc();


    printf("element is sqr: %.2lf [clock]\n", (double)(t2 - t1) / N);

    t1 = rdtsc();
    for (i = 0; i < N; i++) {
        element_sqrt(c, b);
    }
    t2 = rdtsc();

    printf("element sqrt: %.2lf [clock]\n", (double)(t2 - t1) / N);

    element_clear(a);
    element_clear(b);
    element_clear(c);
    element_clear(d);
}

//============================================
//  i/o test
//============================================
void test_io(Field f)
{
    int i;
    unsigned long long int t1, t2;

    char a_str[256];

    size_t blen;
    unsigned char b_str[256];

    Element a, b, c;

    element_init(a, f);
    element_init(b, f);
    element_init(c, f);

    for (i = 0; i < 1000; i++)
    {
        element_random(a);

        element_get_str(a_str, a);
        element_set_str(c, a_str);

        assert(element_cmp(a, c) == 0);
    }

    t1 = rdtsc();
    for (i = 0; i < N; i++) {
        element_get_str(a_str, a);
    }
    t2 = rdtsc();

    printf("element get string: %.2lf [clock]\n", (double)(t2 - t1) / N);

    t1 = rdtsc();
    for (i = 0; i < N; i++) {
        element_set_str(c, a_str);
    }
    t2 = rdtsc();

    printf("element set string: %.2lf [clock]\n", (double)(t2 - t1) / N);

    for (i = 0; i < 10000; i++)
    {
        element_random(b);

        element_to_oct(b_str, &blen, b);
        element_from_oct(c, b_str, blen);

        assert(element_cmp(b, c) == 0);
    }

    t1 = rdtsc();
    for (i = 0; i < N; i++) {
        element_to_oct(b_str, &blen, b);
    }
    t2 = rdtsc();

    printf("element to octet string: %.2lf [clock]\n", (double)(t2 - t1) / N);

    t1 = rdtsc();
    for (i = 0; i < N; i++) {
        element_from_oct(c, b_str, blen);
    }
    t2 = rdtsc();

    printf("element from octet string: %.2lf [clock]\n", (double)(t2 - t1) / N);

    element_clear(a);
    element_clear(b);
    element_clear(c);
}

//============================================
// main program
//============================================
int main(void)
{
    Field fa, fb;

    // test for beuchat's methods
    field_init(fa, "bn254_fpa");
    test_feature(fa);
    test_arithmetic_operation_beuchat(fa);
    test_sqrt(fa);
    test_io(fa);

    // test for aranha's methods
    field_init(fb, "bn254_fpb");
    test_feature(fb);
    test_arithmetic_operation_aranha(fb);
    test_sqrt(fb);
    test_io(fb);

    field_clear(fa);
    field_clear(fb);

    fprintf(stderr, "ok\n");

    return 0;
}
