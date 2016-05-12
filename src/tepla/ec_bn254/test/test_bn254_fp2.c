#include <assert.h>

#include "rdtsc.h"

#include "../ec_bn254_lcl.h"

#define N 1000

//============================================
//   Feature of Field
//============================================
void test_feature(Field f)
{
    fprintf(stdout, "---\n");
    fprintf(stdout, "FieldType: %s\n", field_get_name(f));
    gmp_fprintf(stdout, "   characteristic: %Zx\n", *field_get_char(f));
    gmp_fprintf(stdout, "   order of field: %Zx^%d\n", *field_get_char(f), field_get_degree(f));
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

    char loop[] = "100";

    mpz_t e, exp;

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
    element_set_str(a, "1C12C39A2AD14054EDC9EE504301127AFFEEAADC59A78B50FCFFED87AC6EB8BF 20E1A922384561EA82602CD664D85D442DAC5D391E142ABB3CFEC2A095C22DF9");
    element_set_str(b, "F1B91250A124F268B8239185B23B31EB25179A11A9A0398E61B701F7D4F7265 20D206C5F7D007EDBA34A4B041622289D64F04CA28CEAC490619585AA14F7B2F");
    element_set_str(d, "7BD59BA97A27FBD2AD60CD0173FC358353DE53D5C418EE8649AFDA729BE2B23 1E42B4E392D45A19EE1EB6EE1F557D8C86F922C32EE2D702C497BAFB3711A927");

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
    element_set_str(d, "1D0562FF0AB317FFDE555320A7072D2B29C07077E08996CE5F093BB8E4200B2C 9B04361A24DC7F37C8BD09A7C51A9D8577168AD021BF2B4AC3D67552F481B1A");

    assert(element_cmp(c, d) == 0);

    t1 = rdtsc();
    for (i = 0; i < N; i++) {
        element_mul(c, a, b);
    }
    t2 = rdtsc();

    printf("element mul: %.2lf [clock]\n", (double)(t2 - t1) / N);

    mpz_init_set_str(e, "1B45F16C848B9C476C1D2FF1FD60A0D0C19BBA6F3ECE3CF6C5FCE4FAB7CAD4FF", 16);

    element_pow(c, a, e);
    element_set_str(d, "B40190CE812CB4F668A839952128D19B1748F3BB19E902480D089AF9053A6D2 19DA59F09C3C20472C3BD19A4FC95BCAF266B9D1539AAD23E3C67C4F3A7CA51D");

    assert(element_cmp(c, d) == 0);

    mpz_clear(e);

    //--------------------
    //  sqr
    //--------------------
    element_sqr(c, a);
    element_mul(d, a, a);

    assert(element_cmp(c, d) == 0);

    t1 = rdtsc();
    for (i = 0; i < N; i++) {
        element_sqr(c, a);
    }
    t2 = rdtsc();

    printf("element sqr: %.2lf [clock]\n", (double)(t2 - t1) / N);

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

        assert(element_cmp(b, a) == 0);
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

    char loop[] = "100";

    mpz_t e, exp;

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
    element_set_str(a, "1C12C39A2AD14054EDC9EE504301127AFFEEAADC59A78B50FCFFED87AC6EB8BF 20E1A922384561EA82602CD664D85D442DAC5D391E142ABB3CFEC2A095C22DF9");
    element_set_str(b, "F1B91250A124F268B8239185B23B31EB25179A11A9A0398E61B701F7D4F7265 20D206C5F7D007EDBA34A4B041622289D64F04CA28CEAC490619585AA14F7B2F");
    element_set_str(d, "60AF03CF4E38F79BF17D9E89E24C591511F247D74418ED63C1B5DA729BE2B11 1C904B65F01569D682608406A63A7FC5A2DA620346E2D6F09C181AFB3711A915");

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
    element_set_str(d, "120BDA59946F3C2A90C8D47A712545C17CDF2B23A217F74C3A3C6D0589275F02 170230E6B3E171E932E1EDA09803F6CE20D473BB94BA664A4C3E4EDFB6E9246D");

    assert(element_cmp(c, d) == 0);

    t1 = rdtsc();
    for (i = 0; i < N; i++) {
        element_mul(c, a, b);
    }
    t2 = rdtsc();

    printf("element mul: %.2lf [clock]\n", (double)(t2 - t1) / N);


    mpz_init_set_str(e, "1B45F16C848B9C476C1D2FF1FD60A0D0C19BBA6F3ECE3CF6C5FCE4FAB7CAD4FF", 16);
    element_pow(c, a, e);

    element_set_str(d, "184232E17A5CA89E69C789E94FADA2D941147C49D412C36CA8EF06A92636FDD2 B9AEEBC593CB9131C3F182C5CF70CE5647DBEAE59BFD8E12D5978372AB2A5A6");

    assert(element_cmp(c, d) == 0);

    mpz_clear(e);

    //--------------------
    //  sqr
    //--------------------
    element_sqr(c, a);
    element_mul(d, a, a);

    assert(element_cmp(c, d) == 0);

    t1 = rdtsc();
    for (i = 0; i < N; i++) {
        element_sqr(c, a);
    }
    t2 = rdtsc();

    printf("element sqr: %.2lf [clock]\n", (double)(t2 - t1) / N);

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

        assert(element_cmp(b, a) == 0);
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

    for (i = 0; i < 100; i++)
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
//   Frobenius Map \phi_p
//============================================
void test_frob(Field f)
{
    int i;
    unsigned long long int t1, t2;
    mpz_t p;
    Element a, b, c;

    mpz_init_set(p, *field_get_char(f));

    element_init(a, f);
    element_init(b, f);
    element_init(c, f);

    for (i = 0; i < 100; i++)
    {
        element_random(a);
        element_pow(b, a, p);

        bn254_fp2_frob_p(c, a);

        assert(element_cmp(b, c) == 0);
    }

    t1 = rdtsc();
    for (i = 0; i < N; i++) {
        bn254_fp2_frob_p(c, a);
    }
    t2 = rdtsc();

    printf("element frob: %.2lf [clock]\n", (double)(t2 - t1) / N);

    mpz_clear(p);

    element_clear(a);
    element_clear(b);
    element_clear(c);
}

//============================================
//   i/o test
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

    for (i = 0; i < 1000; i++)
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
//  main program
//============================================
int main(void)
{
    Field fa, fb;

    // test for beuchat's methods
    field_init(fa, "bn254_fp2a");
    test_feature(fa);
    test_arithmetic_operation_beuchat(fa);
    test_sqrt(fa);
    test_frob(fa);
    test_io(fa);

    // test for aranha's methods
    field_init(fb, "bn254_fp2b");
    test_feature(fb);
    test_arithmetic_operation_aranha(fb);
    test_sqrt(fb);
    test_frob(fb);
    test_io(fb);

    field_clear(fa);
    field_clear(fb);

    fprintf(stderr, "ok\n");

    return 0;
}
