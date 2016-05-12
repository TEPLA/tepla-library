#include <assert.h>

#include "rdtsc.h"

#include "../ec_bn254_lcl.h"

#define N 1000
#define M 100

#define t 80            // security level for SHA
#define MAP_STR "LCIS"  // string for map to point

//============================================
//  Feature of Elliptic Curve
//============================================
void test_feature(const EC_GROUP ec)
{
    fprintf(stdout, "---\n");
    fprintf(stdout, "Elliptic Curve Type: %s\n", curve_get_name(ec));
    fprintf(stdout, "   Y^2 = X^3 + aX + b\n");
    fprintf(stdout, "   a: ");
    element_print(ec->a);
    fprintf(stdout, "   b: ");
    element_print(ec->b);
    gmp_fprintf(stdout, "   field order: %Zx\n", ec->field->order);
    fprintf(stdout, "   generator of curve: ");
    point_print(ec->generator);
    gmp_fprintf(stdout, "   order: %Zx\n", ec->order);
    gmp_fprintf(stdout, "   trace: %Zx\n", ec->trace);
    gmp_fprintf(stdout, "   cofactor: %Zx\n", ec->cofactor);
    fprintf(stdout, "---\n");
}

//============================================
//  test for arithmetic operations of EC
//============================================
void test_arithmetic_operation_beuchat(const EC_GROUP ec)
{
    int i, j;
    unsigned long long int t1, t2;

    mpz_t scalar;

    EC_POINT x, y, z, w;

    //-------------------
    //  init
    //-------------------
    point_init(x, ec);
    point_init(y, ec);
    point_init(z, ec);
    point_init(w, ec);

    //-------------------
    //  random
    //-------------------
    for (i = 0; i < 100; i++)
    {
        point_random(x);
        point_random(y);

        assert(point_is_on_curve(x));
        assert(point_is_on_curve(y));
    }

    t1 = rdtsc();
    for (i = 0; i < M; i++) {
        point_random(x);
    }
    t2 = rdtsc();

    printf("point random: %.2lf [clock]\n", (double)(t2 - t1) / M);

    //-------------------
    //  add/dob
    //-------------------
    point_add(w, w, x);
    point_add(z, y, z);

    assert(point_cmp(w, x) == 0);
    assert(point_cmp(z, y) == 0);

    point_set_infinity(w);
    point_dob(z, w);

    assert(point_is_infinity(z));

    point_set_str(x, "[F57944869FBC67AAF4EE00F042A02A152340E5900385E9914C28C52855DEAB,1A48AAA8060DB83BE9D6504E3018BA47F13F9C3CD789EF776A1F59989744E488]");
    point_set_str(y, "[10B12DA9453E305F931B61DAA1F8AEC538242F87864097D568C657F5D689B633,11D8AAB4E05947B76A60CD5779F58DCA603160E9A1A80572DC14FADBE1E47A5F]");
    point_set_str(w, "[BE00F9C60DFFB7B2C60D60CC5972376DF0DB28C30BD532B1CE028C8F23F3BBD,1276586AB64BD85F11757F266CC049D99C30CA6EDAA5D2029F85F064FDCBAAA2]");

    point_add(z, x, y);

    assert(point_cmp(z, w) == 0);

    t1 = rdtsc();
    for (i = 0; i < N; i++) {
        point_add(z, x, y);
    }
    t2 = rdtsc();

    printf("point add: %.2lf [clock]\n", (double)(t2 - t1) / N);

    point_set_str(w, "[17769C18B8240D6A0314ABD9D0D1826805D76C5766B5CB19017DC4B6A7AE6A30,1F9EC14C82F15EEE45D36CCCFB5B8B511D371F6724CAE85CDF19A523841597D4]");

    point_dob(z, x);

    assert(point_cmp(z, w) == 0);

    t1 = rdtsc();
    for (i = 0; i < N; i++) {
        point_dob(z, x);
    }
    t2 = rdtsc();

    printf("point dob: %.2lf [clock]\n", (double)(t2 - t1) / N);

    //------------------
    // neg/sub
    //------------------
    point_neg(z, x);
    point_add(w, x, z);
    point_sub(z, x, x);

    assert(point_is_infinity(w));
    assert(point_is_infinity(z));

    //-------------------
    //  mul
    //-------------------
    mpz_init(scalar);

    mpz_init_set_str(scalar, "2370FB049D410FBDC023FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF", 16);

    ec_bn254_fp_mul(y, scalar, x);

    ec_bn254_fp_point_endomorphism(z, x);

    assert(point_cmp(y, z) == 0);

    point_set_str(x, "[8957DD9D913386BDEE43720A14C1AF42B44C072F78516FBEB12D80DA16A52C9,1520FC7A1090AFE300849C46AA12B287F599C54EC411D498D24E1968B20848A0]");
    point_set_str(y, "[70E32B4C02DB585DEAE420B779B345AE6BDF34A614258B4EF67A7E98B984D14,C53D7896CDE706964C55D1A9F31A74E1ED4C490E16496B8D77C95B4E0A8ED73]");

    mpz_set_str(scalar, "1C26BA60159C76F9C565E8716605BFBF2881FF4D82B2B4081F6EB20604CB1830", 16);

    point_mul(z, scalar, x);

    assert(point_cmp(z, y) == 0);

    t1 = rdtsc();
    for (i = 0; i < M; i++) {
        point_mul(z, scalar, x);
    }
    t2 = rdtsc();

    printf("point mul endomorphism: %.2lf [clock]\n", (double)(t2 - t1) / M);

    for (i = 1; i < 100; i++)
    {
        point_random(x);
        point_set_infinity(w);
        mpz_set_ui(scalar, i);

        for (j = 0; j < i; j++) {
            point_add(w, w, x);
        }

        point_mul(z, scalar, x);

        assert(point_cmp(z, w) == 0);
    }

    mpz_set(scalar, ec->order);

    for (i = 0; i < 100; i++)
    {
        point_random(x);

        ec_bn254_fp_mul(z, scalar, x);

        assert(point_is_infinity(z));
    }

    t1 = rdtsc();
    for (i = 0; i < M; i++) {
        ec_bn254_fp_mul(z, scalar, x);
    }
    t2 = rdtsc();

    printf("point mul normal: %.2lf [clock]\n", (double)(t2 - t1) / M);

    mpz_clear(scalar);

    //-------------------
    //  clear
    //-------------------
    point_clear(x);
    point_clear(y);
    point_clear(z);
    point_clear(w);
}

void test_arithmetic_operation_aranha(const EC_GROUP ec)
{
    int i, j;
    unsigned long long int t1, t2;

    mpz_t scalar;

    EC_POINT x, y, z, w;

    //-------------------
    //  init
    //-------------------
    point_init(x, ec);
    point_init(y, ec);
    point_init(z, ec);
    point_init(w, ec);

    //-------------------
    //  random
    //-------------------
    for (i = 0; i < 100; i++)
    {
        point_random(x);
        point_random(y);

        assert(point_is_on_curve(x));
        assert(point_is_on_curve(y));
    }

    t1 = rdtsc();
    for (i = 0; i < M; i++) {
        point_random(x);
    }
    t2 = rdtsc();

    printf("point random: %.2lf [clock]\n", (double)(t2 - t1) / M);

    //-------------------
    //  add/dob
    //-------------------
    point_add(w, w, x);
    point_add(z, y, z);

    assert(point_cmp(w, x) == 0);
    assert(point_cmp(z, y) == 0);

    point_set_infinity(w);
    point_dob(z, w);

    assert(point_is_infinity(z));
    point_set_str(x, "[CA50EE766C535743EC45F534CBC42A703B0A4C4A6DEB82C1CF41CE76FAFACA7,7E1B2DABD6BE6B6C20438670E2A1D11A0896ADCAFFD539782848DA0C49266ED]");
    point_set_str(y, "[51F6E2F07DB4239C6E9CE8D0CEE6BC46E8F125867AD40A433AC6D3DC66C87DB,18A9843C4000A034C21C2C8C76E9C5F6801BF5DF0F4AC235199984B06E1F83B7]");
    point_set_str(w, "[10C813944B3FD19F0A86EB304AB7C93AC61BA1DC625F295C148744C7664E7DA3,1F6F0F4C9A2C169D98A397322CF54CE812DD096FA914C8DCF03CAD3F6628F43E]");
    point_add(z, x, y);
    assert(point_cmp(z, w) == 0);

    t1 = rdtsc();
    for (i = 0; i < N; i++) {
        point_add(z, x, y);
    }
    t2 = rdtsc();

    printf("point add: %.2lf [clock]\n", (double)(t2 - t1) / N);

    point_set_str(w, "[B3C783B1497247BC18C19C8238928AD36D30F79BF4D6228B214E3B796D3078D,1C15B880392076285F680BB6053CD7D7338B0B246E93A58C78E9B4032A468A28]");

    point_dob(z, x);

    assert(point_cmp(z, w) == 0);

    t1 = rdtsc();
    for (i = 0; i < N; i++) {
        point_dob(z, x);
    }
    t2 = rdtsc();

    printf("point dob: %.2lf [clock]\n", (double)(t2 - t1) / N);

    //------------------
    // neg/sub
    //------------------
    point_neg(z, x);
    point_add(w, x, z);
    point_sub(z, x, x);

    assert(point_is_infinity(w));
    assert(point_is_infinity(z));

    //-------------------
    //  mul
    //-------------------
    mpz_init_set_str(scalar, "252364824000000126CD8900000000024908FFFFFFFFFFFCF9FFFFFFFFFFFFF6", 16);
    ec_bn254_fp_mul(y, scalar, x);
    ec_bn254_fp_point_endomorphism(z, x);

    assert(point_cmp(y, z) == 0);

    point_set_str(x, "[237E5B966EC126BA4117DDAB68326EB05299EBAB87E2C5BD505F13073344C771,94E915680ECEC5BC3BC0990724FA981C0119E2D8E881BFE40FEAAFE4CCE8F6A]");
    point_set_str(y, "[13C9716E0ECB7466584B3F91802047C9E90EAD6963131D2E5CA638908134301A,242929DEFDD844A94DF3F487FCBEC3B95B71D472F012D3E7F0C5991F6D5500A7]");

    mpz_set_str(scalar, "1C26BA60159C76F9C565E8716605BFBF2881FF4D82B2B4081F6EB20604CB1830", 16);

    point_mul(z, scalar, x);

    assert(point_cmp(z, y) == 0);

    t1 = rdtsc();
    for (i = 0; i < M; i++) {
        point_mul(z, scalar, x);
    }
    t2 = rdtsc();

    printf("point mul endomorphism: %.2lf [clock]\n", (double)(t2 - t1) / M);

    for (i = 1; i < 100; i++)
    {
        point_random(x);
        point_set_infinity(w);
        mpz_set_ui(scalar, i);

        for (j = 0; j < i; j++) {
            point_add(w, w, x);
        }

        point_mul(z, scalar, x);

        assert(point_cmp(z, w) == 0);
    }

    mpz_set(scalar, ec->order);

    for (i = 0; i < 100; i++)
    {
        point_random(x);

        ec_bn254_fp_mul(z, scalar, x);

        assert(point_is_infinity(z));
    }

    t1 = rdtsc();
    for (i = 0; i < M; i++) {
        ec_bn254_fp_mul(z, scalar, x);
    }
    t2 = rdtsc();

    printf("point mul normal: %.2lf [clock]\n", (double)(t2 - t1) / M);

    //-------------------
    //  clear
    //-------------------
    point_clear(x);
    point_clear(y);
    point_clear(z);
    point_clear(w);
}
//============================================
//  MAP to POINT test
//============================================
void test_map_to_point(const EC_GROUP ec)
{
    int i;
    unsigned long long int t1, t2;
    EC_POINT P, Q;

    point_init(P, ec);
    point_init(Q, ec);

    point_map_to_point(P, MAP_STR, sizeof(MAP_STR), t);

    assert(point_is_on_curve(P));

    point_map_to_point(Q, MAP_STR, sizeof(MAP_STR), t);

    assert(point_cmp(Q, P) == 0);

    t1 = rdtsc();
    for (i = 0; i < N; i++) {
        point_map_to_point(P, MAP_STR, sizeof(MAP_STR), t);
    }
    t2 = rdtsc();

    printf("point map to point in 128 security: %.2lf [clock]\n", (double)(t2 - t1) / N);

    point_clear(P);
    point_clear(Q);
}

//============================================
//  i/o test of EC
//============================================
void test_io(const EC_GROUP ec)
{
    int i;
    unsigned long long int t1, t2;

    size_t osize;
    unsigned char os[128];

    char str[256];

    EC_POINT P, Q, R;

    point_init(P, ec);
    point_init(Q, ec);
    point_init(R, ec);

    //---------------------
    //  octet string
    //---------------------
    point_set_infinity(R);

    point_to_oct(os, &osize, R);
    point_from_oct(Q, os, osize);

    assert(point_is_infinity(Q));

    for (i = 0; i < 1000; i++)
    {
        point_random(P);

        point_to_oct(os, &osize, P);
        point_from_oct(Q, os, osize);

        assert(point_cmp(P, Q) == 0);
    }

    t1 = rdtsc();
    for (i = 0; i < N; i++) {
        point_to_oct(os, &osize, P);
    }
    t2 = rdtsc();

    printf("point to octet string: %.2lf [clock]\n", (double)(t2 - t1) / N);

    t1 = rdtsc();
    for (i = 0; i < N; i++) {
        point_from_oct(Q, os, osize);
    }
    t2 = rdtsc();

    printf("point from octet string: %.2lf [clock]\n", (double)(t2 - t1) / N);

    //---------------------
    //  string
    //---------------------
    point_set_infinity(R);

    point_get_str(str, R);
    point_set_str(Q, str);

    assert(point_is_infinity(Q));

    for (i = 0; i < 1000; i++)
    {
        point_get_str(str, P);
        point_set_str(Q, str);

        assert(point_cmp(P, Q) == 0);
    }

    t1 = rdtsc();
    for (i = 0; i < N; i++) {
        point_get_str(str, P);
    }
    t2 = rdtsc();

    printf("point get string: %.2lf [clock]\n", (double)(t2 - t1) / N);

    t1 = rdtsc();
    for (i = 0; i < N; i++) {
        point_set_str(Q, str);
    }
    t2 = rdtsc();

    printf("point set string: %.2lf [clock]\n", (double)(t2 - t1) / N);

    point_clear(P);
    point_clear(Q);
    point_clear(R);
}

//============================================
//  main program
//============================================
int main(void)
{
    EC_GROUP eca, ecb;

    // test for beuchat's methods
    curve_init(eca, "ec_bn254_fpa");
    test_feature(eca);
    test_arithmetic_operation_beuchat(eca);
    test_map_to_point(eca);
    test_io(eca);

    // test for aranha's methods
    curve_init(ecb, "ec_bn254_fpb");
    test_feature(ecb);
    test_arithmetic_operation_aranha(ecb);
    test_map_to_point(ecb);
    test_io(ecb);

    curve_clear(eca);
    curve_clear(ecb);

    fprintf(stderr, "ok\n");

    return 0;
}
