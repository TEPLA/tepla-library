#include <assert.h>

#include "rdtsc.h"

#include "../ec_bn254_lcl.h"

#define N 1000
#define M 100

#define t 80            // security Level for SHA
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
    gmp_fprintf(stdout, "   field order: %Zx^%d\n", *field_get_char(ec->field), field_get_degree(ec->field));
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
void test_arithmetic_operation_aranha(const EC_GROUP ec)
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

    assert(point_is_on_curve(a));

    t1 = rdtsc();
    for (i = 0; i < M; i++) {
        point_random(a);
    }
    t2 = rdtsc();

    printf("point random: %.2lf [clock]\n", (double)(t2 - t1) / M);

    //-------------------
    //  add/dob
    //-------------------
    point_add(b, b, a);
    point_add(c, b, c);

    assert(point_cmp(b, a) == 0);
    assert(point_cmp(c, b) == 0);

    point_set_infinity(d);

    point_dob(d, d);

    assert(point_is_infinity(d));

    point_set_str(a, "[CEA56EDB3F5855752214F0FB38227F5F05FEF183D568E3438FFA0C59662D934 1E3B3424F54E3041A2695E2DA8B74008C633E30964CC7E439667FFBA9AD4CA04,10B227CCDA6FE33B195E8486AFD1C260EF9FFB9560F40EA61795A9EC387A2894 BE95E1BD28C3A1664D253012E6342EBE66A41D3E56F38D0C8E2EBDECCE0F2B]");
    point_set_str(b, "[139EF51833157A9BD7D08CBF0733A2A7EFA1E332AF1F89A25E8FBF26366BA0A 10C7EF12DD56C5862AD7D23CD46003B298B44762C83D97C34DA4D06123D15E4A,14DCC68B6147E8195A42AE1DB063CEC65EC79FC5E7A721167725239567160A2C DC6FB2556E0A7894F356C21408AFF3E1F2F9FA4DF69250F8D2FB7A44C02F219]");
    point_set_str(d, "[AA3D89185E6C19038841BB9C57B08877B982CC4692C68471092F09685F06AA5 193FFB7AE74A2EA693F57D50C312E9D8260434B2030D060BC3CAA04B167DC34F,120BB54CF55F14B0B9AD94EEB4132ACC56BFAC1FC5F6632326AFFCA558127D89 18FBC41432A0B4DB601194348389D3C10CEAB66D3539124D256915F3C03F2856]");

    point_add(c, a, b);

    assert(point_cmp(c, d) == 0);

    t1 = rdtsc();
    for (i = 0; i < N; i++) {
        point_add(c, a, b);
    }
    t2 = rdtsc();

    printf("point add: %.2lf [clock]\n", (double)(t2 - t1) / N);

    element_init(dx, ec->field);
    element_init(dy, ec->field);

    element_set_str(dx, "1A78312372D449F8063336935B3162C434414C2E64A9ABABF12FF8047461941B 1ACCC23BE9B1AB7B97EB7AF4B0AD5719CF1D9DDCED78C55B9091651D70AD7452");
    element_set_str(dy, "3AB0B0AE7C78F96AED3E1A18DE1D8C907A394F7E2EBB72FFE2CD39477EA6BCA 1E62CFD510421FFB1AF4906CEFC410FD59E3D447FB697681E530FD14EC36593D");

    point_set_xy(d, dx, dy);

    point_dob(c, a);

    assert(point_cmp(c, d) == 0);

    t1 = rdtsc();
    for (i = 0; i < N; i++) {
        point_dob(c, a);
    }
    t2 = rdtsc();

    printf("point dob: %.2lf [clock]\n", (double)(t2 - t1) / N);

    //------------------
    // Frobenius
    //------------------
    element_set_str(dx, "2D7E6A04FBCF199172BD7093CDE7F4BAD92DFE2FE6B17AD083E9B4DAAA01DC0 1DD91398F9311CDEBFE6C0EE89041D2CCAB24DB9935C4165B1AF6C91CD321329");
    element_set_str(dy, "23A64D31121391EE97150DC1751ACE489378115B43C4980CBB4B7292983BB37F D3A9D9EDAA245ABF148E07D581F731E8B072D791B4AF6D1F239C41C8C7A92C3");
    point_set_xy(d, dx, dy);

    ec_bn254_tw_frob(c, a);
    assert(point_is_on_curve(c));
    assert(point_cmp(c, d) == 0);

    element_set_str(dx, "1F835E7DD2D95DCD62608B01C379BAEE072FC32E2F4D308D6350F2A89C6B13C9 41049BCFAF4DE27009F18491A6A40B3ED5A3D139CC86A23085964F7BA8B184F");
    element_set_str(dy, "14713CB565901CC6A0D5C8F9502E3DA77181046A9F0BF16D8F6A5613C785D77F 2464CEA082D73C6053E7284FED19CBD9A2BA5BE2C1A90C869A71D1421331F0E8");
    point_set_xy(d, dx, dy);

    ec_bn254_tw_frob2(c, a);
    assert(point_is_on_curve(c));
    assert(point_cmp(c, d) == 0);

    element_set_str(dx, "1E3B3424F54E3041A2695E2DA8B74008C633E30964CC7E439667FFBA9AD4CA04 CEA56EDB3F5855752214F0FB38227F5F05FEF183D568E3438FFA0C59662D934");
    element_set_str(dy, "17D17512DEC6E13231F3FBE8AE531BFCDA8EEA4BC3B6806EBB48D6D67C44C94 17E8C6E3655DBA55C8EB6D02A7E08CE9D619D286E4B50941B4C63BE373856D50");
    point_set_xy(d, dx, dy);

    ec_bn254_tw_frob3(c, a);
    assert(point_is_on_curve(c));
    assert(point_cmp(c, d) == 0);

    element_clear(dx);
    element_clear(dy);

    //------------------
    // neg/sub
    //------------------
    point_neg(c, a);
    point_add(d, a, c);
    point_sub(b, a, a);

    assert(point_is_infinity(d));
    assert(point_is_infinity(b));

    //-------------------
    //  mul
    //-------------------

    mpz_init(scalar);

    mpz_set(scalar, ec->order);

    for (i = 0; i < 100; i++)
    {
        point_random(a);

        ec_bn254_fp2_mul(b, scalar, a);	// TEPLA's scalar multiplication

        ec_bn254_fp2_mul_end(c, scalar, a);

        assert(point_is_infinity(b));
        assert(point_cmp(c, b) == 0);
    }

    //-------------------
    // point is on curve
    //-------------------

    //定義式は満たすが,G2に属さない点c
    point_set_str(c, "[1307ea577679b2e61d84508198bb1de6a498cf95c963f1da5fc5e7f99e6844ab ed34cf0a3f89214af672c6bbc5c362c5f48693015fe1a4ebe58ad8adfa66aea,29ea17af844d0d9444f92962d84a5ecacfd1aa3df4e6a54ec53a821c9fcbbb ac581e400ff84f829143eccfb72db63f61b47bade96812f15a079c87bed92cd]");
    assert(point_is_on_curve(a));
    assert(point_is_on_curve(c) == 0);

    t1 = rdtsc();
    for (i = 0; i < M; i++) {
        ec_bn254_fp2_mul_end(b, scalar, a);
    }
    t2 = rdtsc();

    printf("point mul with endomorphism: %.2lf [clock]\n", (double)(t2 - t1) / M);

    t1 = rdtsc();
    for (i = 0; i < M; i++) {
        ec_bn254_fp2_mul(b, scalar, a);
    }
    t2 = rdtsc();

    printf("point mul with binary method: %.2lf [clock]\n", (double)(t2 - t1) / M);

    //-------------------
    //  clear
    //-------------------
    point_clear(a);
    point_clear(b);
    point_clear(c);
    point_clear(d);
}

void test_arithmetic_operation_beuchat(const EC_GROUP ec)
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

    assert(point_is_on_curve(a));

    t1 = rdtsc();
    for (i = 0; i < M; i++) {
        point_random(a);
    }
    t2 = rdtsc();

    printf("point random: %.2lf [clock]\n", (double)(t2 - t1) / M);

    //-------------------
    //  add/dob
    //-------------------
    point_add(b, b, a);
    point_add(c, b, c);

    assert(point_cmp(b, a) == 0);
    assert(point_cmp(c, b) == 0);

    point_set_infinity(d);

    point_dob(d, d);

    assert(point_is_infinity(d));

    point_set_str(a, "[6D2E4115FA177379A504A0EE4EF53767DE51C6364AAB69D4064529EC1FD047A 635B2C858AA4F4A3DB8AA17A588B037CAFFD36678F76E3F3369DFC90C6878C7,193E877C82EFCA81EC2815906630B837BBC6976CC8A7958E6A40D1B190FF2E5F E8A77E88AFCEE9F806DC15BF50EADD138320F1A5A87E78DDE86FA7A867300D]");
    point_set_str(b, "[1A8F5DAB09EE4290F95FE4C824C153E355D55B6CF94B998C6203FEC3D81377CF 15A19F2704C4BDBAAE39A5E26772A3E4E7EC7A9E205651F8822298766DE044FF,1C566EB3917F06B05E0A786BD8030CAFCCDB62864DD0E2A22A9B6817B310FD53 6A0927BB33EB263F45CAB921A20E67A1BD8A791D6EB0415AC92C9B1F74D16D1]");
    point_set_str(d, "[143D414F99AA18C844B331064C9DD66363EBA3D852250CBCF8C9D4B33E0C4C1C 225865D85EC7A34647CA55E026BD1FA201E0C4E8C66F7A43E69AF708F410A0FF,FACA1388C034CF614A72E06EE60DEDC4880CDBD368E5BEC2795130B266FFB9E 1681217E50705AB9A21FEB62E0BF9A5657EB27C3AED3323FE9C57058358735A9]");

    point_add(c, a, b);

    assert(point_cmp(c, d) == 0);

    t1 = rdtsc();
    for (i = 0; i < N; i++) {
        point_add(c, a, b);
    }
    t2 = rdtsc();

    printf("point add: %.2lf [clock]\n", (double)(t2 - t1) / N);

    element_init(dx, ec->field);
    element_init(dy, ec->field);

    element_set_str(dx, "33F550F9A63EF53C786BF7BDFDAB1538CD76A3FCED3C9DBC3307DD4F354775A C814AE99C91C71845F0B51E4349520908E48C70181313D70C05F6ED24EC1F36");
    element_set_str(dy, "3766ED0DD7C988DB76770081A298DAA924D0E3279726F9B5504129AFA3E57B9 520CE8A563F88AF882AB99086BFDBBDCEF9DE65879AB234DFF5AAFD5BEE7E4F");

    point_set_xy(d, dx, dy);

    point_dob(c, a);

    assert(point_cmp(c, d) == 0);

    t1 = rdtsc();
    for (i = 0; i < N; i++) {
        point_dob(c, a);
    }
    t2 = rdtsc();

    printf("point dob: %.2lf [clock]\n", (double)(t2 - t1) / N);

    element_clear(dx);
    element_clear(dy);

    //------------------
    // neg/sub
    //------------------
    point_neg(c, a);
    point_add(d, a, c);
    point_sub(b, a, a);

    assert(point_is_infinity(d));
    assert(point_is_infinity(b));

    //-------------------
    //  mul
    //-------------------
    mpz_init(scalar);

    mpz_set(scalar, ec->order);

    for (i = 0; i < 100; i++)
    {
        point_random(a);

        point_mul(b, scalar, a);

        ec_bn254_fp2_mul_end(c, scalar, a);

        assert(point_is_infinity(b));
        assert(point_cmp(c, b) == 0);
    }

    t1 = rdtsc();
    for (i = 0; i < M; i++) {
        point_mul(b, scalar, a);
    }
    t2 = rdtsc();

    printf("point mul with endomorphism: %.2lf [clock]\n", (double)(t2 - t1) / M);

    t1 = rdtsc();
    for (i = 0; i < M; i++) {
        ec_bn254_fp2_mul(b, scalar, a);
    }
    t2 = rdtsc();

    printf("point mul with binary method: %.2lf [clock]\n", (double)(t2 - t1) / M);

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

    point_mul(Q, ec->order, P);

    assert(point_is_infinity(Q));

    point_map_to_point(Q, MAP_STR, sizeof(MAP_STR), t);

    assert(point_cmp(Q, P) == 0);

    t1 = rdtsc();
    for (i = 0; i < M; i++) {
        point_map_to_point(P, MAP_STR, sizeof(MAP_STR), t);
    }
    t2 = rdtsc();

    printf("point map to point in 128 security: %.2lf [clock]\n", (double)(t2 - t1) / M);

    point_clear(P);
    point_clear(Q);
}

//============================================
//  i/o test
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

    assert(point_is_infinity(Q));

    for (i = 0; i < 100; i++)
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

    for (i = 0; i < 100; i++)
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
    curve_init(eca, "ec_bn254_twa");
    test_feature(eca);
    test_arithmetic_operation_beuchat(eca);
    test_map_to_point(eca);
    test_io(eca);

    // test for aranha's methods
    curve_init(ecb, "ec_bn254_twb");
    test_feature(ecb);
    test_arithmetic_operation_aranha(ecb);
    test_map_to_point(ecb);
    test_io(ecb);

    curve_clear(eca);
    curve_clear(ecb);

    fprintf(stderr, "ok\n");

    return 0;
}
