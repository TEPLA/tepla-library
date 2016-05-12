#include <assert.h>

#include "rdtsc.h"

#include "../ec_bn254_lcl.h"

#define N 1000
#define M 100

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
    Element d1, d2, d3;

    mpz_t exp;

    //--------------------
    //  init
    //--------------------
    element_init(a, f);
    element_init(b, f);
    element_init(c, f);
    element_init(d, f);

    element_init(d1, f->base);
    element_init(d2, f->base);
    element_init(d3, f->base);

    //--------------------
    //  add
    //--------------------
    element_set_str(a, "17A767D1D0F35B9B2CE7CF00A9D036B7E087E24F1CBFFEF8C599F75DFDAD470B 107B87895D05703A57532E47A7BD8CA9C6406B5378C7E81F064749BB848490B2 44C128CE2774FC5C7022584909212BAB973CE5CA964D1B6A0ED9CCED1FA77DD 61C632B1BA796F688FF4A04475242B4116296B8A873B6C535C090278D88A7DB 685D0DA4D199864C925BAC4FD9EDBC9522367F43054E28E4E529A0FBF9BA604 7533624DC55FB89EDD5417E8F7E88E0B4A1F4CABD11353BEB3BD55FB62E7555");
    element_set_str(b, "19A7AA060E14E19A21FEF364C2BA8C3015EE7951DC51B928A55CE8C77E8998D9 43D2B55CCF3C20793C1402EE521F63E0DA933666152D14E6B35670281BB24E7 1A27F4DA83812A83B3189BAA7B9CB92F8401337F9334D9E155F66F31F33DE5DC 189A48DEB41106F112E6CCC46C55110E9BCF02C2B727B373CC6E46F3F80147F5 1A6243EC30C85C55E948EFF298055B8B5991E3CBD578CAB7DF4C1E52427329B7 2331BDFFF829F8D50BC7FD2BEF8607985B8C3F4805BA05F0F5CAB75DAD1E0194");
    element_set_str(d, "DDE16D341C72D770070A7CCE5A5C0A679741C60E111B81FEC7680257C36DFE3 14B8B2DF29F93241EB146E768CDF82E7D3E99EB9DA1AB96D717CB0BE063FB599 1E74076765F87A497A1AC12F0C2ECBEA3D7501DC3C99AB97F6E40C00C5385DB9 1EB6AC09CFB89DE79BE616C8B3A753C2AD31997B5F9B6A39022ED71B8589EFD0 20E814C67DE1F4BAB26EAAB795A43754ABB54BC005CDAD462D9EB862020ECFBB 713F920373EE4A0AB272411F81F8E37932BF4D2AACB3B2B62862CBD634C76E8");

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

    element_set_str(d1, "12B60569C8620CA5D70141D319878E5B060A4AA80413BEC52B5173D61E8A3387 1EF2082B9DCA3ABF201AA99B0230714369969A04C6064A7E9288B9551F95CA36");
    element_set_str(d2, "17B103CB3F5567B2C3A06E07CD31923CD87D72406E23B25B5066E61D029991F1 1791B2FB67A6FE4BC47F7719F2B5720ED863499E582E787A3EE2F513B4BDB4CF");
    element_set_str(d3, "2264428EB266F97E47163B2BE433146D0D5FBD32F6552D6079B4D85719F10635 14418C5689124B593F626719E1A07B7D0A2DA48E5BFEC13BDE16F0B8CC521F9B");

    element_set(((Element *)d->data)[0], d1);
    element_set(((Element *)d->data)[1], d2);
    element_set(((Element *)d->data)[2], d3);

    assert(element_cmp(c, d) == 0);

    t1 = rdtsc();
    for (i = 0; i < N; i++) {
        element_mul(c, a, b);
    }
    t2 = rdtsc();

    printf("element mul: %.2lf [clock]\n", (double)(t2 - t1) / N);

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
        element_inv(c, a);
    }
    t2 = rdtsc();

    printf("element inv: %.2lf [clock]\n", (double)(t2 - t1) / N);

    //--------------------
    //  pow
    //--------------------
    mpz_init(exp);

    mpz_set_str(exp, "AA4A2EE5234E0E95E8BB01F6B4A67F0EE8F2ADC1AA153C48D163AA85F3F534C", 16);

    element_pow(c, a, exp);
    element_set_str(d, "1CE3C3F32AC26429F4B8CDA845C051A4E88297DCCF33466FE03DDBCB06D7F83C 12AB1B1C4204A483CB230C3A25F3D3498CBC006D7DFF3F9D3565544BF09EEECD 16449675C385AA37E54DC18FABD35B84D1714E22991BEBD5CE7EFFBBEB110517 4B8B39A2E4C86D3BBD58B771A220F62E988136ECBB92D1213B3D268152A4E14 5C3CF6A5663F06A4F805674B29A6D2D7B9F0A8DA7E366EC239420E78B5BFAB6 C762C8ECC272929A992D55227F9C2927A04731122265C0E217B9B47AF524AA1");

    assert(element_cmp(c, d) == 0);

    mpz_set(exp, f->order);

    for (i = 0; i < 50; i++)
    {
        element_random(a);

        element_pow(b, a, exp);

        assert(element_cmp(b, a) == 0);
    }

    t1 = rdtsc();
    for (i = 0; i < M; i++) {
        element_pow(b, a, exp);
    }
    t2 = rdtsc();

    printf("element pow with order: %.2lf [clock]\n", (double)(t2 - t1) / M);

    mpz_clear(exp);

    //--------------------
    //  clear
    //--------------------
    element_clear(a);
    element_clear(b);
    element_clear(c);
    element_clear(d);

    element_clear(d1);
    element_clear(d2);
    element_clear(d3);
}

void test_arithmetic_operation_aranha(Field f)
{
    int i;
    unsigned long long int t1, t2;
    Element a, b, c, d;
    Element d1, d2, d3;

    mpz_t exp;

    //--------------------
    //  init
    //--------------------
    element_init(a, f);
    element_init(b, f);
    element_init(c, f);
    element_init(d, f);

    element_init(d1, f->base);
    element_init(d2, f->base);
    element_init(d3, f->base);

    //--------------------
    //  add
    //--------------------
    element_set_str(a, "17A767D1D0F35B9B2CE7CF00A9D036B7E087E24F1CBFFEF8C599F75DFDAD470B 107B87895D05703A57532E47A7BD8CA9C6406B5378C7E81F064749BB848490B2 44C128CE2774FC5C7022584909212BAB973CE5CA964D1B6A0ED9CCED1FA77DD 61C632B1BA796F688FF4A04475242B4116296B8A873B6C535C090278D88A7DB 685D0DA4D199864C925BAC4FD9EDBC9522367F43054E28E4E529A0FBF9BA604 7533624DC55FB89EDD5417E8F7E88E0B4A1F4CABD11353BEB3BD55FB62E7555");
    element_set_str(b, "19A7AA060E14E19A21FEF364C2BA8C3015EE7951DC51B928A55CE8C77E8998D9 43D2B55CCF3C20793C1402EE521F63E0DA933666152D14E6B35670281BB24E7 1A27F4DA83812A83B3189BAA7B9CB92F8401337F9334D9E155F66F31F33DE5DC 189A48DEB41106F112E6CCC46C55110E9BCF02C2B727B373CC6E46F3F80147F5 1A6243EC30C85C55E948EFF298055B8B5991E3CBD578CAB7DF4C1E52427329B7 2331BDFFF829F8D50BC7FD2BEF8607985B8C3F4805BA05F0F5CAB75DAD1E0194");
    element_set_str(d, "C2BAD559F083D3394B274E56C8AC2DF95555BA0F911B80DC3F6E0257C36DFD1 14B8B2DF29F93241EB146E768CDF82E7D3E99EB9DA1AB96D717CB0BE063FB599 1E74076765F87A497A1AC12F0C2ECBEA3D7501DC3C99AB97F6E40C00C5385DB9 1EB6AC09CFB89DE79BE616C8B3A753C2AD31997B5F9B6A39022ED71B8589EFD0 20E814C67DE1F4BAB26EAAB795A43754ABB54BC005CDAD462D9EB862020ECFBB 5618FA2947FF45D3F68F12A7F049070AF0D3412C2CB3B193A068CBD634C76D6");

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

    element_set_str(d1, "146BB4632F1B91494DAE710AC7E835FDFFF75A2707DFBE09295843756DA6346 2155662DCC9D7D44708AB85FD6E9FF81E547D66A8AB6AFF646223502A93E21AD");
    element_set_str(d2, "194E7CC6B62D623D1369E61295147C0A42444DD6F1A59140C700CBD2A18D4AC7 12E1AF5824A251E950F0AF165D847ACA220F73B390A82BFDE8668CB46B01113E");
    element_set_str(d3, "8CCC7470197C4E606A7B0F3711F447F5763BA2CE727905CE9D48082AA674FB2 FB557A5917142AEEE89F4BCE10D91100CCB6BE1AF42E2C3BA8949E172CBF875");

    element_set(((Element *)d->data)[0], d1);
    element_set(((Element *)d->data)[1], d2);
    element_set(((Element *)d->data)[2], d3);

    assert(element_cmp(c, d) == 0);

    t1 = rdtsc();
    for (i = 0; i < N; i++) {
        element_mul(c, a, b);
    }
    t2 = rdtsc();

    printf("element mul: %.2lf [clock]\n", (double)(t2 - t1) / N);

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
        element_inv(c, a);
    }
    t2 = rdtsc();

    printf("element inv: %.2lf [clock]\n", (double)(t2 - t1) / N);

    //--------------------
    //  pow
    //--------------------
    mpz_init(exp);

    mpz_set_str(exp, "AA4A2EE5234E0E95E8BB01F6B4A67F0EE8F2ADC1AA153C48D163AA85F3F534C", 16);

    element_pow(c, a, exp);
    element_set_str(d, "14055E48244EF8B9B68F5478CEFA7AE68A29B0CC5CE90210F8955E9984292DCF 3A60310B2D38A4E0326E5F3CC54CA1A8FE03C8F3F7A594595F4E39ACBCD2FA7 C885D7A1B81658977DB82656240B4525C645397B68AD947D2FAA4875A49512B 1EA6DAA006964358783B5B762E0BDC91A7E006E0A921A6BB18697F690826070B 1DE6A2C60400303328BC818FECF0B93929685CE72D1B6C1046B54587B6FB379C 8458030B512D15609F20153ED0043E8D649EE8DAF3628C8CF75C73A0882CF94");

    assert(element_cmp(c, d) == 0);

    mpz_set(exp, f->order);

    for (i = 0; i < 50; i++)
    {
        element_random(a);

        element_pow(b, a, exp);

        assert(element_cmp(b, a) == 0);
    }

    t1 = rdtsc();
    for (i = 0; i < M; i++) {
        element_pow(b, a, exp);
    }
    t2 = rdtsc();

    printf("element pow with order: %.2lf [clock]\n", (double)(t2 - t1) / M);

    mpz_clear(exp);

    //--------------------
    //  clear
    //--------------------
    element_clear(a);
    element_clear(b);
    element_clear(c);
    element_clear(d);

    element_clear(d1);
    element_clear(d2);
    element_clear(d3);
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

    for (i = 0; i < 50; i++)
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
    for (i = 0; i < M; i++) {
        element_sqrt(c, b);
    }
    t2 = rdtsc();
    printf("element sqrt: %.2lf [clock]\n", (double)(t2 - t1) / M);

    element_clear(a);
    element_clear(b);
    element_clear(c);
    element_clear(d);
}

//============================================
//   i/o test
//============================================
void test_io(Field f)
{
    int i;
    unsigned long long int t1, t2;
    char a_str[390];

    size_t blen;
    unsigned char b_str[192];

    Element a, b, c;

    element_init(a, f);
    element_init(b, f);
    element_init(c, f);

    for (i = 0; i < 100; i++)
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

    for (i = 0; i < 100; i++)
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
    field_init(fa, "bn254_fp6a");
    test_feature(fa);
    test_arithmetic_operation_beuchat(fa);
    test_sqrt(fa);
    test_io(fa);

    // test for aranha's methods
    field_init(fb, "bn254_fp6b");
    test_feature(fb);
    test_arithmetic_operation_aranha(fb);
    test_sqrt(fb);
    test_io(fb);

    field_clear(fa);
    field_clear(fb);

    fprintf(stderr, "ok\n");

    return 0;
}
