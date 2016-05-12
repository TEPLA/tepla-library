#include <assert.h>

#include "../ec_bn254_lcl.h"

#include "rdtsc.h"

#define N 100

//============================================
//  Feature of Finite Field
//============================================
void print_field_feature(const Field f)
{
    fprintf(stdout, "---\n");
    fprintf(stdout, "Field Type: %s\n", field_get_name(f));
    gmp_fprintf(stdout, "   order of field: %Zx^%d\n", *field_get_char(f), field_get_degree(f));
    fprintf(stdout, "---\n");
}

//============================================
//  Feature of Elliptic Curve
//============================================
void print_ec_feature(const EC_GROUP ec)
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
//  Feature of Pairing
//============================================
void test_feature(const EC_PAIRING p)
{
    fprintf(stdout, "---\n");
    fprintf(stdout, "Pairing Type: %s\n", pairing_get_name(p));
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
    EC_PAIRING p1;
    EC_POINT P, P1, P2, Q, Q1, Q2, R, S;
    Element d, e, f, g, h;

    gmp_randstate_t state;
    mpz_t a, b, c, order;

    //-------------------
    //  init
    //-------------------
    point_init(P, p->g1);
    point_init(P1, p->g1);
    point_init(P2, p->g1);
    point_init(R, p->g1);

    point_init(Q, p->g2);
    point_init(Q1, p->g2);
    point_init(Q2, p->g2);
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


    t1 = clock();
    for (i = 0; i < N; i++) {
        pairing_init(p1, p->pairing_name);
    }
    t2 = clock();

    printf("pairing init type (%s): %.5lf [msec]\n", pairing_get_name(p), (double)(t2 - t1) / CLOCKS_PER_SEC / N * 1000);

    for (i = 0; i < 10; i++)
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

        assert(element_is_one(d));
        pairing_map(f, R, Q, p);  // f = e( Q, aP)
        pairing_map(h, P, S, p);  // h = e(bQ,  P)
        pairing_map(e, R, S, p);  // e = e(bQ, aP)

        element_pow(d, g, a);	 // d = e(Q, P)^a

        assert(element_cmp(d, f) == 0);

        element_pow(d, g, b);    // d = e(Q, P)^b

        assert(element_cmp(d, h) == 0);

        element_pow(d, g, c);    // d = e(Q, P)^c

        assert(element_cmp(d, e) == 0);
    }

    for (i = 0; i < 10; i++)
    {
        point_random(P1);
        point_random(P2);
        point_random(Q1);
        point_random(Q2);

        pairing_map(d, P1, Q1, p);
        pairing_map(e, P2, Q2, p);
        element_mul(d, d, e);

        pairing_double_map(e, P1, Q1, P2, Q2, p);

        assert(element_cmp(d, e) == 0);
    }

    // test NTT draft
    if (strcmp(p->pairing_name, "ECBN254a") == 0)
    {
        point_set_str(P, "[A971735A70FBDD0F94D7D6EFBBC81BEA78D2D92A8510F3344038A416419AD97,9456E41754237447752A448282C0873785F724447E1299826F53AC556936D3F]");
        point_set_str(Q, "[115231D7B49901BA97CB93B5227F7F7F438A346532893DD5FAFD518950924AA9 DF12398FB78695A50BB3499B7E23B0D9035989B91A76D13AF7BC64374BFB8A6,51D0E087527BC9F41379FB0272EC91E5F28EE011B183EF7D6712EF3FC9A1A66 107E6654DC6C36E163B7867AECB98E4046084734524DBB562E73E5A811F678A]");
        element_set_str(h, "6A4E0DD1F7FD2F9E5DACAB02CEC9CE8254925C5DC6697E153F05A242CBCA8A8 751037182B5F93BCAB31B115A2C0A0DCC09C6DB7602E0551DD44925F3D364B3 13BE65D47487BF6D96C146C18855C1F87BF994F9F1048524568EA0CB9DC402AD 15F9E3D10B580FF1AB2282EF1DC39A88E06F93A18303E9520D99B86D665F5380 1611153BF02F1CF7985B98C3F3CB641D39283DBA55E22D1C614568F84959C6FC 166BD873D0C65DE66300A168BBDC16F0AB1B57A0809973239F2109A7D25AD349 22A0E22C097AEC1187087B7632C9B963B0E779BC8D09848C44D3EA95CD1C1F8C 4B6BFFB9EB68AD6A99ACF52B8AAD1D17D328847C6313201A6B659C9DAA5CDFE 1202BE31EB2BDCBEF9F3CC00F1B2CC35FADBE1A0D66CCBF40B024ADFA84C77D1 A1C6D26A6D683031D95C4369DB90F5FEE36D5008AA498D2CB6F2DDE6258CDA6 10BEF55B7539743CBEAB13E49116A143302F6F28CCD71A69860CEF5208483809 14D4B5014F840144D03C0C6B6010BB246EE6A69BF704D7542FBAA8F2D2A27308");
        pairing_map(g, P, Q, p);

        assert(element_cmp(g, h) == 0);
    }

    if (strcmp(p->pairing_name, "ECBN254b") == 0)
    {
        point_set_str(P, "[2074A81D4402A0B63B947335C14B2FC3C28FEA2973860F686114BEC4670E4EB7,6A41108087B20038771FC89FB94A82B2006034A6E8D871B3BC284846631CBEB]");
        point_set_str(Q, "[49EEDB108B71A87BFCFC9B65EB5CF1C2F89554E02DF4F8354E4A00F52183C77 1FB93AB676140E87D97226185BA05BF5EC088A9CC76D966697CFB8FA9AA8845D,CD04A1ED14AD3CDF6A1FE4453DA2BB9E686A637FB3FF8E2573644CC1EDF208A 11FF7795CF59D1A1A7D6EE3C3C2DFC765DEF1CAA9F14EA264E71BD7630A43C14]");
        element_set_str(h, "3E1F2693AC6D549898C78897EB158490A4832E296F888D30140500DB7BD3D12 A5A5405542F67384D683A48C281F3676B67554ED5DA1700784169A0B47A57E4 142715D6482BC6FA77377C9CBC2A51C047C16DE88483D5A889C7EF4DF5F03BDB 22371AF975DAE562F686988CDBBD02702C959BBF843A1FB3C7532D07BE3D7A3A 5D259DA3F3AAAA54A6AE5FE8272A5B79D7F4E5BDF3B5E3C815AD781113F7548 13CA93E1377EF0F6DD38FC2F96DBD3E8B0922F60D1F274EAC63DC1AF2EE9754C 1EBC54A76E844EB5D352945226FB103DE9EC1A4FC689B87FAA66EF8ABA79D3ED 48B66DAFCAEE86DB4D46AB71A9FE848443EF81F488D8366A727B39698CF7201 11EE0C12164133041C3DCF312CE111C845B60092818F7B72805D4AFF61427934 4052CA960900684A1B26C434B2776AA70736841474C16208CCD1A7C27927E19 843C37BC5BDBF253E3BCE568F5905A63867D8836855B74CBA0C800D5DC41B71 D467F3DA4FB329A5CB406D0A7B743A3A2FFCD09BF95EE8A856B94AF191D96AF");
        pairing_map(g, P, Q, p);

        assert(element_cmp(g, h) == 0);
    }


    t1 = clock();
    for (i = 0; i < N; i++) {
        p->pairing(g, Q, P, p);
    }
    t2 = clock();
    printf("optimal ate pairing: %.5lf [msec]\n", (double)(t2 - t1) / (CLOCKS_PER_SEC) / N * 1000);


    if (strcmp(p->pairing_name, "ECBN254a") == 0)
    {
        t1 = clock();
        for (i = 0; i < N; i++) {
            ec_bn254_pairing_miller_beuchat(g, Q, P, p);
        }
        t2 = clock();
        printf("optimal ate pairing miller loop: %.5lf [msec]\n", (double)(t2 - t1) / CLOCKS_PER_SEC / N * 1000);
    }
    else
    {
        t1 = clock();
        for (i = 0; i < N; i++) {
            ec_bn254_pairing_miller_aranha_proj(g, Q, P, p);
        }
        t2 = clock();
        printf("optimal ate pairing miller loop (proj): %.5lf [msec]\n", (double)(t2 - t1) / CLOCKS_PER_SEC / N * 1000);

        t1 = clock();
        for (i = 0; i < N; i++) {
            ec_bn254_pairing_miller_aranha_jac(g, Q, P, p);
        }
        t2 = clock();
        printf("optimal ate pairing miller loop (jac): %.5lf [msec]\n", (double)(t2 - t1) / CLOCKS_PER_SEC / N * 1000);
    }
    t1 = clock();
    for (i = 0; i < N; i++) {
        ec_bn254_pairing_finalexp(f, g, p);
    }
    t2 = clock();
    printf("optimal ate pairing final exponentiation: %.5lf [msec]\n", (double)(t2 - t1) / CLOCKS_PER_SEC / N * 1000);

    //-------------------
    //  clear
    //-------------------
    point_clear(P);
    point_clear(P1);
    point_clear(P2);
    point_clear(R);
    point_clear(Q);
    point_clear(Q1);
    point_clear(Q2);
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
//  main program
//============================================
int main(void)
{
    EC_PAIRING pa, pb;

    pairing_init(pa, "ECBN254a");
    test_feature(pa);
    test_pairing(pa);
    pairing_clear(pa);

    pairing_init(pb, "ECBN254b");
    test_feature(pb);
    test_pairing(pb);
    pairing_clear(pb);

    fprintf(stderr, "ok\n");

    return 0;
}
