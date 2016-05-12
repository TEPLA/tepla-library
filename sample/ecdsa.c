#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <tepla/ec.h>
#include <openssl/sha.h>

#define N SHA_DIGEST_LENGTH

void calc_sha1(mpz_t dst, const unsigned char * str)
{
    unsigned char *buf = (unsigned char *)malloc(sizeof(unsigned char *) * N);
    size_t len = strlen((char*)str);
    SHA_CTX c;
    SHA1_Init(&c);
    SHA1_Update(&c, str, len);
    SHA1_Final(buf, &c);

    mpz_import(dst, N, 1, sizeof(buf[0]), 1, 0, buf);
}

void mpz_rand(mpz_t z, const mpz_t order)
{
    static gmp_randstate_t *state = NULL;

    if (state == NULL)
    {
        state = (gmp_randstate_t *)malloc(sizeof(gmp_randstate_t));
        gmp_randinit_default(*state);
        gmp_randseed_ui(*state, (int)time(NULL));
    }
    mpz_urandomm(z, *state, order);
}

void ecdsa(mpz_t r, mpz_t s, const mpz_t d, const unsigned char * str)
{
    EC_GROUP ec;
    EC_POINT kG, G;
    mpz_t hm, k, q, tmp;

    // init
    curve_init(ec, "ec_bn254_fpb");
    point_init(kG, ec);
    point_init(G, ec);
    point_set(G, ec->generator);
    mpz_init(hm);
    mpz_init(k);
    mpz_init_set(q, *curve_get_order(ec));
    mpz_init(tmp);

    // signature
    mpz_rand(k, q);         // k = [1..q]
    point_mul(kG, k, G);    // kG = k * G
    mpz_set(r, kG->x->data);  // r = kG->x

    calc_sha1(hm, str);     // buf = H(m)

    mpz_mul(tmp, d, r);     // tmp = d * r
    mpz_add(s, hm, tmp);    // s = H(m) + d * r
    mpz_invert(tmp, k, q);  // tmp = k^-1 mod q
    mpz_mul(s, s, tmp);     // s = (H(m) + d * r) * k^-1
    mpz_mod(s, s, q);       // s = (H(m) + d * r) * k^-1 mod q

    gmp_printf("Signature:\n(r,s) = (%Zx,%Zx)\n", r, s);

    // free
    mpz_clear(tmp);
    mpz_clear(q);
    mpz_clear(k);
    mpz_clear(hm);
    point_clear(G);
    point_clear(kG);
    curve_clear(ec);
}

void verify(const mpz_t r, const mpz_t s, const EC_POINT Q, const unsigned char * str)
{
    EC_GROUP ec;
    EC_POINT S, T;
    mpz_t r2, hm, tmp1, tmp2;

    // init
    curve_init(ec, "ec_bn254_fpb");
    point_init(S, ec);
    point_set(S, ec->generator);
    point_init(T, ec);
    mpz_init(r2);
    mpz_init(hm);
    mpz_init(tmp1);
    mpz_init(tmp2);

    mpz_invert(r2, s, *curve_get_order(ec));
    calc_sha1(hm, str);
    mpz_mul(tmp1, r2, hm);  // tmp1 = s^-1 * H(m)
    mpz_mul(tmp2, r2, r);   // tmp2 = s^-1 * r
    point_mul(S, tmp1, S);  // S = (s^-1 * H(m)) * G
    point_mul(T, tmp2, Q);  // T = (s^-1 * r) * Q
    point_add(S, S, T);     // S = (s^-1 * H(m)) * G + (s^-1 * r) * Q

    mpz_set(r2, S->x->data);  // r2 = S->x

    gmp_printf("Verification:\n(r,r2) = (%Zx,%Zx)\n", r, r2);

    if (mpz_cmp(r, r2) == 0)
    {
        printf("Successful Verification!\n");
    }
    else
    {
        printf("The Varidation failed\n");
    }

    // free
    mpz_clear(tmp2);
    mpz_clear(tmp1);
    mpz_clear(r2);
    point_clear(T);
    point_clear(S);
    curve_clear(ec);
}

int main()
{
    EC_GROUP ec;
    EC_POINT G, Q;
    mpz_t r, s, d;
    unsigned char *str = (unsigned char*)"LCIS";

    // init
    curve_init(ec, "ec_bn254_fpb");
    point_init(G, ec);
    point_set(G, ec->generator);
    point_init(Q, ec);
    mpz_init(r);
    mpz_init(s);
    mpz_init(d);

    // signature
    mpz_rand(d, *curve_get_order(ec));  // d is a verification key

    point_mul(Q, d, G);     // Q (= d * G) is a back door

    ecdsa(r, s, d, str);    // create ECDSA

    verify(r, s, Q, str);   // verify ECDSA

    // clear
    mpz_clear(d);
    mpz_clear(s);
    mpz_clear(r);

    point_clear(Q);
    point_clear(G);

    curve_clear(ec);

    return 0;
}
