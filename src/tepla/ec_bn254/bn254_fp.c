//==============================================================
//  prime field ( bn254_fp ) implementation with GMP
//--------------------------------------------------------------
//  2015.10.31 created by kanbara
//==============================================================

#include "ec_bn254_lcl.h"

#define rep(x)  (*((mpz_t *)x->data))
#define order(x) (x->field->order)
#define field(x) (x->field)

//-------------------------------------------
//  initialization, clear, set
//-------------------------------------------
void bn254_fp_init(Element x)
{
    x->data = (void *)malloc(sizeof(mpz_t));

    if (x->data == NULL) {
        fprintf(stderr, "fail: allocate in fp init\n");
        exit(100);
    }

    mpz_init(rep(x));
}

void bn254_fp_clear(Element x)
{
    if (x->data != NULL)
    {
        mpz_clear(rep(x));
        free(x->data);
        x->data = NULL;
    }
}

//-------------------------------------------
//  set value
//-------------------------------------------
void bn254_fp_set(Element x, const Element y)
{
    mpz_set(rep(x), rep(y));
}

void bn254_fp_set_str(Element x, const char *s)
{
    mpz_set_str(rep(x), s, 16);
    mpz_mod(rep(x), rep(x), order(x));
}

void bn254_fp_get_str(char *s, const Element x)
{
    mpz_get_str(s, 16, rep(x));
}

void bn254_fp_set_zero(Element x)
{
    mpz_set_ui(rep(x), 0);
}

void bn254_fp_set_one(Element x)
{
    mpz_set_ui(rep(x), 1);
}

//-------------------------------------------
//  arithmetic operation
//-------------------------------------------
void bn254_fp_add(Element z, const Element x, const Element y)
{
    mpz_add(rep(z), rep(x), rep(y));
    if (mpz_cmp(rep(z), order(z)) >= 0)
    {
        mpz_sub(rep(z), rep(z), order(z));
    }
}

void bn254_fp_addn(Element z, const Element x, const Element y)
{
    mpz_add(rep(z), rep(x), rep(y));
}

void bn254_fp_addp(Element z, const Element x)
{
    mpz_add(rep(z), rep(x), order(x));
}

void bn254_fp_add_one(Element z, const Element x)
{
    mpz_add_ui(rep(z), rep(x), 1);
}

void bn254_fp_neg(Element z, const Element x)
{
    mpz_sub(rep(z), order(z), rep(x));
}

void bn254_fp_sub(Element z, const Element x, const Element y)
{
    mpz_sub(rep(z), rep(x), rep(y));
    if (mpz_sgn(rep(z)) < 0)
    {
        mpz_add(rep(z), rep(z), order(z));
    }
}

void bn254_fp_subn(Element z, const Element x, const Element y)
{
    mpz_sub(rep(z), rep(x), rep(y));
}

void bn254_fp_mul(Element z, const Element x, const Element y)
{
    mpz_mul(rep(z), rep(x), rep(y));
    mpz_mod(rep(z), rep(z), order(z));
}

void bn254_fp_muln(Element z, const Element x, const Element y)
{
    mpz_mul(rep(z), rep(x), rep(y));
}

void bn254_fp_mulc(Element z, const Element x, const mpz_t c)
{
    mpz_mul(rep(z), rep(x), c);
}

void bn254_fp_div2(Element z, const Element x)
{
    mpz_set(rep(z), rep(x));
    if (mpz_odd_p(rep(z)))
    { mpz_add(rep(z), rep(z), order(z)); }
    mpz_fdiv_q_2exp(rep(z), rep(z), 1);
}

void bn254_fp_inv(Element z, const Element x)
{
    mpz_invert(rep(z), rep(x), order(z));
}

void bn254_fp_mod(Element z, const Element x)
{
    mpz_mod(rep(z), rep(x), order(z));
}

void bn254_fp_dob(Element z, const Element x)
{
    bn254_fp_add(z, x, x);
}

void bn254_fp_tri(Element z, const Element x)
{
    Element *t = field(z)->tmp;

    bn254_fp_add(t[0], x, x);
    bn254_fp_add(z, t[0], x);
}

void bn254_fp_sqr(Element z, const Element x)
{
    mpz_mul(rep(z), rep(x), rep(x));
    mpz_mod(rep(z), rep(z), order(z));
}

void bn254_fp_pow(Element z, const Element x, const mpz_t exp)
{
    mpz_powm(rep(z), rep(x), exp, order(z));
}

void lucas_sequence(mpz_t Vk, mpz_t Qk, const mpz_t P, const mpz_t Q, const mpz_t n, const mpz_t k)
{
    mpz_t v0, v1, q0, q1, t;

    long r, i;

    mpz_init_set_ui(v0, 2);
    mpz_init_set(v1, P);
    mpz_init_set_ui(q0, 1);
    mpz_init_set_ui(q1, 1);
    mpz_init(t);

    r = mpz_sizeinbase(k, 2);

    for (i = r - 1; i >= 0; i--)
    {
        mpz_mul(q0, q0, q1);
        mpz_mod(q0, q0, n);

        if (mpz_tstbit(k, i))
        {
            mpz_mul(q1, q0, Q);
            mpz_mod(q0, q0, n);	// q0 = q0*q1 mod n
            mpz_mul(v0, v0, v1);
            mpz_mul(t, P, q0);
            mpz_sub(v0, v0, t);
            mpz_mod(v0, v0, n);	// v0 = v0*v1 - P*q0 mod n
            mpz_mul(v1, v1, v1);
            mpz_add(t, q1, q1);
            mpz_sub(v1, v1, t);
            mpz_mod(v1, v1, n);	// v1 = v1^2 - 2*q1 mod n
        }
        else
        {
            mpz_set(q1, q0);	// q1 = q0
            mpz_mul(v1, v1, v0);
            mpz_mul(t, P, q0);
            mpz_sub(v1, v1, t);
            mpz_mod(v1, v1, n);	// v1 = v0*v1 - P*q0 mod n
            mpz_mul(v0, v0, v0);
            mpz_add(t, q0, q0);
            mpz_sub(v0, v0, t);
            mpz_mod(v0, v0, n);	// v0 = v0^2 - 2*q0
        }
    }

    mpz_set(Vk, v0);
    mpz_set(Qk, q0);

    mpz_clear(v0);
    mpz_clear(v1);
    mpz_clear(q0);
    mpz_clear(q1);
    mpz_clear(t);
}

int bn254_fp_sqrt(Element z, const Element x)
{
    mpz_t V, P, Q, Q0, n, k, z1, z2, i2;

    if (!bn254_fp_is_sqr(x)) {
        return FALSE;
    }

    mpz_init_set(Q, rep(x));
    mpz_init_set_ui(P, 0);
    mpz_init_set(n, order(x));

    mpz_init(k);
    mpz_add_ui(k, n, 1);
    mpz_div_ui(k, k, 2);

    mpz_init(V);
    mpz_init(Q0);
    mpz_init(z1);
    mpz_init(z2);
    mpz_init_set_ui(i2, 2);
    mpz_invert(i2, i2, n);

    while (1)
    {
        mpz_add_ui(P, P, 1);
        lucas_sequence(V, Q0, P, Q, n, k);
        mpz_mul(z1, V, i2);
        mpz_mod(z1, z1, n);
        mpz_mul(z2, z1, z1);
        mpz_mod(z2, z2, n);

        if (mpz_cmp(z2, Q) == 0) {
            break;
        }
    }

    mpz_set(rep(z), z1);

    mpz_clear(Q);
    mpz_clear(P);
    mpz_clear(n);
    mpz_clear(k);
    mpz_clear(V);
    mpz_clear(Q0);
    mpz_clear(z1);
    mpz_clear(z2);
    mpz_clear(i2);

    return TRUE;
}

void bn254_fp_OP1_1(Element z, const Element x)
{
    mpz_add(rep(z), rep(x), field(x)->OP1_1);
}

void bn254_fp_OP1_2(Element z, const Element x)
{
    mpz_add(rep(z), rep(x), field(x)->OP1_2);
}

void bn254_fp_OP2(Element z, const Element x)
{
    mpz_add(rep(z), rep(x), field(x)->OP2);
}

//-------------------------------------------
//  comparison operation
//-------------------------------------------
int bn254_fp_is_zero(const Element x)
{
    return (mpz_sgn(rep(x)) == 0);
}

int bn254_fp_is_one(const Element x)
{
    return (mpz_cmp_ui(rep(x), 1) == 0);
}

int bn254_fp_is_sqr(const Element x)
{
    return (mpz_legendre(rep(x), order(x)) == 1);
}

int bn254_fp_cmp(const Element x, const Element y)
{
    return mpz_cmp(rep(x), rep(y));
}

//-------------------------------------------
//  general function for is sqr
//-------------------------------------------
int bn254_fp_is_sqr_general(const Element x)
{
    int hr = FALSE;

    mpz_t q;
    Element t;

    if (element_is_zero(x)) {
        return FALSE;
    }

    mpz_init(q);
    element_init(t, field(x));

    mpz_sub_ui(q, order(x), 1);
    mpz_tdiv_q_2exp(q, q, 1);

    element_pow(t, x, q);

    hr = element_is_one(t);

    mpz_clear(q);
    element_clear(t);

    return hr;
}

//-------------------------------------------
//  generate random element
//-------------------------------------------
void bn254_fp_random(Element z)
{
    static gmp_randstate_t *state = NULL;

    if (state == NULL)
    {
        state = (gmp_randstate_t *)malloc(sizeof(gmp_randstate_t));
        gmp_randinit_default(*state);
        gmp_randseed_ui(*state, (int)time(NULL));
    }
    mpz_urandomm(rep(z), *state, order(z));
}

//-------------------------------------------
//  i/o operation (octet string)
//-------------------------------------------
void bn254_fp_to_oct(unsigned char *os, size_t *size, const Element x)
{
    size_t stmp;

    unsigned char ostmp[32];

    mpz_export(ostmp, &stmp, 1, sizeof(*ostmp), 1, 0, rep(x));

    memset(os, 0x00, 32);

    memcpy(&(os[32 - (int)stmp]), ostmp, stmp);

    (*size) = 32;
}

void bn254_fp_from_oct(Element x, const unsigned char *os, const size_t size)
{
    if (size < 32) {
        fprintf(stderr, "error: please set up the enought buffer for element\n");
        exit(300);
    }

    mpz_import(rep(x), size, 1, sizeof(*os), 1, 0, os);
}
