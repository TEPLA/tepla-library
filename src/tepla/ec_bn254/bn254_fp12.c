//==============================================================
//  prime field ( bn254_fp12 ) implementation with GMP
//--------------------------------------------------------------
//  bn254_fp12a Fp12 := Fp6[z]/(z^2 - y) Beuchat et al.
//  bn254_fp12b Fp12 := Fp6[z]/(z^2 - y) Aranha et al.
//--------------------------------------------------------------
//  2015.10.31 created by kanbara
//==============================================================

#include "ec_bn254_lcl.h"

#define rep(x) (*((mpz_t *)x->data))
#define rep0(x) (((Element *)x->data)[0])
#define rep1(x) (((Element *)x->data)[1])
#define rep2(x) (((Element *)x->data)[2])

#define field(x) (x->field)
#define order(x) (x->field->order)

//-------------------------------------------
//  initialization, clear, set
//-------------------------------------------
void bn254_fp12_init(Element x)
{
    x->data = (void *)malloc(sizeof(Element) * 2);

    if (x->data == NULL) {
        fprintf(stderr, "fail: allocate in bn254_fp12 init\n");
        exit(100);
    }

    element_init(rep0(x), field(x)->base);
    element_init(rep1(x), field(x)->base);
}

void bn254_fp12_clear(Element x)
{
    if (x->data != NULL)
    {
        element_clear(rep0(x));
        element_clear(rep1(x));

        free(x->data);
        x->data = NULL;
    }
}

void bn254_fp12_set(Element x, const Element y)
{
    bn254_fp6_set(rep0(x), rep0(y));
    bn254_fp6_set(rep1(x), rep1(y));
}

void bn254_fp12_set_fp6(Element z, const Element x, const Element y)
{
    bn254_fp6_set(rep0(z), x);
    bn254_fp6_set(rep1(z), y);
}

void bn254_fp12_set_str(Element x, const char *s)
{
    int i = 0;
    int len = strlen(s);

    char msg[790], *p, *c[11];

    if (len > 790) {
        fprintf(stderr, "error: input string is too long, string must be smaller than 780\n");
        exit(200);
    }

    strcpy(msg, s);

    p = msg;

    while ((*p) != '\0')
    {
        if ((*p) == ' ') {
            if (i < 11) {
                c[i] = p;
            }
            i++;
        }
        p++;
    }

    if (i != 11) {
        fprintf(stderr, "error: input string is not correct %d\n", i);
        exit(200);
    }

    (*c[0]) = '\0';
    (*c[1]) = '\0';
    (*c[2]) = '\0';
    (*c[3]) = '\0';
    (*c[4]) = '\0';
    (*c[5]) = '\0';
    (*c[6]) = '\0';
    (*c[7]) = '\0';
    (*c[8]) = '\0';
    (*c[9]) = '\0';
    (*c[10]) = '\0';

    bn254_fp_set_str(rep0(rep0(rep0(x))), msg);
    bn254_fp_set_str(rep0(rep1(rep0(x))), ++c[0]);
    bn254_fp_set_str(rep0(rep2(rep0(x))), ++c[1]);
    bn254_fp_set_str(rep0(rep0(rep1(x))), ++c[2]);
    bn254_fp_set_str(rep0(rep1(rep1(x))), ++c[3]);
    bn254_fp_set_str(rep0(rep2(rep1(x))), ++c[4]);
    bn254_fp_set_str(rep1(rep0(rep0(x))), ++c[5]);
    bn254_fp_set_str(rep1(rep1(rep0(x))), ++c[6]);
    bn254_fp_set_str(rep1(rep2(rep0(x))), ++c[7]);
    bn254_fp_set_str(rep1(rep0(rep1(x))), ++c[8]);
    bn254_fp_set_str(rep1(rep1(rep1(x))), ++c[9]);
    bn254_fp_set_str(rep1(rep2(rep1(x))), ++c[10]);
}

void bn254_fp12_get_str(char *s, const Element x)
{
    char s0[65], s1[65], s2[65], s3[65], s4[65], s5[65], s6[65], s7[65], s8[65], s9[65], s10[65], s11[65];

    bn254_fp_get_str(s0, rep0(rep0(rep0(x))));
    bn254_fp_get_str(s1, rep0(rep1(rep0(x))));
    bn254_fp_get_str(s2, rep0(rep2(rep0(x))));
    bn254_fp_get_str(s3, rep0(rep0(rep1(x))));
    bn254_fp_get_str(s4, rep0(rep1(rep1(x))));
    bn254_fp_get_str(s5, rep0(rep2(rep1(x))));
    bn254_fp_get_str(s6, rep1(rep0(rep0(x))));
    bn254_fp_get_str(s7, rep1(rep1(rep0(x))));
    bn254_fp_get_str(s8, rep1(rep2(rep0(x))));
    bn254_fp_get_str(s9, rep1(rep0(rep1(x))));
    bn254_fp_get_str(s10, rep1(rep1(rep1(x))));
    bn254_fp_get_str(s11, rep1(rep2(rep1(x))));

    sprintf(s, "%s %s %s %s %s %s %s %s %s %s %s %s", s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11);
}

void bn254_fp12_set_zero(Element x)
{
    bn254_fp6_set_zero(rep0(x));
    bn254_fp6_set_zero(rep1(x));
}

void bn254_fp12_set_one(Element x)
{
    bn254_fp6_set_one(rep0(x));
    bn254_fp6_set_zero(rep1(x));
}

//-------------------------------------------
//  arithmetic operation
//-------------------------------------------
void bn254_fp12_add(Element z, const Element x, const Element y)
{
    bn254_fp6_add(rep0(z), rep0(x), rep0(y));
    bn254_fp6_add(rep1(z), rep1(x), rep1(y));
}

void bn254_fp12_neg(Element z, const Element x)
{
    bn254_fp6_neg(rep0(z), rep0(x));
    bn254_fp6_neg(rep1(z), rep1(x));
}

void bn254_fp12_sub(Element z, const Element x, const Element y)
{
    bn254_fp6_sub(rep0(z), rep0(x), rep0(y));
    bn254_fp6_sub(rep1(z), rep1(x), rep1(y));
}

void bn254_fp12_mul(Element z, const Element x, const Element y)
{
    Element *t = field(z)->base->tmp;

    if (strcmp(x->field->field_name, "bn254_fp12a") == 0)
    {
        bn254_fp6_add(t[1], rep0(x), rep1(x));
        bn254_fp6_add(t[2], rep0(y), rep1(y));
        bn254_fp6_mul(t[0], t[1], t[2]);
        bn254_fp6_mul(t[1], rep0(x), rep0(y));
        bn254_fp6_mul(t[2], rep1(x), rep1(y));
        bn254_fp6_sub(t[0], t[0], t[1]);
        bn254_fp6_sub(rep1(z), t[0], t[2]);
        bn254_fp6_gm_mul(t[0], t[2]);
        bn254_fp6_add(rep0(z), t[0], t[1]);
    }

    if (strcmp(x->field->field_name, "bn254_fp12b") == 0)
    {
        bn254_fp6_muln(t[0], rep0(x), rep0(y));
        bn254_fp6_muln(t[1], rep1(x), rep1(y));
        bn254_fp6_add(t[8], rep0(x), rep1(x));
        bn254_fp6_add(t[9], rep0(y), rep1(y));
        bn254_fp6_muln(t[2], t[8], t[9]);
        bn254_fp6_add(t[3], t[0], t[1]);
        bn254_fp6_OP2(t[3], t[3]);
        bn254_fp6_sub(t[2], t[2], t[3]);
        bn254_fp6_OP2(t[2], t[2]);
        bn254_fp6_mod(rep1(z), t[2]);
        bn254_fp6_gm_mul(t[2], t[1]);
        bn254_fp6_add(t[2], t[0], t[2]);
        bn254_fp6_OP2(t[2], t[2]);
        bn254_fp6_mod(rep0(z), t[2]);
    }
}

//----------------------------------------------------------
//  Input  : z in Fp12 and l0, l3, l4 in Fp2
//  Output : z *= { (x0, 0, 0), (x1, x2, 0) } in Fp12
//----------------------------------------------------------
void bn254_fp12_mul_L(Element z, Element x0, Element x1, Element x2)
{
    Element *v = field(z)->base->tmp;

    bn254_fp6_mul_fp2(v[0], rep0(z), x0);
    bn254_fp6_mul_fp2_2(v[1], rep1(z), x1, x2);
    bn254_fp6_add(rep1(z), rep1(z), rep0(z));
    bn254_fp2_add(x0, x0, x1);
    bn254_fp6_mul_fp2_2(rep1(z), rep1(z), x0, x2);
    bn254_fp6_gm_mul(rep0(z), v[1]);
    bn254_fp6_add(rep0(z), rep0(z), v[0]);
    bn254_fp6_sub(rep1(z), rep1(z), v[0]);
    bn254_fp6_sub(rep1(z), rep1(z), v[1]);
}

//----------------------------------------------------------
//  Input  : z in Fp12 and l0, l2, l4 in Fp2
//  Output : z *= { (x0, 0, x1), (0, x2, 0) } in Fp12
//----------------------------------------------------------
void bn254_fp12_mul_L2(Element z, Element x0, Element x1, Element x2)
{
    Element *v = field(z)->base->tmp;

    bn254_fp6_mul_fp2_4(v[0], rep0(z), x0, x1);
    bn254_fp6_mul_fp2_3(v[1], rep1(z), x2);

    bn254_fp6_add(rep1(z), rep1(z), rep0(z));
    bn254_fp6_set_fp2(v[2], x0, x2, x1);
    bn254_fp6_mul(rep1(z), rep1(z), v[2]);

    bn254_fp6_sub(rep1(z), rep1(z), v[0]);
    bn254_fp6_sub(rep1(z), rep1(z), v[1]);
    bn254_fp6_gm_mul(rep0(z), v[1]);
    bn254_fp6_add(rep0(z), rep0(z), v[0]);
}

void bn254_fp12_inv(Element z, const Element x)
{
    Element *t = field(z)->base->tmp;

    bn254_fp6_sqr(t[0], rep0(x));  // t0 = a0^2
    bn254_fp6_sqr(t[1], rep1(x));  // t1 = a1^2

    bn254_fp6_gm_mul(t[2], t[1]);    //
    bn254_fp6_sub(t[0], t[0], t[2]); // t0 = t0 - t1*gamma
    bn254_fp6_inv(t[1], t[0]);       // t1 = t0^-1

    bn254_fp6_mul(rep0(z), rep0(x), t[1]); // c0 = a0*t1
    bn254_fp6_mul(rep1(z), rep1(x), t[1]); // c1 = a1*t1
    bn254_fp6_neg(rep1(z), rep1(z));       // c1 = -c1
}

void bn254_fp12_dob(Element z, const Element x)
{
    bn254_fp6_dob(rep0(z), rep0(x));
    bn254_fp6_dob(rep1(z), rep1(x));
}

void bn254_fp12_tri(Element z, const Element x)
{
    bn254_fp6_tri(rep0(z), rep0(x));
    bn254_fp6_tri(rep1(z), rep1(x));
}

void bn254_fp12_sqr(Element z, const Element x)
{
    Element *t = field(z)->base->tmp;

    if (strcmp(x->field->field_name, "bn254_fp12a") == 0)
    {
        bn254_fp6_sub(t[0], rep0(x), rep1(x));     // t0 = a0 - a1
        bn254_fp6_gm_mul(t[1], rep1(x));           //
        bn254_fp6_sub(t[1], rep0(x), t[1]);        // t1 = a0 - gamma*a1
        bn254_fp6_mul(rep0(z), rep0(x), rep1(x));  // c0 = a0*a1
        bn254_fp6_mul(t[0], t[0], t[1]);           //
        bn254_fp6_add(t[0], t[0], rep0(z));        // t0 = t0*t1 + c0
        bn254_fp6_add(rep1(z), rep0(z), rep0(z));  // c1 = 2*a0*a1
        bn254_fp6_gm_mul(t[1], rep0(z));    // t1 = gamma*t1
        bn254_fp6_add(rep0(z), t[0], t[1]); // c0 = t1 + t0
    }

    if (strcmp(x->field->field_name, "bn254_fp12b") == 0)
    {
        bn254_fp6_add(t[0], rep0(x), rep1(x));
        bn254_fp6_gm_mul(t[1], rep1(x));
        bn254_fp6_add(t[1], t[1], rep0(x));
        bn254_fp6_muln(t[3], rep0(x), rep1(x));
        bn254_fp6_mod(rep1(z), t[3]);
        bn254_fp6_muln(t[4], t[0], t[1]);
        bn254_fp6_mod(t[0], t[4]);
        bn254_fp6_gm_mul(t[1], rep1(z));
        bn254_fp6_add(t[1], t[1], rep1(z));
        bn254_fp6_sub(rep0(z), t[0], t[1]);
        bn254_fp6_dob(rep1(z), rep1(z));
    }
}

//-----------------------------------------------------------
//  exponentiation z = x^exp
//-----------------------------------------------------------
void bn254_fp12_pow(Element z, const Element x, const mpz_t exp)
{
    long t, i;
    Element c;

    element_init(c, field(z));
    element_set(c, x);

    t = mpz_sizeinbase(exp, 2);

    for (i = t - 2; i >= 0; i--)
    {
        element_sqr(c, c);
        if (mpz_tstbit(exp, i)) {
            element_mul(c, c, x);
        }
    }

    element_set(z, c);
    element_clear(c);
}

//-----------------------------------------------------------
//  exponentiation z = x^exp with NAF
//-----------------------------------------------------------
void bn254_fp12_pow_naf(Element z, const Element x, const mpz_t exp)
{
    long t, i;
    Element c, ix;

    int *naf, nlen;

    element_init(c, field(z));
    element_init(ix, field(z));

    element_set(c, x);
    element_inv(ix, x);

    t = mpz_sizeinbase(exp, 2);

    naf = (int *)malloc(sizeof(int) * (t + 1));

    generate_naf(naf, &nlen, exp);

    for (i = nlen - 2; i >= 0; i--)
    {
        element_sqr(c, c);
        if (naf[i])
        {
            if (naf[i] < 0) {
                element_mul(c, c, ix);
            }
            else {
                element_mul(c, c, x);
            }
        }
    }

    element_set(z, c);
    element_clear(c);
    element_clear(ix);

    free(naf);
}

//-----------------------------------------------------------
//  Frobenius Map in Fp12
//-----------------------------------------------------------
//  frob(x) == x^p
//  z = g + h*w  : g = g0 + g1*v + g2*v^2
//               : h = h0 + h1*v + h2*v^2
//-----------------------------------------------------------
void bn254_fp12_frob_p(Element z, const Element x)
{
    field_precomp_frob_p pf;

    pf = ((field_precomp_p)(field(z)->precomp))->pf;

    bn254_fp6_conj(rep0(z), rep0(x));
    bn254_fp6_conj(rep1(z), rep1(x));

    if (strcmp(x->field->field_name, "bn254_fp12a") == 0)
    {
        bn254_fp2_mul_p(rep0(rep1(z)), rep0(rep1(z)), (pf->gamma1)[0]);   //t2 = t2*gamma1
        bn254_fp2_mul_p(rep1(rep0(z)), rep1(rep0(z)), (pf->gamma1)[1]);   //t3 = t3*gamma2
        bn254_fp2_mul_p(rep1(rep1(z)), rep1(rep1(z)), (pf->gamma1)[2]);   //t4 = t4*gamma3
        bn254_fp2_mul_p(rep2(rep0(z)), rep2(rep0(z)), (pf->gamma1)[3]);   //t5 = t5*gamma4
        bn254_fp2_mul_p(rep2(rep1(z)), rep2(rep1(z)), (pf->gamma1)[4]);   //t6 = t6*gamma5
    }

    if (strcmp(x->field->field_name, "bn254_fp12b") == 0)
    {
        bn254_fp2_mul(rep0(rep1(z)), rep0(rep1(z)), (pf->gamma1)[0]);   //t2 = t2*gamma1
        bn254_fp2_mul(rep1(rep0(z)), rep1(rep0(z)), (pf->gamma1)[1]);   //t3 = t3*gamma2
        bn254_fp2_mul(rep1(rep1(z)), rep1(rep1(z)), (pf->gamma1)[2]);   //t4 = t4*gamma3
        bn254_fp2_mul(rep2(rep0(z)), rep2(rep0(z)), (pf->gamma1)[3]);   //t5 = t5*gamma4
        bn254_fp2_mul(rep2(rep1(z)), rep2(rep1(z)), (pf->gamma1)[4]);   //t6 = t6*gamma5
    }
}

void bn254_fp12_frob_p2(Element z, const Element x)
{
    field_precomp_frob_p pf;

    pf = ((field_precomp_p)(field(z)->precomp))->pf;

    bn254_fp2_set(rep0(rep0(z)), rep0(rep0(x)));   //t1 = g0;

    if (strcmp(x->field->field_name, "bn254_fp12a") == 0)
    {
        bn254_fp2_mul_p(rep0(rep1(z)), rep0(rep1(x)), (pf->gamma2)[0]);   //t2 = h0*gamma1
        bn254_fp2_mul_p(rep1(rep0(z)), rep1(rep0(x)), (pf->gamma2)[1]);   //t3 = g1*gamma2
        bn254_fp2_mul_p(rep1(rep1(z)), rep1(rep1(x)), (pf->gamma2)[2]);   //t4 = h1*gamma3
        bn254_fp2_mul_p(rep2(rep0(z)), rep2(rep0(x)), (pf->gamma2)[3]);   //t5 = g2*gamma4
        bn254_fp2_mul_p(rep2(rep1(z)), rep2(rep1(x)), (pf->gamma2)[4]);   //t6 = h2*gamma5
    }

    if (strcmp(x->field->field_name, "bn254_fp12b") == 0)
    {
        bn254_fp2_mul(rep0(rep1(z)), rep0(rep1(x)), (pf->gamma2)[0]);   //t2 = h0*gamma1
        bn254_fp2_mul(rep1(rep0(z)), rep1(rep0(x)), (pf->gamma2)[1]);   //t3 = g1*gamma2
        bn254_fp2_mul(rep1(rep1(z)), rep1(rep1(x)), (pf->gamma2)[2]);   //t4 = h1*gamma3
        bn254_fp2_mul(rep2(rep0(z)), rep2(rep0(x)), (pf->gamma2)[3]);   //t5 = g2*gamma4
        bn254_fp2_mul(rep2(rep1(z)), rep2(rep1(x)), (pf->gamma2)[4]);   //t6 = h2*gamma5
    }
}

void bn254_fp12_frob_p3(Element z, const Element x)
{
    field_precomp_frob_p pf;

    pf = ((field_precomp_p)(field(z)->precomp))->pf;

    bn254_fp6_conj(rep0(z), rep0(x));
    bn254_fp6_conj(rep1(z), rep1(x));

    if (strcmp(x->field->field_name, "bn254_fp12a") == 0)
    {
        bn254_fp2_mul_p(rep0(rep1(z)), rep0(rep1(z)), (pf->gamma3)[0]);   //t2 = t2*gamma1
        bn254_fp2_mul_p(rep1(rep0(z)), rep1(rep0(z)), (pf->gamma3)[1]);   //t3 = t3*gamma2
        bn254_fp2_mul_p(rep1(rep1(z)), rep1(rep1(z)), (pf->gamma3)[2]);   //t4 = t4*gamma3
        bn254_fp2_mul_p(rep2(rep0(z)), rep2(rep0(z)), (pf->gamma3)[3]);   //t5 = t5*gamma4
        bn254_fp2_mul_p(rep2(rep1(z)), rep2(rep1(z)), (pf->gamma3)[4]);   //t6 = t6*gamma5
    }

    if (strcmp(x->field->field_name, "bn254_fp12b") == 0)
    {
        bn254_fp2_mul(rep0(rep1(z)), rep0(rep1(z)), (pf->gamma3)[0]);   //t2 = t2*gamma1
        bn254_fp2_mul(rep1(rep0(z)), rep1(rep0(z)), (pf->gamma3)[1]);   //t3 = t3*gamma2
        bn254_fp2_mul(rep1(rep1(z)), rep1(rep1(z)), (pf->gamma3)[2]);   //t4 = t4*gamma3
        bn254_fp2_mul(rep2(rep0(z)), rep2(rep0(z)), (pf->gamma3)[3]);   //t5 = t5*gamma4
        bn254_fp2_mul(rep2(rep1(z)), rep2(rep1(z)), (pf->gamma3)[4]);   //t6 = t6*gamma5
    }
}

void bn254_fp12_conj(Element z, const Element x)
{
    bn254_fp6_set(rep0(z), rep0(x));
    bn254_fp6_neg(rep1(z), rep1(x));
}

//------------------------------------------------------------
//  square operation of Fp^4 (for bn254_sqr_forpairing )
//------------------------------------------------------------
void bn254_fp4_sqr(Element c0, Element c1, const Element a0, const Element a1)
{
    Element *t = field(c0)->tmp;

    bn254_fp2_sqr(t[0], a0);   //t0 = a0^2
    bn254_fp2_sqr(t[1], a1);   //t1 = a1^2

    bn254_fp2_xi_mul(c0, t[1]);  //c0 = t1*xi
    bn254_fp2_add(c0, c0, t[0]); //c0 = c0 + t0

    bn254_fp2_add(c1, a0, a1);   //c1 = a0 + a1
    bn254_fp2_sqr(c1, c1);       //c1 = c1^2
    bn254_fp2_sub(c1, c1, t[0]); //c1 = c1 - t0
    bn254_fp2_sub(c1, c1, t[1]); //c1 = c1 - t1
}

//------------------------------------------------------------------
//  Special exponentiation in Fp12 by Karabina
//    input  : element x
//    output : z = x^t : t = -(2^62+2^55+1)
//-------------------------------------------------------------------

void bn254_fp12_sqr_forpairing_karabina(Element z, const Element x)
{

    Element *T = field(z)->base->base->tmp;
    //------------------------------------------------
    //	g = (g0+g1*s)+(g2+g3*s)*t+(g4+g5*s)t^2
    //	g^2 = h = (h0+h1*s)+(h2+h3*s)*t+(h4+h5*s)t^2
    //	C(g) = [g2,g3,g4,g5]
    //	D(g^2) = (h0+h1*s)+(h2+h3*s)*t+(h4+h5*s)t^2
    //	if h2 neq 0
    //	h1 = (h5^2*xi+3*h4^2-2*h3)/4x2, x0 = (2*h1^2+h2h5-h3h4)*xi+1
    //	if h2 eq 0
    //	h1 = 2*h4h5/h3, h0 = (2*h1^2-3*h3h4)*xi+1
    // 	h2 = 2*g2+3(S4,5-S4,S5)*xi
    // 	h3 = 3(S4+S5*xi)-2*g3
    //	h4 = 3(S2+S3*xi)-2*g4
    // 	h5 = 2*g5+3(S2,3-S2-S3)
    //	Si,j = (gi+gj)^2
    //	Si = gi^2
    //------------------------------------------------

    // g2 = rep0(rep1(x))
    // g3 = rep2(rep0(x))
    // g4 = rep1(rep0(x))
    // g5 = rep2(rep1(x))

    bn254_fp2_sqr(T[0], rep1(rep0(x)));						// T0 = g4^2
    bn254_fp2_sqr(T[1], rep2(rep1(x)));						// T1 = g5^2
    bn254_fp2_add(T[5], rep1(rep0(x)), rep2(rep1(x)));		// t0 = g4+g5
    bn254_fp2_sqr(T[2], T[5]);								// T2 = t0^2
    bn254_fp2_add(T[3], T[0], T[1]);						// T3 = T0+T1
    bn254_fp2_sub(T[3], T[2], T[3]);						// T3 = T2-T3
    bn254_fp2_mod(T[5], T[3]);								// t0 = T3 mod p
    bn254_fp2_add(T[6], rep0(rep1(x)), rep2(rep0(x)));		// t1 = g2+g3
    bn254_fp2_sqr(T[3], T[6]);								// T3 = t1^2
    bn254_fp2_sqr(T[2], rep0(rep1(x)));						// T2 = g2^2
    bn254_fp2_xi_mul(T[6], T[5]);							// t1 = t0*xi
    bn254_fp2_add(T[5], T[6], rep0(rep1(x)));				// t0 = t1+g2
    bn254_fp2_dob(T[5], T[5]);								// t0 = 2*t0
    bn254_fp2_add(rep0(rep1(z)), T[5], T[6]);				// c2 = t0+t1

    bn254_fp2_xi_mul(T[4], T[1]);							// T4 = T1*xi
    bn254_fp2_add(T[4], T[0], T[4]);						// T4 = T0+T4
    bn254_fp2_mod(T[5], T[4]);								// t0 = T4 mod p
    bn254_fp2_sub(T[6], T[5], rep2(rep0(x)));				// t1 = t0-g3
    bn254_fp2_dob(T[6], T[6]);								// t1 = 2*t1

    bn254_fp2_sqr(T[1], rep2(rep0(x)));						// T1 = g3^2
    bn254_fp2_add(rep2(rep0(z)), T[6], T[5]);				// c3 = t1+t0

    bn254_fp2_xi_mul(T[4], T[1]);							// T4 = T1*xi
    bn254_fp2_add(T[4], T[2], T[4]);						// T4 = T2+T4
    bn254_fp2_mod(T[5], T[4]);								// t0 = T4 mod p
    bn254_fp2_sub(T[6], T[5], rep1(rep0(x)));				// t1 = t0-g4
    bn254_fp2_dob(T[6], T[6]);								// t1 = 2*t1
    bn254_fp2_add(rep1(rep0(z)), T[6], T[5]);				// c4 = t1+t0

    bn254_fp2_addn(T[0], T[2], T[1]);						// T0 = T2+T1
    bn254_fp2_sub(T[3], T[3], T[0]);						// T3 = T3-T0
    bn254_fp2_mod(T[5], T[3]);								// t0 = T3 mod p
    bn254_fp2_add(T[6], T[5], rep2(rep1(x)));				// t1 = t0+g5
    bn254_fp2_dob(T[6], T[6]);								// t1 = 2*t1
    bn254_fp2_add(rep2(rep1(z)), T[6], T[5]);				// c5 = t1+t0
}

void bn254_fp12_decompose_forpairing_karabina(Element z, const Element x)
{
    Element *t = field(z)->base->base->tmp;

    // g2 = rep0(rep1(x))
    // g3 = rep2(rep0(x))
    // g4 = rep1(rep0(x))
    // g5 = rep2(rep1(x))

    if (bn254_fp2_is_zero(rep0(rep1(x))))
    {
        bn254_fp2_mul(t[1], rep1(rep0(x)), rep2(rep1(x)));	// t1 = g4*g5
        bn254_fp2_dob(t[1], t[1]);							// t1 = 2g4g5
        bn254_fp2_inv(t[2], rep2(rep0(x)));					// t2 = 1/g3
        bn254_fp2_mul(t[7], t[1], t[2]); 					// c1 = 2g4g5/g3

        bn254_fp2_sqr(t[2], t[1]);							// t2 = g1^2
        bn254_fp2_dob(t[2], t[2]);							// t2 = 2t2
        bn254_fp2_mul(t[3], rep2(rep0(x)), rep1(rep0(x)));	// t3 = g3*g4
        bn254_fp2_dob(t[4], t[3]);							// t4 = 2t3
        bn254_fp2_add(t[4], t[4], t[3]);					// t4 = 3t3
        bn254_fp2_sub(t[0], t[3], t[4]);					// g0 = t3-t4
        bn254_fp2_xi_mul(t[0], t[0]);						// t0 = t0*xi
        bn254_fp2_add_one(t[0], t[0]);						// t0 = t0+1
    }
    else
    {
        bn254_fp2_sqr(t[5], rep2(rep1(x)));					// t5 = g5^2
        bn254_fp2_xi_mul(t[5], t[5]);						// t5 = g5^2*xi
        bn254_fp2_sqr(t[4], rep1(rep0(x)));					// t4 = g4^2
        bn254_fp2_dob(t[6], t[4]);							// t4 = 2*g4^2
        bn254_fp2_add(t[4], t[6], t[4]);					// t4 = 3*g4^2
        bn254_fp2_dob(t[3], rep2(rep0(x)));					// t3 = 2*g3
        bn254_fp2_dob(t[2], rep0(rep1(x)));					// t2 = 2*g2
        bn254_fp2_dob(t[2], t[2]);							// t2 = 4*g2
        bn254_fp2_inv(t[2], t[2]);							// t2 = 1/4*g2
        bn254_fp2_add(t[1], t[5], t[4]);					// t1 = t5+t4
        bn254_fp2_sub(t[1], t[1], t[3]);					// t1 = t5+t4-t3
        bn254_fp2_mul(t[7], t[1], t[2]);					// c1 = (t5+t4+t3)*t2

        bn254_fp2_sqr(t[1], t[7]);							// t1 = t1^2
        bn254_fp2_dob(t[1], t[1]);							// t1 = 2*t1
        bn254_fp2_mul(t[2], rep0(rep1(x)), rep2(rep1(x)));	// t2 = g2*g5
        bn254_fp2_mul(t[3], rep2(rep0(x)), rep1(rep0(x)));	// t3 = g3*g4
        bn254_fp2_dob(t[4], t[3]);							// t4 = 2*t3
        bn254_fp2_add(t[3], t[4], t[3]);					// t3 = 3*t3
        bn254_fp2_add(t[0], t[1], t[2]);					// t0 = t1+t2
        bn254_fp2_sub(t[0], t[0], t[3]); 					// t0 = t0+t3
        bn254_fp2_xi_mul(t[0], t[0]);						// t0 = t0*xi
        bn254_fp2_add_one(t[0], t[0]); 						// c0 = t0+1
    }

    bn254_fp2_set(rep0(rep0(z)), t[0]);					// z0 = c0
    bn254_fp2_set(rep1(rep0(z)), rep1(rep0(x)));		// z1 = c4
    bn254_fp2_set(rep2(rep0(z)), rep2(rep0(x)));		// z2 = c3
    bn254_fp2_set(rep0(rep1(z)), rep0(rep1(x)));		// z3 = c2
    bn254_fp2_set(rep1(rep1(z)), t[7]);					// z4 = c1
    bn254_fp2_set(rep2(rep1(z)), rep2(rep1(x)));		// z5 = c5

}

void bn254_fp12_pow_forpairing_karabina(Element z, const Element x, const int *t, int tlen)
{
    int i;

    Element z55, z62;
    element_init(z55, field(z));
    element_init(z62, field(z));

    if (z == x)
    {
        fprintf(stderr, "error: bad input for bn254_fp12_pow_forpairing\n");
        exit(200);
    }

    element_set(z, x);
    for (i = 0; i < tlen - 1; i++)
    {
        bn254_fp12_sqr_forpairing_karabina(z, z);
        if (i == 55 - 1) {
            bn254_fp12_decompose_forpairing_karabina(z55, z);
        }
    }
    bn254_fp12_decompose_forpairing_karabina(z62, z);

    bn254_fp12_mul(z, x, z55);
    bn254_fp12_mul(z, z, z62);

    bn254_fp12_conj(z, z);

    element_clear(z55);
    element_clear(z62);
}

//------------------------------------------------------------------
//  Special exponentiation in Fp12 by Beuchat
//    input  : element x
//    output : z = x^t : t = 2^62-2^54+2^44
//-------------------------------------------------------------------

void bn254_fp12_sqr_forpairing_beuchat(Element z, const Element x)
{
    Element *t = field(z)->base->base->tmp;
    Element *c = field(z)->base->tmp;

    //------------------------
    // z = g + h*w
    // g = g0 + g1*v + g2*v^2
    // h = h0 + h1*v + h2*v^2
    //------------------------
    bn254_fp4_sqr(t[2], t[3], rep0(rep0(x)), rep1(rep1(x)));   //t00, t11 = (g0 + h1*V)^2
    bn254_fp4_sqr(t[4], t[5], rep0(rep1(x)), rep2(rep0(x)));   //t01, t12 = (h0 + g2*V)^2
    bn254_fp4_sqr(t[6], t[7], rep1(rep0(x)), rep2(rep1(x)));   //t02, t10 = (g1 + h2*V)^2

    bn254_fp2_xi_mul(t[7], t[7]);       // t10 = t10*xi

    bn254_fp2_tri(t[2], t[2]);          // c00 = 3*t00
    bn254_fp2_dob(t[8], rep0(rep0(x))); // tmp = 2*g0
    bn254_fp2_sub(t[2], t[2], t[8]);    // c00 = c00 - tmp

    bn254_fp2_tri(t[4], t[4]);          // c01 = 3*t01
    bn254_fp2_dob(t[8], rep1(rep0(x))); // tmp = 2*g1
    bn254_fp2_sub(t[4], t[4], t[8]);    // c01 = c01 - tmp

    bn254_fp2_tri(t[6], t[6]);          // c02 = 3*t02
    bn254_fp2_dob(t[8], rep2(rep0(x))); // tmp = 2*g2
    bn254_fp2_sub(t[6], t[6], t[8]);    // c02 = c02 - tmp

    bn254_fp2_tri(t[7], t[7]);          // c10 = 3*t10
    bn254_fp2_dob(t[8], rep0(rep1(x))); // tmp = 2*h0
    bn254_fp2_add(t[7], t[7], t[8]);    // c10 = c10 + tmp

    bn254_fp2_tri(t[3], t[3]);          // c11 = 3*t11
    bn254_fp2_dob(t[8], rep1(rep1(x))); // tmp = 2*h1
    bn254_fp2_add(t[3], t[3], t[8]);    // c11 = c11 + tmp

    bn254_fp2_tri(t[5], t[5]);          // c12 = 3*t12
    bn254_fp2_dob(t[8], rep2(rep1(x))); // tmp = 2*h2
    bn254_fp2_add(t[5], t[5], t[8]);    // c12 = c12 + tmp

    bn254_fp6_set_fp2(c[0], t[2], t[4], t[6]); // c0 = c00 + c01*v + c02*v^2
    bn254_fp6_set_fp2(c[1], t[7], t[3], t[5]); // c1 = c10 + c11*v + c12*v^2

    bn254_fp12_set_fp6(z, c[0], c[1]);   // z = c0 + c1*w
}

void bn254_fp12_pow_forpairing_beuchat(Element z, const Element x, const int *t, int tlen)
{
    int i;
    Element ix;

    if (z == x)
    {
        fprintf(stderr, "error: bad input for bn254_fp12_pow_forpairing\n");
        exit(200);
    }

    element_init(ix, field(x));

    bn254_fp12_set(z, x);
    bn254_fp12_conj(ix, x);

    for (i = tlen - 2; i >= 0; i--)
    {
        bn254_fp12_sqr_forpairing_beuchat(z, z);

        if (t[i])
        {
            if (t[i] < 0) {
                bn254_fp12_mul(z, z, ix);
            }
            else {
                bn254_fp12_mul(z, z, x);
            }
        }
    }

    element_clear(ix);
}

//---------------------------------------------------------
//  precomputation for Fp12 frobenius
//---------------------------------------------------------
void bn254_fp12_precomp_frob_aranha(field_precomp_frob_p pf, const Field f)
{
    int i;
    mpz_t  t, u;
    Element xi, tmp1, tmp2;

    //---------------------------------
    //  allocate buffer
    //---------------------------------
    Element *g1 = (Element *)malloc(sizeof(Element) * 5);
    Element *g2 = (Element *)malloc(sizeof(Element) * 5);
    Element *g3 = (Element *)malloc(sizeof(Element) * 5);

    struct ec_field_st *fp2 = f->base->base;

    //---------------------------------
    //  initialization
    //---------------------------------
    for (i = 0; i < 5; i++)
    {
        element_init(g1[i], fp2);
        element_init(g2[i], fp2);
        element_init(g3[i], fp2);
        element_init(tmp2, fp2);
    }

    //---------------------------------
    //  xi = xi^{(p-1)/6}
    //---------------------------------
    mpz_init(t);
    mpz_init(u);
    mpz_sub_ui(t, fp2->base->order, 1);
    mpz_fdiv_q_ui(t, t, 6);

    element_init(xi, fp2);
    element_init(tmp1, fp2);
    element_set_str(xi, "1 1");
    element_pow(tmp1, xi, t);

    //---------------------------------
    //  set gamma 1, 2, 3
    //---------------------------------
    element_set(g1[0], tmp1);

    for (i = 1; i < 5; i++) {
        element_mul(g1[i], g1[i - 1], g1[0]);
    }
    for (i = 0; i < 5; i++)
    {
        bn254_fp2_conj(tmp2, g1[i]);
        element_mul(g2[i], g1[i], tmp2);
    }
    for (i = 0; i < 5; i++) {
        element_mul(g3[i], g1[i], g2[i]);
    }

    pf->gamma1 = g1;
    pf->gamma2 = g2;
    pf->gamma3 = g3;

    pf->glen1 = pf->glen2 = pf->glen3 = 5;

    mpz_clear(t);
    mpz_clear(u);

    element_clear(xi);
    element_clear(tmp1);
    element_clear(tmp2);
}

void bn254_fp12_precomp_frob_beuchat(field_precomp_frob_p pf, const Field f)
{
    int i;
    mpz_t  t, u;
    Element xi, tmp;

    //---------------------------------
    //  allocate buffer
    //---------------------------------
    Element *g1 = (Element *)malloc(sizeof(Element) * 5);
    Element *g2 = (Element *)malloc(sizeof(Element) * 5);
    Element *g3 = (Element *)malloc(sizeof(Element) * 5);

    struct ec_field_st *fp  = f->base->base->base;
    struct ec_field_st *fp2 = f->base->base;

    //---------------------------------
    //  initialization
    //---------------------------------
    for (i = 0; i < 5; i++)
    {
        element_init(g1[i], fp);
        element_init(g2[i], fp);
        element_init(g3[i], fp);
    }

    //---------------------------------
    //  xi = xi^{(p-1)/6}
    //---------------------------------
    mpz_init(t);
    mpz_init(u);
    mpz_sub_ui(t, fp->order, 1);
    mpz_fdiv_q_ui(t, t, 6);

    element_init(xi, fp2);
    element_init(tmp, fp2);
    element_set_str(xi, "0 1");
    element_pow(tmp, xi, t);

    //---------------------------------
    //  set gamma 1, 2, 3
    //---------------------------------
    element_set(g1[0], ((Element*)tmp->data)[0]);

    for (i = 1; i < 5; i++) {
        element_mul(g1[i], g1[i - 1], g1[0]);
    }
    for (i = 0; i < 5; i++) {
        element_mul(g2[i], g1[i], g1[i]);
    }
    for (i = 0; i < 5; i++) {
        element_mul(g3[i], g1[i], g2[i]);
    }

    pf->gamma1 = g1;
    pf->gamma2 = g2;
    pf->gamma3 = g3;

    pf->glen1 = pf->glen2 = pf->glen3 = 5;

    mpz_clear(t);
    mpz_clear(u);

    element_clear(xi);
    element_clear(tmp);
}

//---------------------------------------------------------
// precomputation for Fp12 operation
//---------------------------------------------------------
void bn254_fp12_precomp(Field f)
{
    field_precomp_p precomp = NULL;

    precomp = (field_precomp_p)malloc(sizeof(struct ec_field_precomp_st));

    precomp->ps = (field_precomp_sqrt_p)malloc(sizeof(struct ec_field_precomp_sqrt_st));
    bn254_fp2_precomp_sqrt(precomp->ps, f);

    precomp->pf = (field_precomp_frob_p)malloc(sizeof(struct ec_field_precomp_frob_st));

    if (strcmp(f->field_name, "bn254_fp12a") == 0)
    {
        bn254_fp12_precomp_frob_beuchat(precomp->pf, f);
    }

    if (strcmp(f->field_name, "bn254_fp12b") == 0)
    {
        bn254_fp12_precomp_frob_aranha(precomp->pf, f);
    }

    f->precomp = (void *)precomp;
}

//-------------------------------------------
//  comparison operation
//-------------------------------------------
int bn254_fp12_is_zero(const Element x)
{
    if (bn254_fp6_is_zero(rep1(x)))
    {
        if (bn254_fp6_is_zero(rep0(x))) {
            return TRUE;
        }
    }
    return FALSE;
}

int bn254_fp12_is_one(const Element x)
{
    if (bn254_fp6_is_zero(rep1(x)))
    {
        if (bn254_fp6_is_one(rep0(x))) {
            return TRUE;
        }
    }
    return FALSE;
}

int bn254_fp12_cmp(const Element x, const Element y)
{
    if (bn254_fp6_cmp(rep1(x), rep1(y)) == 0)
    {
        if (bn254_fp6_cmp(rep0(x), rep0(y)) == 0) {
            return 0;
        }
    }
    return 1;
}

int bn254_fp12_is_sqr(const Element x)
{
    int hr = FALSE;

    Element *t = field(x)->base->tmp;

    if (element_is_zero(x)) {
        return FALSE;
    }

    bn254_fp6_inv(t[0], rep1(x));
    bn254_fp6_mul(t[0], t[0], rep0(x));
    bn254_fp6_sqr(t[0], t[0]);
    bn254_fp6_add(t[0], t[0], field(x)->irre_poly[0]);

    hr = bn254_fp6_is_sqr(t[0]);

    return hr;
}

//-------------------------------------------
//  generate random element
//-------------------------------------------
void bn254_fp12_random(Element z)
{
    bn254_fp6_random(rep0(z));
    bn254_fp6_random(rep1(z));
}

//-------------------------------------------
//  i/o operation (octet string)
//-------------------------------------------
void bn254_fp12_to_mpz(mpz_t a, const Element x)
{
    mpz_mul(a, rep(rep1(rep2(rep1(x)))), field(x)->base->base->base->order);   // a = rep121*p
    mpz_add(a, a, rep(rep1(rep1(rep1(x)))));   // a = a + rep111
    mpz_mul(a, a, field(x)->base->base->base->order);   //a = a*p
    mpz_add(a, a, rep(rep1(rep0(rep1(x)))));
    mpz_mul(a, a, field(x)->base->base->base->order);
    mpz_add(a, a, rep(rep1(rep2(rep0(x)))));
    mpz_mul(a, a, field(x)->base->base->base->order);
    mpz_add(a, a, rep(rep1(rep1(rep0(x)))));
    mpz_mul(a, a, field(x)->base->base->base->order);
    mpz_add(a, a, rep(rep1(rep0(rep0(x)))));
    mpz_mul(a, a, field(x)->base->base->base->order);
    mpz_add(a, a, rep(rep0(rep2(rep1(x)))));
    mpz_mul(a, a, field(x)->base->base->base->order);
    mpz_add(a, a, rep(rep0(rep1(rep1(x)))));
    mpz_mul(a, a, field(x)->base->base->base->order);
    mpz_add(a, a, rep(rep0(rep0(rep1(x)))));
    mpz_mul(a, a, field(x)->base->base->base->order);
    mpz_add(a, a, rep(rep0(rep2(rep0(x)))));
    mpz_mul(a, a, field(x)->base->base->base->order);
    mpz_add(a, a, rep(rep0(rep1(rep0(x)))));
    mpz_mul(a, a, field(x)->base->base->base->order);
    mpz_add(a, a, rep(rep0(rep0(rep0(x)))));
}

void bn254_fp12_to_oct(unsigned char *os, size_t *size, const Element x)
{
    size_t s0;

    unsigned char b0[380];
    mpz_t z;

    mpz_init(z);

    bn254_fp12_to_mpz(z, x);
    mpz_export(b0, &s0, 1, sizeof(*b0), 1, 0, z);

    memset(os, 0x00, 380);

    memcpy(&os[380 - (int)s0], b0, s0);

    (*size) = 380;

    mpz_clear(z);
}

void bn254_fp12_from_oct(Element x, const unsigned char *os, const size_t size)
{
    mpz_t quo, rem;

    if (size < 380) {
        fprintf(stderr, "error: please set up the enought buffer for element\n");
        exit(300);
    }

    mpz_init(quo);
    mpz_init(rem);

    mpz_import(quo, size, 1, sizeof(*os), 1, 0, os);

    mpz_tdiv_qr(quo, rem, quo, field(x)->base->base->base->order);
    mpz_set(rep(rep0(rep0(rep0(x)))), rem);
    mpz_tdiv_qr(quo, rem, quo, field(x)->base->base->base->order);
    mpz_set(rep(rep0(rep1(rep0(x)))), rem);
    mpz_tdiv_qr(quo, rem, quo, field(x)->base->base->base->order);
    mpz_set(rep(rep0(rep2(rep0(x)))), rem);
    mpz_tdiv_qr(quo, rem, quo, field(x)->base->base->base->order);
    mpz_set(rep(rep0(rep0(rep1(x)))), rem);
    mpz_tdiv_qr(quo, rem, quo, field(x)->base->base->base->order);
    mpz_set(rep(rep0(rep1(rep1(x)))), rem);
    mpz_tdiv_qr(quo, rem, quo, field(x)->base->base->base->order);
    mpz_set(rep(rep0(rep2(rep1(x)))), rem);
    mpz_tdiv_qr(quo, rem, quo, field(x)->base->base->base->order);
    mpz_set(rep(rep1(rep0(rep0(x)))), rem);
    mpz_tdiv_qr(quo, rem, quo, field(x)->base->base->base->order);
    mpz_set(rep(rep1(rep1(rep0(x)))), rem);
    mpz_tdiv_qr(quo, rem, quo, field(x)->base->base->base->order);
    mpz_set(rep(rep1(rep2(rep0(x)))), rem);
    mpz_tdiv_qr(quo, rem, quo, field(x)->base->base->base->order);
    mpz_set(rep(rep1(rep0(rep1(x)))), rem);
    mpz_tdiv_qr(quo, rem, quo, field(x)->base->base->base->order);
    mpz_set(rep(rep1(rep1(rep1(x)))), rem);
    mpz_tdiv_qr(quo, rem, quo, field(x)->base->base->base->order);
    mpz_set(rep(rep1(rep2(rep1(x)))), rem);

    mpz_clear(quo);
    mpz_clear(rem);
}
