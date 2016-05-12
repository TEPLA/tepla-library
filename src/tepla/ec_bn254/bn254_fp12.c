//==============================================================
//  prime field ( bn254_fp12 ) implementation with GMP
//--------------------------------------------------------------
//  bn254_fp12 Fp12 := Fp6[z]/(z^2 - y)
//--------------------------------------------------------------
//  2012.09.23 created by ishii
//==============================================================

#include "ec_bn254_lcl.h"

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
	x->data = (void *)malloc(sizeof(Element)*2);

	if( x->data == NULL ){ fprintf(stderr, "fail: allocate in bn254_fp12 init\n"); exit(100); }

	element_init(rep0(x), field(x)->base);
	element_init(rep1(x), field(x)->base);
}

void bn254_fp12_clear(Element x)
{
	if( x->data != NULL )
	{
		element_clear(rep0(x));
		element_clear(rep1(x));

		free(x->data); x->data = NULL;
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
	int  i=0;
	char msg[780], *p, *c[11];

	strcpy(msg, s);

	p = msg;

	while( (*p) != '\0' )
	{
		if( (*p)==' ' ){ if( i<11 ){ c[i]=p; } i++; }
		p++;
	}

	if( i != 11 ){ fprintf(stderr,"error: input string is not correct %d\n", i); exit(200); }

	(*c[5]) = '\0';

	bn254_fp6_set_str(rep0(x), msg);
	bn254_fp6_set_str(rep1(x), ++c[5]);
}

void bn254_fp12_get_str(char *s, const Element x)
{
	char s1[390], s2[390];

	bn254_fp6_get_str(s1, rep0(x));
	bn254_fp6_get_str(s2, rep1(x));

	sprintf(s, "%s %s", s1, s2);
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

	for (i=t-2; i>=0; i--)
	{
		element_sqr(c, c);
		if( mpz_tstbit(exp, i) ){ element_mul(c, c, x); }
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

	naf = (int *)malloc(sizeof(int)*(t+1));

	generate_naf(naf, &nlen, exp);

	for (i=nlen-2; i>=0; i--)
	{
		element_sqr(c, c);
		if( naf[i] )
		{
			if( naf[i] < 0 ){ element_mul(c, c, ix); }
			else { element_mul(c, c, x); }
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

	bn254_fp2_mul_p(rep0(rep1(z)), rep0(rep1(z)), (pf->gamma1)[0]);   //t2 = t2*gamma1
	bn254_fp2_mul_p(rep1(rep0(z)), rep1(rep0(z)), (pf->gamma1)[1]);   //t3 = t3*gamma2
	bn254_fp2_mul_p(rep1(rep1(z)), rep1(rep1(z)), (pf->gamma1)[2]);   //t4 = t4*gamma3
	bn254_fp2_mul_p(rep2(rep0(z)), rep2(rep0(z)), (pf->gamma1)[3]);   //t5 = t5*gamma4
	bn254_fp2_mul_p(rep2(rep1(z)), rep2(rep1(z)), (pf->gamma1)[4]);   //t6 = t6*gamma5
}

void bn254_fp12_frob_p2(Element z, const Element x)
{
	field_precomp_frob_p pf;
	
	pf = ((field_precomp_p)(field(z)->precomp))->pf;

	bn254_fp2_set(rep0(rep0(z)), rep0(rep0(x)));   //t1 = g0;
	bn254_fp2_mul_p(rep0(rep1(z)), rep0(rep1(x)), (pf->gamma2)[0]);   //t2 = h0*gamma1
	bn254_fp2_mul_p(rep1(rep0(z)), rep1(rep0(x)), (pf->gamma2)[1]);   //t3 = g1*gamma2
	bn254_fp2_mul_p(rep1(rep1(z)), rep1(rep1(x)), (pf->gamma2)[2]);   //t4 = h1*gamma3
	bn254_fp2_mul_p(rep2(rep0(z)), rep2(rep0(x)), (pf->gamma2)[3]);   //t5 = g2*gamma4
	bn254_fp2_mul_p(rep2(rep1(z)), rep2(rep1(x)), (pf->gamma2)[4]);   //t6 = h2*gamma5
}

void bn254_fp12_frob_p3(Element z, const Element x)
{
	field_precomp_frob_p pf;
	
	pf = ((field_precomp_p)(field(z)->precomp))->pf;

	bn254_fp6_conj(rep0(z), rep0(x));
	bn254_fp6_conj(rep1(z), rep1(x));

	bn254_fp2_mul_p(rep0(rep1(z)), rep0(rep1(z)), (pf->gamma3)[0]);   //t2 = t2*gamma1
	bn254_fp2_mul_p(rep1(rep0(z)), rep1(rep0(z)), (pf->gamma3)[1]);   //t3 = t3*gamma2
	bn254_fp2_mul_p(rep1(rep1(z)), rep1(rep1(z)), (pf->gamma3)[2]);   //t4 = t4*gamma3
	bn254_fp2_mul_p(rep2(rep0(z)), rep2(rep0(z)), (pf->gamma3)[3]);   //t5 = t5*gamma4
	bn254_fp2_mul_p(rep2(rep1(z)), rep2(rep1(z)), (pf->gamma3)[4]);   //t6 = t6*gamma5
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

void bn254_fp12_sqr_forpairing(Element z, const Element x)
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

//------------------------------------------------------------------
//  Special exponentiation in Fp12
//    input  : element x
//    output : z = x^t : t = 2^62-2^54+2^44
//-------------------------------------------------------------------
void bn254_fp12_pow_forpairing(Element z, const Element x, const int *t, int tlen)
{
	int i;
	Element ix;

	if( z == x )
	{
		fprintf(stderr, "error: bad input for bn254_fp12_pow_forpairing\n");
		exit(200);
	}

	element_init(ix, field(x));

	bn254_fp12_set(z, x);
	bn254_fp12_conj(ix, x);

	for (i=tlen-2; i>=0; i--)
	{
		bn254_fp12_sqr_forpairing(z, z);

		if ( t[i] )
		{
			if ( t[i] < 0 ) { bn254_fp12_mul(z, z, ix); }
			else { bn254_fp12_mul(z, z, x); }
		}
	}

	element_clear(ix);
}

//---------------------------------------------------------
//  precomputation for Fp12 frobenius
//---------------------------------------------------------
void bn254_fp12_precomp_frob(field_precomp_frob_p pf, const Field f)
{
	int i;
	mpz_t  t, u;
	Element xi, tmp;

	//---------------------------------
	//  allocate buffer
	//---------------------------------
	Element *g1 = (Element *)malloc(sizeof(Element)*5);
	Element *g2 = (Element *)malloc(sizeof(Element)*5);
	Element *g3 = (Element *)malloc(sizeof(Element)*5);

	struct ec_field_st *fp  = f->base->base->base;
	struct ec_field_st *fp2 = f->base->base;

	//---------------------------------
	//  initialization
	//---------------------------------
	for(i=0;i<5;i++)
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

	for(i=1;i<5;i++) { element_mul(g1[i], g1[i-1], g1[0]); }
	for(i=0;i<5;i++) { element_mul(g2[i], g1[i], g1[i]); }
	for(i=0;i<5;i++) { element_mul(g3[i], g1[i], g2[i]); }

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
	field_precomp_p precomp=NULL;

	precomp = (field_precomp_p)malloc(sizeof(struct ec_field_precomp_st));

	precomp->ps = (field_precomp_sqrt_p)malloc(sizeof(struct ec_field_precomp_sqrt_st));
	bn254_fp2_precomp_sqrt(precomp->ps, f);

	precomp->pf = (field_precomp_frob_p)malloc(sizeof(struct ec_field_precomp_frob_st));
	bn254_fp12_precomp_frob(precomp->pf, f);

	f->precomp = (void *)precomp;
}

//-------------------------------------------
//  comparison operation
//-------------------------------------------
int bn254_fp12_is_zero(const Element x)
{
	if ( bn254_fp6_is_zero(rep1(x)) )
	{
		if( bn254_fp6_is_zero(rep0(x)) ){ return TRUE; }
	}
	return FALSE;
}

int bn254_fp12_is_one(const Element x)
{
	if ( bn254_fp6_is_zero(rep1(x)) )
	{
		if( bn254_fp6_is_one(rep0(x)) ){ return TRUE; }
	}
	return FALSE;
}

int bn254_fp12_cmp(const Element x, const Element y)
{
	if ( bn254_fp6_cmp(rep1(x), rep1(y)) == 0 )
	{
		if( bn254_fp6_cmp(rep0(x), rep0(y)) == 0 ){ return 0; }
	}
	return 1;
}

int bn254_fp12_is_sqr(const Element x)
{
	int hr = FALSE;

	Element *t = field(x)->base->tmp;

	if( element_is_zero(x) ){ return FALSE; }

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
void bn254_fp12_to_oct(unsigned char *os, size_t *size, const Element x)
{
	size_t s0, s1;

	unsigned char b0[192];
	unsigned char b1[192];

	bn254_fp6_to_oct(b0, &s0, rep0(x));
	bn254_fp6_to_oct(b1, &s1, rep1(x));

	memset(os, 0x00, 384);

	memcpy(&(os[0]),   b0, s0);
	memcpy(&(os[192]), b1, s1);

	(*size) = 384;
}

void bn254_fp12_from_oct(Element x, const unsigned char *os, const size_t size)
{
	if( size < 384 ){ fprintf(stderr, "error: please set up the enough buffer for element\n"); exit(300); }

	bn254_fp6_from_oct(rep0(x), &(os[0]),   192);
	bn254_fp6_from_oct(rep1(x), &(os[192]), 192);
}
