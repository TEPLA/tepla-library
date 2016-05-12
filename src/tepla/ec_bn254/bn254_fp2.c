//==============================================================
//  extension field ( bn254_fp2 ) implementation with GMP
//--------------------------------------------------------------
//  bn254_fp2 Fp2 := Fp[x]/(x^2 + 5)
//--------------------------------------------------------------
//  2012.09.23 created by ishii
//==============================================================

#include "ec_bn254_lcl.h"

#define rep0(x) (((Element *)x->data)[0])
#define rep1(x) (((Element *)x->data)[1])

#define field(x) (x->field)
#define order(x) (x->field->order)
#define beta(x) ((x->field->irre_poly)[0])

#define MODOPT

//-------------------------------------------
//  initialization, clear, set
//-------------------------------------------
void bn254_fp2_init(Element x)
{
	x->data = (void *)malloc(sizeof(Element)*2);

	if( x->data == NULL ) { fprintf(stderr, "fail: allocate in fp2 init\n"); exit(100); }

	element_init(rep0(x), field(x)->base);
	element_init(rep1(x), field(x)->base);
}

void bn254_fp2_clear(Element x)
{
	if( x->data != NULL )
	{
		element_clear(rep0(x));
		element_clear(rep1(x));

		free(x->data); x->data = NULL;
	}
}

void bn254_fp2_set(Element x, const Element y)
{
	bn254_fp_set(rep0(x), rep0(y));
	bn254_fp_set(rep1(x), rep1(y));
}

void bn254_fp2_set_fp(Element z, const Element x, const Element y)
{
	bn254_fp_set(rep0(z), x);
	bn254_fp_set(rep1(z), y);
}

void bn254_fp2_set_str(Element x, const char *s)
{
	int i=0;

	char msg[140], *p, *c=NULL;

	strcpy(msg, s);

	p = msg;

	while( (*p) != '\0' ){ if((*p)==' '){ if(i==0){ c=p; } i++; } p++; }

	if( i != 1 ){ fprintf(stderr,"error: input string is not correct\n"); exit(200); }

	(*c) = '\0';

	bn254_fp_set_str(rep0(x), msg);
	bn254_fp_set_str(rep1(x), ++c);
}

void bn254_fp2_get_str(char *s, const Element x)
{
	char s1[65], s2[65];

	bn254_fp_get_str(s1, rep0(x));
	bn254_fp_get_str(s2, rep1(x));

	sprintf(s, "%s %s", s1, s2);
}

void bn254_fp2_set_zero(Element x)
{
	bn254_fp_set_zero(rep0(x));
	bn254_fp_set_zero(rep1(x));
}

void bn254_fp2_set_one(Element x)
{
	bn254_fp_set_one(rep0(x));
	bn254_fp_set_zero(rep1(x));
}

//-------------------------------------------
//  arithmetic operation
//-------------------------------------------
void bn254_fp2_add(Element z, const Element x, const Element y)
{
	bn254_fp_add(rep0(z), rep0(x), rep0(y));
	bn254_fp_add(rep1(z), rep1(x), rep1(y));
}

void bn254_fp2_neg(Element z, const Element x)
{
	bn254_fp_neg(rep0(z), rep0(x));
	bn254_fp_neg(rep1(z), rep1(x));
}

void bn254_fp2_sub(Element z, const Element x, const Element y)
{
	bn254_fp_sub(rep0(z), rep0(x), rep0(y));
	bn254_fp_sub(rep1(z), rep1(x), rep1(y));
}

//------------------------------------------------------
//  multipulication is implemented by Karatsuba method.
//------------------------------------------------------
void bn254_fp2_mul(Element z, const Element x, const Element y)
{
	Element* t = field(z)->base->tmp;

	bn254_fp_addn(t[1], rep0(x), rep1(x)); // t1 = x0 + x1
	bn254_fp_addn(t[2], rep0(y), rep1(y)); // t2 = y0 + y1
	bn254_fp_muln(t[0], t[1], t[2]);       // t0 = t1 * t2
	bn254_fp_muln(t[1], rep0(x), rep0(y)); // t1 = x0 * y0
	bn254_fp_muln(t[2], rep1(x), rep1(y)); // t2 = x1 * y1
	bn254_fp_subn(t[0], t[0], t[1]);       //
	bn254_fp_subn(rep1(z), t[0], t[2]);    // z1 = x0*y1 + y0*x1
	bn254_fp_mod(rep1(z), rep1(z));        //
	bn254_fp_addn(rep0(z), t[2], t[2]);    //
	bn254_fp_addn(rep0(z), rep0(z), rep0(z));
	bn254_fp_addn(rep0(z), rep0(z), t[2]); //
	bn254_fp_subn(rep0(z), t[1], rep0(z)); //
	bn254_fp_mod(rep0(z), rep0(z));        // z0 = t1 - 5*t2
}

//--------------------------------------------------------
//  multiplication of element of fp and element of fp^2
//--------------------------------------------------------
void bn254_fp2_mul_p(Element z, const Element x, const Element y)
{
	if ( field(x)->ID == bn254_fp )
	{
		bn254_fp_mul(rep0(z), x, rep0(y));
		bn254_fp_mul(rep1(z), x, rep1(y));
	}
	else if( field(y)->ID == bn254_fp )
	{
		bn254_fp_mul(rep0(z), rep0(x), y);
		bn254_fp_mul(rep1(z), rep1(x), y);
	}
	else
	{
		fprintf(stderr,"error: input should be element in bn254_fp2\n");
		exit(200);
	}
}

void bn254_fp2_inv(Element z, const Element x)
{
	Element* t = field(z)->base->tmp;

	bn254_fp_muln(t[1], rep1(x), rep1(x)); // t1 = a1^2
	bn254_fp_addn(t[0], t[1], t[1]);
	bn254_fp_addn(t[0], t[0], t[0]);
	bn254_fp_addn(t[1], t[1], t[0]);      // t1 = 5*a1^2
	bn254_fp_muln(t[0], rep0(x), rep0(x));// t0 = a0^2
	bn254_fp_addn(t[0], t[0], t[1]);      // t0 = t0 - t1
	bn254_fp_inv(t[1], t[0]);             // t1 = t0^-1
	bn254_fp_mul(rep0(z), rep0(x), t[1]); // c0 = a0*t1
	bn254_fp_mul(rep1(z), rep1(x), t[1]); // c1 = a1*t1
	bn254_fp_neg(rep1(z), rep1(z));       // c1 = -1*a1*t1
}

void bn254_fp2_dob(Element z, const Element x)
{
	bn254_fp_dob(rep0(z), rep0(x));
	bn254_fp_dob(rep1(z), rep1(x));
}

void bn254_fp2_tri(Element z, const Element x)
{
	bn254_fp_tri(rep0(z), rep0(x));
	bn254_fp_tri(rep1(z), rep1(x));
}

void bn254_fp2_xi_mul(Element z, const Element x)
{
	Element* t = field(z)->base->tmp;

	bn254_fp_add(t[0], rep1(x), rep1(x));
	bn254_fp_add(t[0], t[0], t[0]);
	bn254_fp_add(t[0], t[0], rep1(x));
	bn254_fp_set(rep1(z), rep0(x));
	bn254_fp_neg(rep0(z), t[0]);
}

void bn254_fp2_sqr(Element z, const Element x)
{
	Element* t = field(z)->base->tmp;

	bn254_fp_addn(t[0], rep1(x), rep1(x)); //
	bn254_fp_muln(t[0], t[0], rep0(x));    // t0 = 2*x1*x0
	bn254_fp_addp(t[1], rep0(x));
	bn254_fp_subn(t[1], t[1], rep1(x));    // t1 = x0-x1
	bn254_fp_addn(t[2], rep1(x), rep1(x)); //
	bn254_fp_addn(t[2], t[2], t[2]);       //
	bn254_fp_addn(t[2], t[2], rep1(x));    //
	bn254_fp_addn(t[2], t[2], rep0(x));    // t2 = 5*x1 + x0
	bn254_fp_muln(t[1], t[1], t[2]);       // t1 = t1 * t2
	bn254_fp_mod(rep1(z), t[0]);           // c1 = t0
	bn254_fp_addn(t[0], t[0], t[0]);       //
	bn254_fp_subn(t[1], t[1], t[0]);       // t1 = 2*t0*t1
	bn254_fp_mod(rep0(z), t[1]);           // c0 = t1
}

void bn254_fp2_pow(Element z, const Element x, const mpz_t exp)
{
	long t, i;
	Element c;

	element_init(c, field(z));
	element_set(c, x);

	t = (long)mpz_sizeinbase(exp, 2);

	for (i=t-2; i>=0; i--)
	{
		element_sqr(c, c);
		if( mpz_tstbit(exp, i) ){ element_mul(c, c, x); }
	}

	element_set(z, c);
	element_clear(c);
}

void bn254_fp2_frob_p(Element z, const Element x)
{
	bn254_fp_set(rep0(z), rep0(x));
	bn254_fp_neg(rep1(z), rep1(x));
}

void bn254_fp2_conj(Element z, const Element x)
{
	bn254_fp_set(rep0(z), rep0(x));
	bn254_fp_neg(rep1(z), rep1(x));
}

//---------------------------------------------------------
//  precomputation for sqrt
//---------------------------------------------------------
void bn254_fp2_precomp_sqrt(field_precomp_sqrt_p ps, const Field f)
{
	//-----------------------------
	//  decompose of value
	//    (p^2-1) = 2^e * v
	//-----------------------------
	mpz_init_set(ps->v, f->order);
	mpz_sub_ui(ps->v, ps->v, 1);
	ps->e = (int)mpz_scan1(ps->v, 0);
	mpz_fdiv_q_2exp(ps->v, ps->v, ps->e);

	//----------------------------------
	//  n_v = n^v :
	//    n is some integer (n/p) = -1
	//----------------------------------
	element_init(ps->n_v, f);
	do { element_random(ps->n_v); }
	while( element_is_sqr(ps->n_v) );
	element_pow(ps->n_v, ps->n_v, ps->v);
}

//---------------------------------------------------------
//  precomputation for Fp2 operation
//---------------------------------------------------------
void bn254_fp2_precomp(Field f)
{
	field_precomp_p precomp=NULL;

	precomp = (field_precomp_p)malloc(sizeof(struct ec_field_precomp_st));

	precomp->ps = (field_precomp_sqrt_p)malloc(sizeof(struct ec_field_precomp_sqrt_st));
	bn254_fp2_precomp_sqrt(precomp->ps, f);

	precomp->pf = NULL;

	f->precomp = (void *)precomp;
}

//--------------------------------------------------
//  square root in extended Fp
//--------------------------------------------------
int bn254_fp2_sqrt(Element z, const Element x)
{
	mpz_t _v;
	int m, r, i;
	
	field_precomp_sqrt_p ps;

	Element *t = field(z)->tmp;

	if( !element_is_sqr(x) ){ return FALSE; }

	ps = ((field_precomp_p)(field(x)->precomp))->ps;

	element_set(t[0], ps->n_v);  // t0 = n^v

	r = ps->e; // r = e

	mpz_init_set(_v, ps->v);
	mpz_sub_ui(_v, _v, 1);
	mpz_tdiv_q_2exp(_v, _v, 1);

	element_pow(t[1], x, _v);   // t1 = x^{(v-1)/2}
	element_sqr(t[2], t[1]);
	element_mul(t[2], t[2], x); // t2 = x*t1^2
	element_mul(t[1], x, t[1]); // t1 = x*t1

	mpz_clear(_v);

	while ( !element_is_one(t[2]) )
	{
		m = 0;
		element_set(t[3], t[2]);

		do { element_sqr(t[3], t[3]); m++; }
		while ( !element_is_one(t[3]) && m < r );

		r = r-m-1;
		element_set(t[3], t[0]);
		for (i=1; i<=r; i++) { element_sqr(t[3], t[3]); } // t3 = t2^{r-m-1}
		element_sqr(t[0], t[3]);      // t0 = t3^2
		r = m;
		element_mul(t[1], t[1], t[3]);// t1 = t1*t3
		element_mul(t[2], t[2], t[0]);// t2 = t2*t0
	}

	element_set(z, t[1]);

	return TRUE;
}

//-------------------------------------------
//  comparison operation
//-------------------------------------------
int bn254_fp2_is_zero(const Element x)
{
	return ( bn254_fp_is_zero(rep1(x)) && bn254_fp_is_zero(rep0(x)) );
}

int bn254_fp2_is_one(const Element x)
{
	if ( bn254_fp_is_zero(rep1(x)) )
	{
		return ( bn254_fp_is_one(rep0(x)) );
	}
	return FALSE;
}

int bn254_fp2_cmp(const Element x, const Element y)
{
	if ( bn254_fp_cmp(rep1(x), rep1(y)) == 0 )
	{
		if( bn254_fp_cmp(rep1(x), rep1(y)) == 0 ){ return 0; }
	}
	return 1;
}

int bn254_fp2_is_sqr(const Element x)
{
	int hr = FALSE;

	Element *t = field(x)->base->tmp;

	if( element_is_zero(x) ){ return FALSE; }

	bn254_fp_inv(t[0], rep1(x));
	bn254_fp_mul(t[0], t[0], rep0(x));
	bn254_fp_sqr(t[0], t[0]);
	bn254_fp_add(t[0], t[0], field(x)->irre_poly[0]);

	hr = bn254_fp_is_sqr(t[0]);

	return hr;
}

//-------------------------------------------
//  generate random element
//-------------------------------------------
void bn254_fp2_random(Element z)
{
	element_random(rep0(z));
	element_random(rep1(z));
}

//-------------------------------------------
//  i/o operation (octet string)
//-------------------------------------------
void bn254_fp2_to_oct(unsigned char *os, size_t *size, const Element x)
{
	size_t s0, s1;

	unsigned char os0[32];
	unsigned char os1[32];

	bn254_fp_to_oct(os0, &s0, rep0(x));
	bn254_fp_to_oct(os1, &s1, rep1(x));

	memset(os, 0x00, 64);

	memcpy(&(os[0]), os0, s0);
	memcpy(&(os[32]), os1, s1);

	(*size) = 64;
}

void bn254_fp2_from_oct(Element x, const unsigned char *os, const size_t size)
{
	if( size < 64 ){ fprintf(stderr, "error: please set up the enought buffer for element\n"); exit(300); }

	bn254_fp_from_oct(rep0(x), &(os[0]), 32);
	bn254_fp_from_oct(rep1(x), &(os[32]), 32);
}
