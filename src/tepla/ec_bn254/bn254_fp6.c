//==============================================================
//  extension field ( bn254_fp6 ) implementation with GMP
//--------------------------------------------------------------
//  bn254_fp6 Fp6 := Fp2[y]/(y^3 - x)
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
void bn254_fp6_init(Element x)
{
	x->data = (void *)malloc(sizeof(Element)*3);

	if( x->data == NULL ){ fprintf(stderr, "fail: allocate in bn254_fp6 init\n"); exit(100); }

	element_init(rep0(x), field(x)->base);
	element_init(rep1(x), field(x)->base);
	element_init(rep2(x), field(x)->base);
}

void bn254_fp6_clear(Element x)
{
	if( x->data != NULL )
	{
		element_clear(rep0(x));
		element_clear(rep1(x));
		element_clear(rep2(x));

		free(x->data); x->data = NULL;
	}
}

void bn254_fp6_set(Element x, const Element y)
{
	bn254_fp2_set(rep0(x), rep0(y));
	bn254_fp2_set(rep1(x), rep1(y));
	bn254_fp2_set(rep2(x), rep2(y));
}

void bn254_fp6_set_fp2(Element z, const Element w, const Element x, const Element y)
{
	bn254_fp2_set(rep0(z), w);
	bn254_fp2_set(rep1(z), x);
	bn254_fp2_set(rep2(z), y);
}

void bn254_fp6_set_str(Element x, const char *s)
{
	int  i=0;

	char msg[400], *p, *c[5];

	strcpy(msg, s);

	p = msg;

	while( (*p) != '\0' )
	{
		if( (*p)==' ' ){ if( i<5 ){ c[i]=p; } i++; }
		p++;
	}

	if( i != 5 ) { fprintf(stderr,"error: input string is not correct\n"); exit(200); }

	(*c[1]) = '\0';
	(*c[3]) = '\0';

	bn254_fp2_set_str(rep0(x), msg);
	bn254_fp2_set_str(rep1(x), ++c[1]);
	bn254_fp2_set_str(rep2(x), ++c[3]);
}

void bn254_fp6_get_str(char *s, const Element x)
{
	char s0[130], s1[130], s2[130];

	bn254_fp2_get_str(s0, rep0(x));
	bn254_fp2_get_str(s1, rep1(x));
	bn254_fp2_get_str(s2, rep2(x));

	sprintf(s, "%s %s %s", s0, s1, s2);
}

void bn254_fp6_set_zero(Element x)
{
	bn254_fp2_set_zero(rep0(x));
	bn254_fp2_set_zero(rep1(x));
	bn254_fp2_set_zero(rep2(x));
}

void bn254_fp6_set_one(Element x)
{
	bn254_fp2_set_one(rep0(x));
	bn254_fp2_set_zero(rep1(x));
	bn254_fp2_set_zero(rep2(x));
}

//-------------------------------------------
//  arithmetic operation
//-------------------------------------------
void bn254_fp6_add(Element z, const Element x, const Element y)
{
	bn254_fp2_add(rep0(z), rep0(x), rep0(y));
	bn254_fp2_add(rep1(z), rep1(x), rep1(y));
	bn254_fp2_add(rep2(z), rep2(x), rep2(y));
}

void bn254_fp6_neg(Element z, const Element x)
{
	bn254_fp2_neg(rep0(z), rep0(x));
	bn254_fp2_neg(rep1(z), rep1(x));
	bn254_fp2_neg(rep2(z), rep2(x));
}

void bn254_fp6_sub(Element z, const Element x, const Element y)
{
	bn254_fp2_sub(rep0(z), rep0(x), rep0(y));
	bn254_fp2_sub(rep1(z), rep1(x), rep1(y));
	bn254_fp2_sub(rep2(z), rep2(x), rep2(y));
}

void bn254_fp6_mul(Element z, const Element x, const Element y)
{
	Element *t = field(z)->base->tmp;

	bn254_fp2_mul(t[0], rep0(x), rep0(y));  // t0 = a0*b0
	bn254_fp2_mul(t[1], rep1(x), rep1(y));  // t1 = a1*b1
	bn254_fp2_mul(t[2], rep2(x), rep2(y));  // t2 = a2*b2

	//------------------------------------------
	//  c0 = ((a1+a2)*(b1+b2)-t1-t2)*xi + t0
	//------------------------------------------
	bn254_fp2_add(t[4], rep1(x), rep2(x));  // c0 = a1+a2
	bn254_fp2_add(t[3], rep1(y), rep2(y));  // t3 = b1+b2
	bn254_fp2_mul(t[4], t[4], t[3]);        // c0 = c0*t3
	bn254_fp2_sub(t[4], t[4], t[1]);        //
	bn254_fp2_sub(t[4], t[4], t[2]);        // c0 = c0-t1-t2
	bn254_fp2_xi_mul(t[4], t[4]);           // c0 = c0*xi
	bn254_fp2_add(t[4], t[4], t[0]);        // c0 = c0 + t0

	//------------------------------------------
	//  c1 = (a0+a1)*(b0+b1) - t0 - t1 + xi*t2
	//------------------------------------------
	bn254_fp2_add(t[3], rep0(x), rep1(x));    // t3 = a0+a1
	bn254_fp2_add(rep1(z), rep0(y), rep1(y)); // c1 = b0+b1
	bn254_fp2_mul(t[3], t[3], rep1(z));       // t3 = t3*c1
	bn254_fp2_sub(t[3], t[3], t[0]);          //
	bn254_fp2_sub(t[3], t[3], t[1]);          // t3 = t3 - t0 - t1
	bn254_fp2_xi_mul(rep1(z), t[2]);          // c1 = xi*t2
	bn254_fp2_add(rep1(z), rep1(z), t[3]);    // c0 = c0 + t3

	//------------------------------------------
	//  c2 = (a0+a2)*(b0+b2) - t0 - t2 + t1
	//------------------------------------------
	bn254_fp2_add(t[3], rep0(x), rep2(x));    // t3 = a0+a2
	bn254_fp2_add(rep2(z), rep0(y), rep2(y)); // c2 = b0+b2
	bn254_fp2_mul(rep2(z), rep2(z), t[3]);    // c2 = c2*t3
	bn254_fp2_sub(rep2(z), rep2(z), t[0]);    //
	bn254_fp2_sub(rep2(z), rep2(z), t[2]);    // c2 = c2 - t0 - t2
	bn254_fp2_add(rep2(z), rep2(z), t[1]);    // c2 = c2 + t1

	bn254_fp2_set(rep0(z), t[4]);
}

//--------------------------------------------------------
//   z = x * gamma ( Fp12 : Fp6[x]/x^2-gamma )
//--------------------------------------------------------
void bn254_fp6_gm_mul(Element z, const Element x)
{
	if( z == x ){ fprintf(stderr, "fail gm mul\n"); exit(500); }

	bn254_fp2_xi_mul(rep0(z), rep2(x));
	bn254_fp2_set(rep1(z), rep0(x));
	bn254_fp2_set(rep2(z), rep1(x));
}

//--------------------------------------------------------
//   z = x * (y, 0, 0)
//--------------------------------------------------------
void bn254_fp6_mul_fp2(Element z, const Element x, const Element y)
{
	if( field(y)->ID != bn254_fp2 )
	{
		fprintf(stderr, "error: input should be element in bn254_fp6\n");
		exit(200);
	}

	bn254_fp2_mul(rep0(z), rep0(x), y);
	bn254_fp2_mul(rep1(z), rep1(x), y);
	bn254_fp2_mul(rep2(z), rep2(x), y);
}

//--------------------------------------------------------
//   z = x * (y1, y2, 0)
//--------------------------------------------------------
void bn254_fp6_mul_fp2_2(Element z, const Element x, const Element y1, const Element y2)
{
	Element *v = field(z)->base->tmp;

	if( field(y1)->ID != bn254_fp2 || field(y2)->ID != bn254_fp2 )
	{
		fprintf(stderr, "error: input should be element in bn254_fp6\n");
		exit(200);
	}

	bn254_fp2_mul(v[0], rep0(x), y1);
	bn254_fp2_mul(v[1], rep1(x), y2);
	bn254_fp2_add(v[2], rep1(x), rep2(x));
	bn254_fp2_mul(v[2], v[2], y2);
	bn254_fp2_sub(v[2], v[2], v[1]);
	bn254_fp2_xi_mul(v[2], v[2]);
	bn254_fp2_add(v[2], v[2], v[0]);
	bn254_fp2_add(rep1(z), rep0(x), rep1(x));
	bn254_fp2_add(v[3], y1, y2);
	bn254_fp2_mul(rep1(z), rep1(z), v[3]);
	bn254_fp2_sub(rep1(z), rep1(z), v[0]);
	bn254_fp2_sub(rep1(z), rep1(z), v[1]);
	bn254_fp2_add(rep2(z), rep0(x), rep2(x));
	bn254_fp2_mul(rep2(z), rep2(z), y1);
	bn254_fp2_sub(rep2(z), rep2(z), v[0]);
	bn254_fp2_add(rep2(z), rep2(z), v[1]);
	bn254_fp2_set(rep0(z), v[2]);
}

void bn254_fp6_inv(Element z, const Element x)
{
	Element *t = field(z)->base->tmp;

	bn254_fp2_sqr(t[0], rep0(x));   // t0 = a0^2
	bn254_fp2_sqr(t[1], rep1(x));   // t1 = a1^2
	bn254_fp2_sqr(t[2], rep2(x));   // t2 = a2^2
	bn254_fp2_mul(t[3], rep0(x), rep1(x));   // t3 = a0*a1
	bn254_fp2_mul(t[4], rep0(x), rep2(x));   // t4 = a0*a2
	bn254_fp2_mul(t[5], rep1(x), rep2(x));   // t5 = a1*a2

	//-------------------------
	// c0 = t0 - xi*t5
	//-------------------------
	bn254_fp2_xi_mul(t[5], t[5]);     // c0 = xi*t5
	bn254_fp2_sub(t[0], t[0], t[5]);  // t0 = t0 - c0

	//-------------------------
	// c1 = xi*t2 - t3
	//-------------------------
	bn254_fp2_xi_mul(t[2], t[2]);     //
	bn254_fp2_sub(t[2], t[2], t[3]);  // t2 = xi*ts - t3

	//-------------------------
	// c2 = t1*t4
	//-------------------------
	bn254_fp2_sub(t[1], t[1], t[4]);  // t1 = t1-t4

	bn254_fp2_mul(t[4], rep0(x), t[0]); // t4 = a0*c0
	bn254_fp2_mul(t[3], rep2(x), t[2]); // t3 = a2*c1
	bn254_fp2_xi_mul(t[3], t[3]);       // t3 = t3*xi
	bn254_fp2_add(t[4], t[4], t[3]);    // t4 = t4+t3
	bn254_fp2_mul(t[3], rep1(x), t[1]); // t3 = a1*c2
	bn254_fp2_xi_mul(t[3], t[3]);       // t3 = t3*xi
	bn254_fp2_add(t[4], t[4], t[3]);    // t4 = t4+t3
	bn254_fp2_inv(t[4], t[4]);          // t4 = t4^-1

	bn254_fp2_mul(rep0(z), t[0], t[4]);
	bn254_fp2_mul(rep1(z), t[2], t[4]);
	bn254_fp2_mul(rep2(z), t[1], t[4]);
}

void bn254_fp6_dob(Element z, const Element x)
{
	bn254_fp2_dob(rep0(z), rep0(x));
	bn254_fp2_dob(rep1(z), rep1(x));
	bn254_fp2_dob(rep2(z), rep2(x));
}

void bn254_fp6_tri(Element z, const Element x)
{
	bn254_fp2_tri(rep0(z), rep0(x));
	bn254_fp2_tri(rep1(z), rep1(x));
	bn254_fp2_tri(rep2(z), rep2(x));
}

void bn254_fp6_sqr(Element z, const Element x)
{
	Element *t = field(z)->base->tmp;

	bn254_fp2_mul(t[0], rep0(x), rep1(x));
	bn254_fp2_add(t[0], t[0], t[0]);  // t0 = 2*a0*a1
	bn254_fp2_sqr(t[1], rep2(x));     // t1 = a2^2

	//-------------------------
	// c1 = t1*xi + t0
	//-------------------------
	bn254_fp2_xi_mul(t[2], t[1]);    //
	bn254_fp2_add(t[2], t[2], t[0]); // t2 = t1*xi + t0

	//-------------------------
	// c2 = t0 - t1
	//-------------------------
	bn254_fp2_sub(t[0], t[0], t[1]); // t0 = t0 - t1

	bn254_fp2_sqr(t[1], rep0(x));              // t1 = a0^2
	bn254_fp2_sub(rep0(z), rep0(x), rep1(x));  //
	bn254_fp2_add(rep0(z), rep0(z), rep2(x));  // v0 = a0 - a1 + a2
	bn254_fp2_mul(rep1(z), rep1(x), rep2(x));  //
	bn254_fp2_add(rep1(z), rep1(z), rep1(z));  // v1 = 2*a1*a2
	bn254_fp2_sqr(rep2(z), rep0(z));           // v2 = v0^2

	//-------------------------
	// c0 = v1*xi + t2
	//-------------------------
	bn254_fp2_xi_mul(rep0(z), rep1(z));    // c0 = v1*xi
	bn254_fp2_add(rep0(z), rep0(z), t[1]); // c0 = c0 + t1

	//-------------------------
	// c2 = c2 + t0 + t1 - t2
	//-------------------------
	bn254_fp2_add(rep2(z), rep2(z), t[0]);    // c2 = v2 + t0
	bn254_fp2_add(rep2(z), rep2(z), rep1(z)); // c2 = c2 + v1
	bn254_fp2_sub(rep2(z), rep2(z), t[1]);    // c2 = c2 - t2

	bn254_fp2_set(rep1(z), t[2]);
}

/*
void bn254_fp6_frob_p(Element z, const Element x)
{
}
*/

void bn254_fp6_conj(Element z, const Element x)
{
	bn254_fp2_conj(rep0(z), rep0(x));
	bn254_fp2_conj(rep1(z), rep1(x));
	bn254_fp2_conj(rep2(z), rep2(x));
}

//---------------------------------------------------------
//  precomputation for Fp6 operation
//---------------------------------------------------------
void bn254_fp6_precomp(Field f)
{
	field_precomp_p precomp=NULL;

	precomp = (field_precomp_p)malloc(sizeof(struct ec_field_precomp_st));

	precomp->ps = (field_precomp_sqrt_p)malloc(sizeof(struct ec_field_precomp_sqrt_st));
	bn254_fp2_precomp_sqrt(precomp->ps, f);

	precomp->pf = NULL;

	f->precomp = (void *)precomp;
}

//-------------------------------------------
//  comparison operation
//-------------------------------------------
int bn254_fp6_is_zero(const Element x)
{
	if ( bn254_fp2_is_zero(rep2(x)) )
	{
		if ( bn254_fp2_is_zero(rep1(x)) )
		{
			if ( bn254_fp2_is_zero(rep0(x)) ){ return TRUE; }
		}
	}
	return FALSE;
}

int bn254_fp6_is_one(const Element x)
{
	if ( bn254_fp2_is_zero(rep2(x)) )
	{
		if ( bn254_fp2_is_zero(rep1(x)) )
		{
			if( bn254_fp2_is_one(rep0(x)) ){ return TRUE; }
		}
	}
	return FALSE;
}

int bn254_fp6_cmp(const Element x, const Element y)
{
	if ( bn254_fp2_cmp(rep2(x), rep2(y)) == 0 )
	{
		if ( bn254_fp2_cmp(rep1(x), rep1(y)) == 0 )
		{
			if( bn254_fp2_cmp(rep0(x), rep0(y)) == 0 ){ return 0; }
		}
	}
	return 1;
}

int bn254_fp6_is_sqr(const Element x)
{
	int k = 1;

	Element *t = field(x)->base->tmp;

	if( element_is_zero(x) ){ return FALSE; }

	k *= bn254_fp2_is_sqr(rep2(x))? 1: -1;

	bn254_fp2_sqr(t[1], rep1(x));
	bn254_fp2_mul(t[2], rep0(x), rep2(x));
	bn254_fp2_sub(t[1], t[1], t[2]);      // t1 = x1^2-x0*x2
	bn254_fp2_mul(t[2], rep0(x), rep1(x));
	bn254_fp2_sqr(t[3], rep2(x));
	bn254_fp2_xi_mul(t[3], t[3]);
	bn254_fp2_sub(t[2], t[2], t[3]);      // t2 = x0*x1-x2^2*xi
	bn254_fp2_inv(t[1], t[1]);
	bn254_fp2_mul(t[1], t[1], t[2]);      // t1 = t2 / t1

	bn254_fp2_inv(t[2], rep2(x));
	bn254_fp2_mul(t[3], t[2], rep1(x));
	bn254_fp2_sub(t[3], t[3], t[1]);
	bn254_fp2_mul(t[3], t[3], t[1]);      // t3 = ((x1/x2)-t1)t1

	bn254_fp2_mul(t[2], t[2], rep0(x));
	bn254_fp2_sub(t[2], t[2], t[3]);      // t2 = (x0/x2)-t3

	k *= bn254_fp2_is_sqr(t[2])? 1: -1;

	return (k == 1);
}

//-------------------------------------------
//  generate random element
//-------------------------------------------
void bn254_fp6_random(Element z)
{
	bn254_fp2_random(rep0(z));
	bn254_fp2_random(rep1(z));
	bn254_fp2_random(rep2(z));
}

//-------------------------------------------
//  i/o operation (octet string)
//-------------------------------------------
void bn254_fp6_to_oct(unsigned char *os, size_t *size, const Element x)
{
	size_t s0, s1, s2;

	unsigned char b0[64];
	unsigned char b1[64];
	unsigned char b2[64];

	bn254_fp2_to_oct(b0, &s0, rep0(x));
	bn254_fp2_to_oct(b1, &s1, rep1(x));
	bn254_fp2_to_oct(b2, &s2, rep2(x));

	memset(os, 0x00, 192);

	memcpy(&(os[0]),   b0, s0);
	memcpy(&(os[64]),  b1, s1);
	memcpy(&(os[128]), b2, s2);

	(*size) = 192;
}

void bn254_fp6_from_oct(Element x, const unsigned char *os, const size_t size)
{
	if( size < 192 ){ fprintf(stderr, "error: please set up the enought buffer for element\n"); exit(300); }

	bn254_fp2_from_oct(rep0(x), &(os[0]),   64);
	bn254_fp2_from_oct(rep1(x), &(os[64]),  64);
	bn254_fp2_from_oct(rep2(x), &(os[128]), 64);
}
