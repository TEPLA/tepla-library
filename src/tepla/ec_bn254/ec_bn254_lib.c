//==================================================
// Function which making instance of structures
//==================================================
#include "ec_bn254_lcl.h"

#define TMP_NUM 10

#define SAFE_FREE(p) { if(p){ free(p); (p)=NULL; } }

//-----------------------------------------------------
// utility for set field/curve name
//-----------------------------------------------------
void set_field_name(Field f, const char* name)
{
	int len = strlen(name) + 1;

	f->field_name = (char*)malloc(sizeof(char)*len);

	strcpy(f->field_name, name);
}

void set_curve_name(EC_GROUP ec, const char* name)
{
	int len = strlen(name) + 1;

	ec->curve_name = (char*)malloc(sizeof(char)*len);

	strcpy(ec->curve_name, name);
}

//----------------------------------------------
//  function creating field fp254_p
//----------------------------------------------
void ec_bn254_fp_new(Field f)
{
	int i;

	f->type = Field_fp;

	set_field_name(f, "bn254_fp");

	f->ID = bn254_fp;

	f->str_len = 65;
	f->oct_len = 32;

	//------------------------------
	//  set pointer of function
	//------------------------------
	f->init     = bn254_fp_init;
	f->clear    = bn254_fp_clear;
	f->set      = bn254_fp_set;
	f->set_str  = bn254_fp_set_str;
	f->get_str  = bn254_fp_get_str;
	f->set_zero = bn254_fp_set_zero;
	f->set_one  = bn254_fp_set_one;

	f->add  = bn254_fp_add;
	f->sub  = bn254_fp_sub;
	f->neg  = bn254_fp_neg;
	f->mul  = bn254_fp_mul;
	f->sqr  = bn254_fp_sqr;
	f->inv  = bn254_fp_inv;
	f->pow  = bn254_fp_pow;
	f->sqrt = bn254_fp_sqrt;

	f->is_zero = bn254_fp_is_zero;
	f->is_one  = bn254_fp_is_one;
	f->is_sqr  = bn254_fp_is_sqr;
	f->cmp     = bn254_fp_cmp;

	f->random = bn254_fp_random;

	f->to_oct   = bn254_fp_to_oct;
	f->from_oct = bn254_fp_from_oct;

	//-----------------------------------------
	//  base field
	//-----------------------------------------
	f->base = NULL;

	//-----------------------------------------
	//  characteristic of prime field
	//-----------------------------------------
	//    p = 36t^4 + 36t^3 + 24t^2 + 6t + 1
	//    t = 2^62 - 2^54 + 2^44
	//-----------------------------------------
	mpz_init_set_str(f->order, "2370FB049D410FBE4E761A9886E502417D023F40180000017E80600000000001", 16);

	//-----------------------------------------
	//  irreducible polynomial
	//-----------------------------------------
	f->irre_poly_num = 0;
	f->irre_poly_deg = 1;
	f->irre_poly = NULL;

	//-----------------------------------------
	//  pre-computation for square root
	//-----------------------------------------
	f->precomp = NULL;

	//----------------------------------
	//  temporary element init
	//----------------------------------
	f->tmp = (Element *)malloc(sizeof(Element)*TMP_NUM);
	for(i=0;i<TMP_NUM;i++){ element_init(f->tmp[i], f); }

	return;
}

//----------------------------------------------
//  function creating field bn254_fp2
//----------------------------------------------
void ec_bn254_fp2_new(Field f)
{
	int i;

	f->type = Field_fpn;

	set_field_name(f, "bn254_fp2");

	f->ID = bn254_fp2;

	f->str_len = 130;
	f->oct_len = 64;

	//------------------------------
	//  set pointer of function
	//------------------------------
	f->init     = bn254_fp2_init;
	f->clear    = bn254_fp2_clear;
	f->set      = bn254_fp2_set;
	f->set_str  = bn254_fp2_set_str;
	f->get_str  = bn254_fp2_get_str;
	f->set_zero = bn254_fp2_set_zero;
	f->set_one  = bn254_fp2_set_one;

	f->add  = bn254_fp2_add;
	f->sub  = bn254_fp2_sub;
	f->neg  = bn254_fp2_neg;
	f->mul  = bn254_fp2_mul;
	f->sqr  = bn254_fp2_sqr;
	f->inv  = bn254_fp2_inv;
	f->pow  = bn254_fp2_pow;
	f->sqrt = bn254_fp2_sqrt;

	f->is_zero = bn254_fp2_is_zero;
	f->is_one  = bn254_fp2_is_one;
	f->is_sqr  = bn254_fp2_is_sqr;
	f->cmp     = bn254_fp2_cmp;

	f->random = bn254_fp2_random;

	f->to_oct   = bn254_fp2_to_oct;
	f->from_oct = bn254_fp2_from_oct;

	//-----------------------------------------
	//  set base field
	//-----------------------------------------
	f->base = (struct ec_field_st *)malloc(sizeof(struct ec_field_st));

	field_init(f->base, "bn254_fp");

	//-----------------------------------------
	//  characteristic of prime field
	//-----------------------------------------
	//    p = 36t^4 + 36t^3 + 24t^2 + 6t + 1
	//    t = 2^62 - 2^54 + 2^44
	//-----------------------------------------
	mpz_init(f->order);
	mpz_mul(f->order, f->base->order, f->base->order);

	//-----------------------------------------
	//  Irreducible polynomial: x^2 + 5
	//-----------------------------------------
	f->irre_poly_num = 1;
	f->irre_poly_deg = 2;

	f->irre_poly = (Element *)malloc(sizeof(Element));
	element_init(f->irre_poly[0], f->base);
	element_set_str(f->irre_poly[0], "5");

	//-----------------------------------------
	//  pre-computation for square root
	//-----------------------------------------
	bn254_fp2_precomp(f);

	//----------------------------------
	//  temporary element init
	//----------------------------------
	f->tmp = (Element *)malloc(sizeof(Element)*TMP_NUM);
	for(i=0;i<TMP_NUM;i++){ element_init(f->tmp[i], f); }

	return;
}

//----------------------------------------------
//  function creating field fp254_6
//----------------------------------------------
void ec_bn254_fp6_new(Field f)
{
	int i;

	f->type = Field_fpn;

	set_field_name(f, "bn254_fp6");

	f->ID = bn254_fp6;

	f->str_len = 390;
	f->oct_len = 192;

	//------------------------------
	//  set pointer of function
	//------------------------------
	f->init     = bn254_fp6_init;
	f->clear    = bn254_fp6_clear;
	f->set      = bn254_fp6_set;
	f->set_str  = bn254_fp6_set_str;
	f->get_str  = bn254_fp6_get_str;
	f->set_zero = bn254_fp6_set_zero;
	f->set_one  = bn254_fp6_set_one;

	f->add  = bn254_fp6_add;
	f->sub  = bn254_fp6_sub;
	f->neg  = bn254_fp6_neg;
	f->mul  = bn254_fp6_mul;
	f->sqr  = bn254_fp6_sqr;
	f->inv  = bn254_fp6_inv;
	f->pow  = bn254_fp2_pow;
	f->sqrt = bn254_fp2_sqrt;

	f->is_zero = bn254_fp6_is_zero;
	f->is_one  = bn254_fp6_is_one;
	f->is_sqr  = bn254_fp6_is_sqr;

	f->cmp = bn254_fp6_cmp;

	f->random = bn254_fp6_random;

	f->to_oct   = bn254_fp6_to_oct;
	f->from_oct = bn254_fp6_from_oct;

	//-----------------------------------------
	//  set base field
	//-----------------------------------------
	f->base = (struct ec_field_st *)malloc(sizeof(struct ec_field_st));

	field_init(f->base, "bn254_fp2");

	//-----------------------------------------
	//  characteristic of prime field
	//-----------------------------------------
	//    p = 36t^4 + 36t^3 + 24t^2 + 6t + 1
	//    t = 2^62 - 2^54 + 2^44
	//-----------------------------------------
	mpz_init(f->order);
	mpz_mul(f->order, f->base->order, f->base->order);
	mpz_mul(f->order, f->order, f->base->order);

	//-----------------------------------------
	//  irreducible polynomial: y^3 - x
	//-----------------------------------------
	f->irre_poly_num = 1;
	f->irre_poly_deg = 3;

	f->irre_poly = (Element *)malloc(sizeof(Element));

	element_init(f->irre_poly[0], f->base);
	element_set_str(f->irre_poly[0], "0 -1");

	//-----------------------------------------
	//  pre-computation for square root
	//-----------------------------------------
	bn254_fp6_precomp(f);

	//----------------------------------
	//  temporary element init
	//----------------------------------
	f->tmp = (Element *)malloc(sizeof(Element)*TMP_NUM);
	for(i=0;i<TMP_NUM;i++){ element_init(f->tmp[i], f); }

	return;
}

//----------------------------------------------
//  function creating field fp254_6
//----------------------------------------------
void ec_bn254_fp12_new(Field f)
{
	int i;

	f->type = Field_fpn;

	set_field_name(f, "bn254_fp12");

	f->ID = bn254_fp12;

	f->str_len = 780;
	f->oct_len = 384;

	//------------------------------
	//  set pointer of function
	//------------------------------
	f->init     = bn254_fp12_init;
	f->clear    = bn254_fp12_clear;
	f->set      = bn254_fp12_set;
	f->set_str  = bn254_fp12_set_str;
	f->get_str  = bn254_fp12_get_str;
	f->set_zero = bn254_fp12_set_zero;
	f->set_one  = bn254_fp12_set_one;

	f->add  = bn254_fp12_add;
	f->sub  = bn254_fp12_sub;
	f->neg  = bn254_fp12_neg;
	f->mul  = bn254_fp12_mul;
	f->sqr  = bn254_fp12_sqr;
	f->inv  = bn254_fp12_inv;
	f->pow  = bn254_fp12_pow_naf;
	f->sqrt = bn254_fp2_sqrt;

	f->is_zero = bn254_fp12_is_zero;
	f->is_one  = bn254_fp12_is_one;
	f->is_sqr  = bn254_fp12_is_sqr;
	f->cmp     = bn254_fp12_cmp;

	f->random = bn254_fp12_random;

	f->to_oct   = bn254_fp12_to_oct;
	f->from_oct = bn254_fp12_from_oct;

	//-----------------------------------------
	//  set base field
	//-----------------------------------------
	f->base = (struct ec_field_st *)malloc(sizeof(struct ec_field_st));

	field_init(f->base, "bn254_fp6");

	//-----------------------------------------
	//  characteristic of prime field
	//-----------------------------------------
	//    p = 36t^4 + 36t^3 + 24t^2 + 6t + 1
	//    t = 2^62 - 2^54 + 2^44
	//-----------------------------------------
	mpz_init(f->order);
	mpz_mul(f->order, f->base->order, f->base->order);

	//-----------------------------------------
	//  irreducible polynomial: z^2 - y
	//-----------------------------------------
	f->irre_poly_num = 1;
	f->irre_poly_deg = 2;

	f->irre_poly = (Element *)malloc(sizeof(Element));

	element_init(f->irre_poly[0], f->base);
	element_set_str(f->irre_poly[0], "0 0 -1 0 0 0");

	//-----------------------------------------
	//  pre-computation for square root
	//-----------------------------------------
	bn254_fp12_precomp(f);

	//----------------------------------
	//  temporary element init
	//----------------------------------
	f->tmp = (Element *)malloc(sizeof(Element)*TMP_NUM);
	for(i=0;i<TMP_NUM;i++){ element_init(f->tmp[i], f); }

	return;
}

//----------------------------------------------
//  function release field "bn254"
//----------------------------------------------
void ec_bn254_field_clear(Field f)
{
	unsigned int i, j;

	if( f->precomp != NULL )
	{
		field_precomp_p precomp = (field_precomp_p)(f->precomp);

		field_precomp_sqrt_p ps = precomp->ps;
		field_precomp_frob_p pf = precomp->pf;

		if( ps != NULL )
		{
			element_clear(ps->n_v);
			mpz_clear(ps->v);
			free(ps);
		}

		if( pf != NULL )
		{
			for(i=0;i<pf->glen1;i++){ element_clear(pf->gamma1[i]); }
			for(i=0;i<pf->glen2;i++){ element_clear(pf->gamma2[i]); }
			for(i=0;i<pf->glen3;i++){ element_clear(pf->gamma3[i]); }
			free(pf->gamma1);
			free(pf->gamma2);
			free(pf->gamma3);
			free(pf);
		}
		SAFE_FREE(f->precomp);
	}

	if( f->irre_poly != NULL )
	{
		for(i=0, j=f->irre_poly_num; i<j; i++) { element_clear(f->irre_poly[i]); }
		f->irre_poly_num = 0;
		f->irre_poly_deg = 0;
		SAFE_FREE(f->irre_poly);
	}

	if( f->tmp != NULL )
	{
		for(i=0;i<TMP_NUM;i++){ element_clear(f->tmp[i]); }
		SAFE_FREE(f->tmp);
	}

	if( f->base != NULL )
	{
		field_clear(f->base);
		SAFE_FREE(f->base);
	}

	mpz_clear(f->order);

	SAFE_FREE(f->field_name);

	f->str_len = 0;
	f->oct_len = 0;

	f->init = NULL;
	f->clear = NULL;
	f->set = NULL;
	f->set_str = NULL;
	f->get_str = NULL;
	f->set_zero = NULL;
	f->set_one = NULL;
	f->add = NULL;
	f->sub = NULL;
	f->neg = NULL;
	f->mul = NULL;
	f->sqr = NULL;
	f->inv = NULL;
	f->pow = NULL;
	f->sqrt = NULL;
	f->is_zero = NULL;
	f->is_one = NULL;
	f->is_sqr = NULL;
	f->cmp = NULL;
	f->random = NULL;
	f->to_oct = NULL;
	f->from_oct = NULL;
}

//----------------------------------------------
//  function generating elliptic curve method
//----------------------------------------------
void ec_bn254_fp_method_new(EC_METHOD method)
{
	method->point_init  = ec_bn254_fp_point_init;
	method->point_clear = ec_bn254_fp_point_clear;
	method->point_set = ec_bn254_fp_point_set;
	method->point_set_infinity = ec_bn254_fp_point_set_infinity;
	method->point_set_str = ec_bn254_fp_point_set_str;
	method->point_get_str = ec_bn254_fp_point_get_str;

	method->add = ec_bn254_fp_add;
	method->dob = ec_bn254_fp_dob;
	method->neg = ec_bn254_fp_neg;
	method->sub = ec_bn254_fp_sub;

#ifdef ENABLE_FASTALG
	method->mul = ec_bn254_fp_mul_end;
#else
	method->mul = ec_bn254_fp_mul_naf;
#endif

	method->is_infinity = ec_bn254_fp_is_infinity;
	method->is_on_curve = ec_bn254_fp_is_on_curve;
	method->cmp = ec_bn254_fp_cmp;

	method->make_affine = ec_bn254_fp_make_affine;
	method->map_to_point = ec_bn254_fp_map_to_point;
	method->random = ec_bn254_fp_random;
	method->to_oct = ec_bn254_fp_to_oct;
	method->from_oct = ec_bn254_fp_from_oct;
}

//----------------------------------------------
//  function generating elliptic curve group
//----------------------------------------------
void ec_bn254_fp_group_new(EC_GROUP ec)
{
	ec->type = Curve_BN;

	set_curve_name(ec, "ec_bn254_p");

	ec->ID = ec_bn254_fp;

	ec->str_len = 132;
	ec->oct_len = 65;

	ec->field = (struct ec_field_st *)malloc(sizeof(struct ec_field_st));
	field_init(ec->field, "bn254_fp");

	ec->method = (struct ec_method_st *)malloc(sizeof(struct ec_method_st));
	ec_bn254_fp_method_new(ec->method);

	point_init(ec->generator, ec);

	element_init(ec->a, ec->field);
	element_init(ec->b, ec->field);

	element_set_zero(ec->a);
	element_set_str(ec->b, "5");

	mpz_init_set_str(ec->order, "2370FB049D410FBE4E761A9886E502411DC1AF70120000017E80600000000001", 16);
	mpz_init_set_str(ec->trace, "5F408FD0060000000000000000000001", 16);
	mpz_init_set_str(ec->cofactor, "1", 16);

	ec_bn254_fp_init_ec_data(ec);
}

//----------------------------------------------
//  function generating elliptic curve method
//----------------------------------------------
void ec_bn254_tw_method_new(EC_METHOD method)
{
	method->point_init  = ec_bn254_fp_point_init;
	method->point_clear = ec_bn254_fp_point_clear;
	method->point_set = ec_bn254_fp_point_set;
	method->point_set_infinity = ec_bn254_fp_point_set_infinity;
	method->point_set_str = ec_bn254_fp2_point_set_str;
	method->point_get_str = ec_bn254_fp2_point_get_str;

	method->add = ec_bn254_fp2_add;
	method->dob = ec_bn254_fp2_dob;
	method->neg = ec_bn254_fp2_neg;
	method->sub = ec_bn254_fp2_sub;

#ifdef ENABLE_FASTALG
	method->mul = ec_bn254_fp2_mul_end;
#else
	method->mul = ec_bn254_fp2_mul;
#endif

	method->is_infinity = ec_bn254_fp_is_infinity;
	method->is_on_curve = ec_bn254_fp_is_on_curve;
	method->cmp = ec_bn254_fp_cmp;

	method->make_affine = ec_bn254_fp2_make_affine;
	method->map_to_point = ec_bn254_fp2_map_to_point;
	method->random = ec_bn254_fp2_random;
	method->to_oct = ec_bn254_fp2_to_oct;
	method->from_oct = ec_bn254_fp2_from_oct;
}

//----------------------------------------------
//  function generating elliptic curve group
//----------------------------------------------
void ec_bn254_tw_group_new(EC_GROUP ec)
{
	ec->type = Curve_BN;

	set_curve_name(ec, "ec_bn254_tw");

	ec->ID = ec_bn254_fp2;

	ec->str_len = 262;
	ec->oct_len = 129;

	ec->field = (struct ec_field_st *)malloc(sizeof(struct ec_field_st));
	field_init(ec->field, "bn254_fp2");

	ec->method = (struct ec_method_st *)malloc(sizeof(struct ec_method_st));
	ec_bn254_tw_method_new(ec->method);

	point_init(ec->generator, ec);

	element_init(ec->a, ec->field);
	element_init(ec->b, ec->field);

	element_set_zero(ec->a);
	element_set_str(ec->b, "0 -1");

	mpz_init_set_str(ec->order, "2370FB049D410FBE4E761A9886E502411DC1AF70120000017E80600000000001", 16);
	mpz_init_set_str(ec->trace, "2370FB049D410FBDC02400000000000000000000000000000000000000000001", 16);
	mpz_init_set_str(ec->cofactor, "2370FB049D410FBE4E761A9886E50241DC42CF101E0000017E80600000000001", 16);

	ec_bn254_fp2_init_ec_data(ec);
}

//----------------------------------------------
//  clear curve group : ec_bn254
//----------------------------------------------
void ec_bn254_group_clear(EC_GROUP ec)
{
	point_clear(ec->generator);

	element_clear(ec->a);
	element_clear(ec->b);

	mpz_clear(ec->order);
	mpz_clear(ec->trace);
	mpz_clear(ec->cofactor);

	if( ec->ID == ec_bn254_fp )
	{
		ec_bn254_fp_clear_ec_data(ec);
	}
	else if( ec->ID == ec_bn254_fp2 )
	{
		ec_bn254_fp2_clear_ec_data(ec);
	}

	SAFE_FREE(ec->method);

	field_clear(ec->field);
	SAFE_FREE(ec->field);

	SAFE_FREE(ec->curve_name);

	ec->str_len = 0;
	ec->oct_len = 0;
}

//-------------------------------------------
// pairing group : Init, Clear
//-------------------------------------------
void ec_bn254_pairing_new(EC_PAIRING p)
{
	p->type = Pairing_ECBN254;

	p->pairing = ec_bn254_pairing;

	curve_init(p->g1, "ec_bn254_fp");
	curve_init(p->g2, "ec_bn254_tw");

	field_init(p->g3, "bn254_fp12");

	ec_bn254_pairing_precomp(p);
}

void ec_bn254_pairing_clear(EC_PAIRING p)
{
	if( p->precomp != NULL )
	{
		pairing_precomp_p precomp = (pairing_precomp_p)(p->precomp);

		free(precomp->si);
		free(precomp->ti);
		free(precomp);
		p->precomp = NULL;
	}

	p->pairing = NULL;

	curve_clear(p->g1);
	curve_clear(p->g2);
	field_clear(p->g3);
}
