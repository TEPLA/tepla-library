//----------------------------------------------------
//  Header file for ec_bn254
//----------------------------------------------------
#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <tepla/ec.h>
#include <tepla/hash.h>

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif

#ifndef NULL
#define NULL 0
#endif

//---------------------------------------------------
// Finite Field (BN254) ID
//---------------------------------------------------
typedef enum
{
	bn254_fp,
	bn254_fp2,
	bn254_fp6,
	bn254_fp12,

} BN254_FieldType;

//---------------------------------------------------
//  precomputation values for sqrt
//---------------------------------------------------
//    (p^m-1) = 2^e *v
//    n_v = n^v : (n/p) = -1
//---------------------------------------------------
typedef struct ec_field_precomp_sqrt_st
{
	int   e;
	mpz_t v;
	Element  n_v;

} *field_precomp_sqrt_p;

//---------------------------------------------------
//  precomputation values for frobenius map
//---------------------------------------------------
typedef struct ec_field_precomp_frob_st
{
	size_t glen1;
	size_t glen2;
	size_t glen3;

	Element *gamma1;
	Element *gamma2;
	Element *gamma3;

} *field_precomp_frob_p;

//---------------------------------------------------
// structure for precomputation values
//---------------------------------------------------
typedef struct ec_field_precomp_st
{
	field_precomp_sqrt_p  ps;
	field_precomp_frob_p  pf;

} *field_precomp_p;

//---------------------------------------------------
// Elliptic Curve (BN254) ID
//---------------------------------------------------
typedef enum
{
	ec_bn254_fp,
	ec_bn254_fp2,

} BN254_CurveType;

//---------------------------------------------------
// structure for Elliptic Curve values
//---------------------------------------------------
typedef struct ec_bn254_fp_ec_data_st
{
	Element beta;

	mpz_t n;
	mpz_t n2;
	mpz_t a1, a2;
	mpz_t b1, b2;

} *ec_data_fp;

//---------------------------------------------------
// structure for Elliptic Curve values
//---------------------------------------------------
typedef struct ec_bn254_fp2_ec_data_st
{
	mpz_t _6x;
	mpz_t _6x2;

	Element vfrobx, vfroby;
	Element vfrobx2, vfroby2;
	Element vfrobx3, vfroby3;

} *ec_data_fp2;

//---------------------------------------------------
//  structure for precomputation value
//---------------------------------------------------
typedef struct ec_pairing_precomp_st
{
	size_t slen;
	int *si;

	size_t tlen; // for calculating f^t
	int *ti;     // for calculating f^t

} *pairing_precomp_p;

//----------------------------------------------
// declaration function of field bn254_fp
//----------------------------------------------
void bn254_fp_init(Element x);
void bn254_fp_clear(Element x);
void bn254_fp_set(Element x, const Element y);
void bn254_fp_set_str(Element x, const char *str);
void bn254_fp_get_str(char *str, const Element x);
void bn254_fp_set_zero(Element x);
void bn254_fp_set_one(Element x);
void bn254_fp_mod(Element z, const Element x);
void bn254_fp_add(Element z, const Element x, const Element y);
void bn254_fp_addn(Element z, const Element x, const Element y);
void bn254_fp_addp(Element z, const Element x);
void bn254_fp_dob(Element z, const Element x);
void bn254_fp_tri(Element z, const Element x);
void bn254_fp_neg(Element z, const Element x);
void bn254_fp_sub(Element z, const Element x, const Element y);
void bn254_fp_subn(Element z, const Element x, const Element y);
void bn254_fp_mul(Element z, const Element x, const Element y);
void bn254_fp_muln(Element z, const Element x, const Element y);
void bn254_fp_sqr(Element z, const Element x);
void bn254_fp_inv(Element z, const Element x);
void bn254_fp_pow(Element z, const Element x, const mpz_t exp);
int  bn254_fp_sqrt(Element z, const Element x);
int  bn254_fp_is_zero(const Element x);
int  bn254_fp_is_one(const Element x);
int  bn254_fp_is_sqr(const Element x);
int  bn254_fp_is_sqr_general(const Element x);
int  bn254_fp_cmp(const Element x, const Element y);
void bn254_fp_random(Element z);
void bn254_fp_to_oct(unsigned char *os, size_t *size, const Element x);
void bn254_fp_from_oct(Element z, const unsigned char *os, const size_t size);

//----------------------------------------------
// declaration function of field bn254_fp2
//----------------------------------------------
void bn254_fp2_init(Element x);
void bn254_fp2_clear(Element x);
void bn254_fp2_set(Element x, const Element y);
void bn254_fp2_set_fp(Element z, const Element x, const Element y);
void bn254_fp2_set_str(Element x, const char *s);
void bn254_fp2_get_str(char *s, const Element x);
void bn254_fp2_set_zero(Element x);
void bn254_fp2_set_one(Element x);
void bn254_fp2_add(Element z, const Element x, const Element y);
void bn254_fp2_dob(Element z, const Element x);
void bn254_fp2_tri(Element z, const Element x);
void bn254_fp2_neg(Element z, const Element x);
void bn254_fp2_sub(Element z, const Element x, const Element y);
void bn254_fp2_mul(Element z, const Element x, const Element y);
void bn254_fp2_mul_p(Element z, const Element x, const Element y);
void bn254_fp2_sqr(Element z, const Element x);
void bn254_fp2_xi_mul(Element z, const Element x);
void bn254_fp2_inv(Element z, const Element x);
void bn254_fp2_pow(Element z, const Element x, const mpz_t exp);
int  bn254_fp2_sqrt(Element z, const Element x);
void bn254_fp2_frob_p(Element z, const Element x);
void bn254_fp2_conj(Element z, const Element x);
int  bn254_fp2_is_zero(const Element x);
int  bn254_fp2_is_one(const Element x);
int  bn254_fp2_is_sqr(const Element x);
int  bn254_fp2_cmp(const Element x, const Element y);
void bn254_fp2_precomp_sqrt(field_precomp_sqrt_p ps, const Field f);
void bn254_fp2_precomp(Field f);
void bn254_fp2_random(Element z);
void bn254_fp2_to_oct(unsigned char *os, size_t *size, const Element x);
void bn254_fp2_from_oct(Element z, const unsigned char *os, const size_t size);

//----------------------------------------------
// declaration function of field bn254_fp6
//----------------------------------------------
void bn254_fp6_init(Element x);
void bn254_fp6_clear(Element x);
void bn254_fp6_set(Element x, const Element y);
void bn254_fp6_set_fp2(Element z, const Element w, const Element x, const Element y);
void bn254_fp6_set_str(Element x, const char *s);
void bn254_fp6_get_str(char *s, const Element x);
void bn254_fp6_set_zero(Element x);
void bn254_fp6_set_one(Element x);
void bn254_fp6_add(Element z, const Element x, const Element y);
void bn254_fp6_dob(Element z, const Element x);
void bn254_fp6_tri(Element z, const Element x);
void bn254_fp6_neg(Element z, const Element x);
void bn254_fp6_sub(Element z, const Element x, const Element y);
void bn254_fp6_mul(Element z, const Element x, const Element y);
void bn254_fp6_gm_mul(Element z, const Element x);
void bn254_fp6_mul_fp2(Element z, const Element x, const Element y);
void bn254_fp6_mul_fp2_2(Element z, const Element y, const Element x1, const Element x2);
void bn254_fp6_conj(Element z, const Element x);
void bn254_fp6_sqr(Element z, const Element x);
void bn254_fp6_inv(Element z, const Element x);
void bn254_fp6_frob_p(Element z, const Element x);
int  bn254_fp6_is_zero(const Element x);
int  bn254_fp6_is_one(const Element x);
int  bn254_fp6_is_sqr(const Element x);
int  bn254_fp6_cmp(const Element x, const Element y);
void bn254_fp6_precomp(Field f);
void bn254_fp6_random(Element z);
void bn254_fp6_to_oct(unsigned char *os, size_t *size, const Element x);
void bn254_fp6_from_oct(Element z, const unsigned char *os, const size_t size);

//----------------------------------------------
// declaration function of field bn254_fp12
//----------------------------------------------
void bn254_fp12_init(Element x);
void bn254_fp12_clear(Element x);
void bn254_fp12_set(Element x, const Element y);
void bn254_fp12_set_fp6(Element z, const Element x, const Element y);
void bn254_fp12_set_str(Element x, const char *s);
void bn254_fp12_get_str(char *s, const Element x);
void bn254_fp12_set_zero(Element x);
void bn254_fp12_set_one(Element x);
void bn254_fp12_add(Element z, const Element x, const Element y);
void bn254_fp12_dob(Element z, const Element x);
void bn254_fp12_tri(Element z, const Element x);
void bn254_fp12_neg(Element z, const Element x);
void bn254_fp12_sub(Element z, const Element x, const Element y);
void bn254_fp12_mul(Element z, const Element x, const Element y);
void bn254_fp12_mul_L(Element z, Element x0, Element x1, Element x2);
void bn254_fp12_sqr(Element z, const Element x);
void bn254_fp12_inv(Element z, const Element x);
void bn254_fp12_pow(Element z, const Element x, const mpz_t exp);
void bn254_fp12_pow_naf(Element z, const Element x, const mpz_t exp);
void bn254_fp12_frob_p(Element z, const Element x);
void bn254_fp12_frob_p2(Element z, const Element x);
void bn254_fp12_frob_p3(Element z, const Element x);
void bn254_fp12_conj(Element z, const Element x);
void bn254_fp4_sqr(Element c0, Element c1, const Element a0, const Element a1);
void bn254_fp12_sqr_forpairing(Element z, const Element x);
void bn254_fp12_pow_forpairing(Element z, const Element x, const int *t, int tlen);
int  bn254_fp12_is_zero(const Element x);
int  bn254_fp12_is_one(const Element x);
int  bn254_fp12_is_sqr(const Element x);
int  bn254_fp12_cmp(const Element x, const Element y);
void bn254_fp12_precomp(Field f);
void bn254_fp12_random(Element z);
void bn254_fp12_to_oct(unsigned char *os, size_t *size, const Element x);
void bn254_fp12_from_oct(Element z, const unsigned char *os, const size_t size);

//----------------------------------------------
// declaration function of elliptic curve
//----------------------------------------------
void ec_bn254_fp_point_init(EC_POINT p);
void ec_bn254_fp_point_clear(EC_POINT p);
void ec_bn254_fp_point_set(EC_POINT z, const EC_POINT x);
void ec_bn254_fp_point_set_str(EC_POINT z, const char* s);
void ec_bn254_fp_point_set_infinity(EC_POINT z);
void ec_bn254_fp_point_get_str(char *s, const EC_POINT z);
void ec_bn254_fp_add(EC_POINT z, const EC_POINT x, const EC_POINT y);
void ec_bn254_fp_dob(EC_POINT z, const EC_POINT x);
void ec_bn254_fp_add_formul(EC_POINT z, const EC_POINT x, const EC_POINT y);
void ec_bn254_fp_dob_formul(EC_POINT z, const EC_POINT x);
void ec_bn254_fp_neg(EC_POINT z, const EC_POINT x);
void ec_bn254_fp_sub(EC_POINT z, const EC_POINT x, const EC_POINT y);
void ec_bn254_fp_mul_affine(EC_POINT z, const mpz_t s, const EC_POINT x);
void ec_bn254_fp_mul(EC_POINT z, const mpz_t s, const EC_POINT x);
void ec_bn254_fp_mul_naf(EC_POINT z, const mpz_t s, const EC_POINT x);
void ec_bn254_fp_mul_end(EC_POINT z, const mpz_t s, const EC_POINT x);
int  ec_bn254_fp_is_infinity(const EC_POINT P);
int  ec_bn254_fp_is_on_curve(const EC_POINT P);
int  ec_bn254_fp_cmp(const EC_POINT x, const EC_POINT y);
void ec_bn254_fp_make_affine(EC_POINT z, const EC_POINT x);
void ec_bn254_fp_map_to_point(EC_POINT z, const char *s, size_t slen, int t);
void ec_bn254_fp_point_endomorphism(EC_POINT Q, const EC_POINT P);
void ec_bn254_fp_random(EC_POINT z);
void ec_bn254_fp_to_oct(unsigned char *os, size_t *size, const EC_POINT z);
void ec_bn254_fp_from_oct(EC_POINT z, const unsigned char *os, size_t size);

void generate_naf(int *naf, int *len, const mpz_t s);
void cat_int_str(unsigned char *os, size_t *oslen, const mpz_t i, const unsigned char *s, const size_t slen);

void ec_bn254_fp_decompose_scalar_init(mpz_t a1, mpz_t a2, mpz_t b1, mpz_t b2, const mpz_t n, const mpz_t l);

void ec_bn254_fp_init_ec_data(EC_GROUP ec);
void ec_bn254_fp_clear_ec_data(EC_GROUP ec);

//----------------------------------------------
// declaration function of elliptic curve
//----------------------------------------------
void ec_bn254_fp2_point_init(EC_POINT p);
void ec_bn254_fp2_point_clear(EC_POINT p);
void ec_bn254_fp2_point_set(EC_POINT z, const EC_POINT x);
void ec_bn254_fp2_point_set_str(EC_POINT z, const char* s);
void ec_bn254_fp2_point_set_infinity(EC_POINT z);
void ec_bn254_fp2_point_get_str(char *s, const EC_POINT z);
void ec_bn254_fp2_add(EC_POINT z, const EC_POINT x, const EC_POINT y);
void ec_bn254_fp2_dob(EC_POINT z, const EC_POINT x);
void ec_bn254_fp2_add_formul(EC_POINT z, const EC_POINT x, const EC_POINT y);
void ec_bn254_fp2_dob_formul(EC_POINT z, const EC_POINT x);
void ec_bn254_fp2_neg(EC_POINT z, const EC_POINT x);
void ec_bn254_fp2_sub(EC_POINT z, const EC_POINT x, const EC_POINT y);
void ec_bn254_fp2_mul(EC_POINT z, const mpz_t s, const EC_POINT x);
void ec_bn254_fp2_mul_naf(EC_POINT z, const mpz_t s, const EC_POINT x);
void ec_bn254_fp2_mul_end(EC_POINT z, const mpz_t s, const EC_POINT x);
void ec_bn254_fp2_frob_p(EC_POINT Q, const EC_POINT P);
int  ec_bn254_fp2_is_infinity(const EC_POINT P);
int  ec_bn254_fp2_is_on_curve(const EC_POINT P);
int  ec_bn254_fp2_cmp(const EC_POINT x, const EC_POINT y);
void ec_bn254_fp2_make_affine(EC_POINT z, const EC_POINT x);
void ec_bn254_fp2_map_to_point(EC_POINT z, const char *s, size_t slen,int t);
void ec_bn254_fp2_random(EC_POINT z);
void ec_bn254_fp2_to_oct(unsigned char *os, size_t *size, const EC_POINT z);
void ec_bn254_fp2_from_oct(EC_POINT z, const unsigned char *os, size_t size);

void ec_bn254_tw_frob(EC_POINT Q, const EC_POINT P);
void ec_bn254_tw_frob2(EC_POINT Q, const EC_POINT P);
void ec_bn254_tw_frob3(EC_POINT Q, const EC_POINT P);

void ec_bn254_fp2_init_ec_data(EC_GROUP ec);
void ec_bn254_fp2_clear_ec_data(EC_GROUP ec);

//----------------------------------------------
// declaration function of pairing
//----------------------------------------------
void ec_bn254_pairing_precomp(EC_PAIRING p);
void ec_bn254_pairing_dob(EC_POINT T, Element l0, Element l3, Element l4, const EC_POINT P);
void ec_bn254_pairing_add(EC_POINT T, Element l0, Element l3, Element l4, const EC_POINT P, const EC_POINT Q);
void ec_bn254_pairing_miller(Element z, const EC_POINT Q, const EC_POINT P, const EC_PAIRING p);
void ec_bn254_pairing_finalexp(Element z, const Element x, const EC_PAIRING p);
void ec_bn254_pairing(Element z, const EC_POINT Q, const EC_POINT P, const EC_PAIRING p);
