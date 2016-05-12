//===================================================
//  Elliptic curve operation
//
//   2015.xx.xx
//===================================================
#pragma once

#include <tepla/ec/field.h>

#ifdef __cplusplus
extern "C" {
#endif

//---------------------------------------------------
// Elliptic Curve Type
//---------------------------------------------------
typedef enum
{
    Curve_BN,

} ECType;

//---------------------------------------------------
// EC point structure
//---------------------------------------------------
typedef struct ec_point_st
{
    const struct ec_group_st *ec;

    Element x;
    Element y;
    Element z;

    int isinfinity;

} EC_POINT[1];

//---------------------------------------------------
// EC group structure E/K: Y^2 = X^3 + aX + b
//---------------------------------------------------
typedef struct ec_group_st
{
    ECType type;

    char* curve_name;

    int ID;

    int str_len;
    int oct_len;

    void (*curve_init)(struct ec_group_st *ec);
    void (*curve_clear)(struct ec_group_st *ec);

    struct ec_field_st *field;
    struct ec_method_st *method;

    EC_POINT generator; // point

    Element a;  // coefficient a
    Element b;  // coefficient b

    mpz_t   order;    // prime order of curve : order | #E(K)
    mpz_t   trace;    // trace of frobenius
    mpz_t   cofactor; // cofactor for order : #E(K) = {cofactor} * {order}

    void *ec_data; // some values for computation

} EC_GROUP[1];

//---------------------------------------------------
// EC method structure
//---------------------------------------------------
typedef struct ec_method_st
{
    void (*point_init)(EC_POINT p);
    void (*point_clear)(EC_POINT p);
    void (*point_set)(EC_POINT z, const EC_POINT x);
    void (*point_set_str)(EC_POINT p, const char* s);
    void (*point_set_infinity)(EC_POINT z);
    void (*point_get_str)(char *s, const EC_POINT p);

    void (*add)(EC_POINT R, const EC_POINT P, const EC_POINT Q);
    void (*dob)(EC_POINT Q, const EC_POINT P);
    void (*neg)(EC_POINT Q, const EC_POINT P);
    void (*sub)(EC_POINT R, const EC_POINT P, const EC_POINT Q);
    void (*mul)(EC_POINT Q, const mpz_t s, const EC_POINT P);

    int (*is_infinity)(const EC_POINT P);
    int (*is_on_curve)(const EC_POINT P);
    int (*cmp)(const EC_POINT P, const EC_POINT Q);

    void (*make_affine)(EC_POINT Q, const EC_POINT P);
    void (*map_to_point)(EC_POINT P, const char *s, size_t slen, int t);
    void (*random)(EC_POINT P);
    void (*to_oct)(unsigned char* os, size_t *size, const EC_POINT P);
    void (*from_oct)(EC_POINT P, const unsigned char* os, size_t size);

} EC_METHOD[1];

//---------------------------------------------------
// Pairing type
//---------------------------------------------------
typedef enum
{
    Pairing_ECBN254a,
    Pairing_ECBN254b,

} PairingType;

//---------------------------------------------------
// pairing structure
//---------------------------------------------------
typedef struct ec_pairing_st
{
    PairingType type;

    char* pairing_name;

    void (*pairing)(Element z, const EC_POINT x, const EC_POINT y, const struct ec_pairing_st* p);
    void (*pairing_double)(Element z, const EC_POINT x1, const EC_POINT y1, const EC_POINT x2, const EC_POINT y2, const struct ec_pairing_st* p);

    EC_GROUP g1;
    EC_GROUP g2;

    Field    g3;

    void* precomp;

} EC_PAIRING[1];

//---------------------------------------------------
//  functions for Elliptic Curve
//---------------------------------------------------
void curve_init(EC_GROUP ec, const char *param);
void curve_clear(EC_GROUP ec);

const char*  curve_get_name(const EC_GROUP ec);
const mpz_t* curve_get_order(const EC_GROUP ec);

//---------------------------------------------------
//  functions for Point on Elliptic Curve
//---------------------------------------------------
void point_init(EC_POINT p, const EC_GROUP ec);
void point_clear(EC_POINT p);

void point_set(EC_POINT P, const EC_POINT Q);
void point_set_str(EC_POINT P, const char *s);
void point_set_xy(EC_POINT P, const Element x, const Element y);
void point_set_infinity(EC_POINT P);
void point_get_str(char *s, const EC_POINT P);

void point_add(EC_POINT R, const EC_POINT P, const EC_POINT Q);
void point_dob(EC_POINT Q, const EC_POINT P);
void point_neg(EC_POINT Q, const EC_POINT P);
void point_sub(EC_POINT R, const EC_POINT P, const EC_POINT Q);
void point_mul(EC_POINT Q, const mpz_t s, const EC_POINT P);

int  point_is_infinity(const EC_POINT P);
int  point_is_on_curve(const EC_POINT P);
int  point_cmp(const EC_POINT P, const EC_POINT Q);

void point_make_affine(EC_POINT Q, const EC_POINT P);
void point_map_to_point(EC_POINT P, const char *s, size_t slen, int t);
void point_random(EC_POINT P);
void point_to_oct(unsigned char* os, size_t *size, EC_POINT P);
void point_from_oct(EC_POINT P, const unsigned char *os, size_t size);

int  point_get_str_length(const EC_POINT P);
int  point_get_oct_length(const EC_POINT P);

void point_print(const EC_POINT P);

//---------------------------------------------------
//  functions for pairing
//---------------------------------------------------
void pairing_init(EC_PAIRING p, char *param);
void pairing_clear(EC_PAIRING p);

void pairing_map(Element g, const EC_POINT P, const EC_POINT Q, const EC_PAIRING p);
void pairing_double_map(Element g, const EC_POINT P1, const EC_POINT Q1, const EC_POINT P2, const EC_POINT Q2, const EC_PAIRING p);

const mpz_t* pairing_get_order(const EC_PAIRING p);
const char* pairing_get_name(const EC_PAIRING p);

#ifdef __cplusplus
}
#endif
