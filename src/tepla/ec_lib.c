//==================================================
// Function which making instance of structures
//==================================================
#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <tepla/ec.h>

#include "ec_bn254/ec.h"

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif

#ifndef NULL
#define NULL 0
#endif

#define Field(x)  (x->field)
#define Curve(p)  (p->ec->method)

//============================================
//  Field : Init, Clear
//============================================
void field_init(Field f, const char *param)
{
    // Beuchat's parameter
    if (strcmp(param, "bn254_fpa") == 0)
    {
        f->field_init  = ec_bn254_fpa_new;
        f->field_clear = ec_bn254_field_clear;
    }
    else if (strcmp(param, "bn254_fp2a") == 0)
    {
        f->field_init  = ec_bn254_fp2a_new;
        f->field_clear = ec_bn254_field_clear;
    }
    else if (strcmp(param, "bn254_fp6a") == 0)
    {
        f->field_init  = ec_bn254_fp6a_new;
        f->field_clear = ec_bn254_field_clear;
    }
    else if (strcmp(param, "bn254_fp12a") == 0)
    {
        f->field_init  = ec_bn254_fp12a_new;
        f->field_clear = ec_bn254_field_clear;
    }

    // Aranha's parameter
    else if (strcmp(param, "bn254_fpb") == 0)
    {
        f->field_init  = ec_bn254_fpb_new;
        f->field_clear = ec_bn254_field_clear;
    }
    else if (strcmp(param, "bn254_fp2b") == 0)
    {
        f->field_init  = ec_bn254_fp2b_new;
        f->field_clear = ec_bn254_field_clear;
    }
    else if (strcmp(param, "bn254_fp6b") == 0)
    {
        f->field_init  = ec_bn254_fp6b_new;
        f->field_clear = ec_bn254_field_clear;
    }
    else if (strcmp(param, "bn254_fp12b") == 0)
    {
        f->field_init  = ec_bn254_fp12b_new;
        f->field_clear = ec_bn254_field_clear;
    }

    else
    {
        fprintf(stderr, "We do not support that identity : %s\n", param);
        exit(200);
    }

    f->field_init(f);
}

void field_clear(Field f)
{
    f->field_clear(f);
}

//============================================
//  Get Field Parameters
//============================================
const mpz_t* field_get_char(const Field f)
{
    if (f->base != NULL) {
        return field_get_char(f->base);
    }
    return &(f->order);
}

int field_get_degree(const Field f)
{
    if (f->base != NULL) {
        return f->irre_poly_deg * field_get_degree(f->base);
    }
    return f->irre_poly_deg;
}

const char* field_get_name(const Field f)
{
    return f->field_name;
}

//============================================
//  Element Operation
//============================================
void element_init(Element x, const Field f)
{
    x->field = f;
    f->init(x);
}

void element_clear(Element x)
{
    x->field->clear(x);
    x->field = NULL;
}

void element_set(Element x, const Element y)
{
    Field(x)->set(x, y);
}

void element_set_str(Element x, const char *str)
{
    Field(x)->set_str(x, str);
}

void element_get_str(char *str, const Element x)
{
    Field(x)->get_str(str, x);
}

void element_set_zero(Element x)
{
    Field(x)->set_zero(x);
}

void element_set_one(Element x)
{
    Field(x)->set_one(x);
}

void element_add(Element z, const Element x, const Element y)
{
    Field(x)->add(z, x, y);
}

void element_neg(Element z, const Element x)
{
    Field(x)->neg(z, x);
}

void element_sub(Element z, const Element x, const Element y)
{
    Field(x)->sub(z, x, y);
}

void element_mul(Element z, const Element x, const Element y)
{
    Field(x)->mul(z, x, y);
}

void element_sqr(Element z, const Element x)
{
    Field(x)->sqr(z, x);
}

void element_inv(Element z, const Element x)
{
    Field(x)->inv(z, x);
}

void element_pow(Element z, const Element x, const mpz_t exp)
{
    Field(x)->pow(z, x, exp);
}

int element_sqrt(Element z, const Element x)
{
    return Field(x)->sqrt(z, x);
}

int element_is_zero(const Element x)
{
    return Field(x)->is_zero(x);
}

int element_is_one(const Element x)
{
    return Field(x)->is_one(x);
}

int element_is_sqr(const Element x)
{
    return Field(x)->is_sqr(x);
}

int element_cmp(const Element x, const Element y)
{
    return Field(x)->cmp(x, y);
}

void element_random(Element x)
{
    Field(x)->random(x);
}

void element_to_oct(unsigned char *os, size_t *size, Element x)
{
    Field(x)->to_oct(os, size, x);
}

void element_from_oct(Element z, const unsigned char *os, size_t size)
{
    Field(z)->from_oct(z, os, size);
}

int element_get_str_length(const Element x)
{
    return Field(x)->str_len;
}

int element_get_oct_length(const Element x)
{
    return Field(x)->oct_len;
}

void element_print(const Element x)
{
    int len = element_get_str_length(x);

    char *s = (char*)malloc(sizeof(char) * len);

    element_get_str(s, x);

    printf("element: %s\n", s);

    free(s);
}

//============================================
//  Initialization of Elliptic Curve
//============================================
void curve_init(EC_GROUP ec, const char *param)
{
    if (strcmp(param, "ec_bn254_fpa") == 0)
    {
        ec->curve_init = ec_bn254_fpa_group_new;
        ec->curve_clear = ec_bn254_group_clear;
    }
    else if (strcmp(param, "ec_bn254_twa") == 0)
    {
        ec->curve_init = ec_bn254_twa_group_new;
        ec->curve_clear = ec_bn254_group_clear;
    }
    else if (strcmp(param, "ec_bn254_fpb") == 0)
    {
        ec->curve_init = ec_bn254_fpb_group_new;
        ec->curve_clear = ec_bn254_group_clear;
    }
    else if (strcmp(param, "ec_bn254_twb") == 0)
    {
        ec->curve_init = ec_bn254_twb_group_new;
        ec->curve_clear = ec_bn254_group_clear;
    }
    else
    {
        fprintf(stderr, "We donot support the identity : %s\n", param);
        exit(200);
    }

    ec->curve_init(ec);
}

void curve_clear(EC_GROUP ec)
{
    ec->curve_clear(ec);
}

const char* curve_get_name(const EC_GROUP ec)
{
    return ec->curve_name;
}

const mpz_t* curve_get_order(const EC_GROUP ec)
{
    return &(ec->order);
}

//============================================
//  Point Operation
//============================================
void point_init(EC_POINT p, const EC_GROUP ec)
{
    p->ec = ec;
    ec->method->point_init(p);
}

void point_clear(EC_POINT p)
{
    Curve(p)->point_clear(p);
    p->ec = NULL;
}

void point_set(EC_POINT p, const EC_POINT q)
{
    Curve(p)->point_set(p, q);
}

void point_set_str(EC_POINT p, const char *s)
{
    Curve(p)->point_set_str(p, s);
}

void point_set_xy(EC_POINT p, const Element x, const Element y)
{
    element_set(p->x, x);
    element_set(p->y, y);
    element_set_one(p->z);

    p->isinfinity = FALSE;
}

void point_set_infinity(EC_POINT p)
{
    Curve(p)->point_set_infinity(p);
}

void point_get_str(char *s, const EC_POINT p)
{
    Curve(p)->point_get_str(s, p);
}

void point_add(EC_POINT z, const EC_POINT x, const EC_POINT y)
{
    Curve(x)->add(z, x, y);
}

void point_dob(EC_POINT Q, const EC_POINT P)
{
    Curve(P)->dob(Q, P);
}

void point_neg(EC_POINT Q, const EC_POINT P)
{
    Curve(P)->neg(Q, P);
}

void point_sub(EC_POINT z, const EC_POINT x, const EC_POINT y)
{
    Curve(x)->sub(z, x, y);
}

void point_mul(EC_POINT z, const mpz_t s, const EC_POINT x)
{
    Curve(x)->mul(z, s, x);
}

int point_is_infinity(const EC_POINT x)
{
    return x->isinfinity;
}

int point_is_on_curve(const EC_POINT x)
{
    return Curve(x)->is_on_curve(x);
}

int point_cmp(const EC_POINT x, const EC_POINT y)
{
    return Curve(x)->cmp(x, y);
}

void point_make_affine(EC_POINT z, const EC_POINT x)
{
    Curve(x)->make_affine(z, x);
}

void point_map_to_point(EC_POINT z, const char *s, size_t slen, int t)
{
    Curve(z)->map_to_point(z, s, slen, t);
}

void point_random(EC_POINT P)
{
    Curve(P)->random(P);
}

void point_to_oct(unsigned char *os, size_t *size, EC_POINT P)
{
    Curve(P)->to_oct(os, size, P);
}

void point_from_oct(EC_POINT P, const unsigned char *os, size_t size)
{
    Curve(P)->from_oct(P, os, size);
}

int point_get_str_length(const EC_POINT P)
{
    return P->ec->str_len;
}

int point_get_oct_length(const EC_POINT P)
{
    return P->ec->oct_len;
}

void point_print(const EC_POINT P)
{
    int len = point_get_str_length(P);

    char *s = (char*)malloc(sizeof(char) * len);

    point_get_str(s, P);

    printf("point: %s\n", s);

    free(s);
}

//============================================
//  Pairing operation
//============================================
void pairing_init(EC_PAIRING p, char *param)
{
    if (strcmp(param, "ECBN254a") == 0)
    {
        ec_bn254_pairing_a_new(p);
    }
    else if (strcmp(param, "ECBN254b") == 0)
    {
        ec_bn254_pairing_b_new(p);
    }
    else
    {
        fprintf(stderr, "We donot suppoert the identity : %s\n", param);
        exit(200);
    }
}

void pairing_clear(EC_PAIRING p)
{
    if (p->type == Pairing_ECBN254a || p->type == Pairing_ECBN254b)
    {
        ec_bn254_pairing_clear(p);
    }
}

void pairing_map(Element g, const EC_POINT P, const EC_POINT Q, const EC_PAIRING p)
{
    p->pairing(g, Q, P, p);
}

void pairing_double_map(Element g, const EC_POINT P1, const EC_POINT Q1, const EC_POINT P2, const EC_POINT Q2, const EC_PAIRING p)
{
    p->pairing_double(g, Q1, P1, Q2, P2, p);
}

const mpz_t* pairing_get_order(const EC_PAIRING p)
{
    return &(p->g1->order);
}

const char* pairing_get_name(const EC_PAIRING p)
{
    return p->pairing_name;
}
