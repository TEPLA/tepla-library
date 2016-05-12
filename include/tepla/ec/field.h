//===================================================
//  Elliptic curve operation
//
//   2015.xx.xx
//===================================================
#pragma once

#include <gmp.h>

#ifdef __cplusplus
extern "C" {
#endif

//---------------------------------------------------
// Field type
//---------------------------------------------------
typedef enum fieldtype
{
    Field_fp,   // prime field
    Field_fpn,  // extension field

} FieldType;

//---------------------------------------------------
// Element structure
//---------------------------------------------------
typedef struct ec_element_st
{
    const struct ec_field_st *field;
    void* data;

} Element[1];

//---------------------------------------------------
// Field structure
//---------------------------------------------------
typedef struct ec_field_st
{
    FieldType type;

    char* field_name;

    int ID;

    int str_len;
    int oct_len;

    void (*field_init)(struct ec_field_st *f);
    void (*field_clear)(struct ec_field_st *f);

    void (*init)(Element x);
    void (*clear)(Element x);

    void (*set)(Element x, const Element y);
    void (*set_str)(Element x, const char *str);
    void (*get_str)(char *str, const Element x);
    void (*set_zero)(Element x);
    void (*set_one)(Element x);

    void (*add)(Element z, const Element x, const Element y);
    void (*neg)(Element z, const Element x);
    void (*sub)(Element z, const Element x, const Element y);
    void (*mul)(Element z, const Element x, const Element y);
    void (*sqr)(Element z, const Element x);
    void (*inv)(Element z, const Element x);
    void (*pow)(Element z, const Element x, const mpz_t exp);
    int (*sqrt)(Element z, const Element x);

    int (*is_zero)(const Element x);
    int (*is_one)(const Element x);
    int (*is_sqr)(const Element x);
    int (*cmp)(const Element x, const Element y);

    void (*random)(Element x);

    void (*to_oct)(unsigned char *os, size_t *size, const Element x);
    void (*from_oct)(Element z, const unsigned char *os, const size_t size);

    mpz_t order; // order of field (prime or prime power)

    mpz_t OP1_1;
    mpz_t OP1_2;
    mpz_t OP2;

    int irre_poly_deg; // degree of irreducible polynomial
    int irre_poly_num; // number of irreducible polynomial element

    Element *irre_poly;	// represent irreducible polynomial

    void *precomp; // values of precomputation

    Element *tmp;  // temporary elements for computation

    struct ec_field_st *base; // pointer of base field

} Field[1];

//---------------------------------------------------
//  functions for Field
//---------------------------------------------------
void field_init(Field f, const char *param);
void field_clear(Field f);

const char* field_get_name(const Field f);
const mpz_t* field_get_char(const Field f);

int field_get_degree(const Field f);

//---------------------------------------------------
//  functions for element in Field
//---------------------------------------------------
void element_init(Element x, const Field f);
void element_clear(Element x);

void element_set(Element x, const Element y);
void element_set_str(Element x, const char *str);
void element_get_str(char *str, const Element x);
void element_set_zero(Element x);
void element_set_one(Element x);
void element_add(Element z, const Element x, const Element y);
void element_neg(Element z, const Element x);
void element_sub(Element z, const Element x, const Element y);
void element_mul(Element z, const Element x, const Element y);
void element_sqr(Element z, const Element x);
void element_inv(Element z, const Element x);
void element_pow(Element z, const Element x, const mpz_t exp);
int  element_sqrt(Element z, const Element x);
int  element_is_zero(const Element x);
int  element_is_one(const Element x);
int  element_is_sqr(const Element x);
int  element_cmp(const Element x, const Element y);
void element_random(Element x);
void element_to_oct(unsigned char *os, size_t *size, Element x);
void element_from_oct(Element z, const unsigned char *os, size_t size);

int  element_get_str_length(const Element x);
int  element_get_oct_length(const Element x);

void element_print(const Element x);

#ifdef __cplusplus
}
#endif
