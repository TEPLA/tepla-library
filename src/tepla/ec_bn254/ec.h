#include <tepla/ec.h>

//----------------------------------------------
//  finite field : init clear
//----------------------------------------------
void ec_bn254_fp_new(Field f);
void ec_bn254_fp2_new(Field f);
void ec_bn254_fp6_new(Field f);
void ec_bn254_fp12_new(Field f);

void ec_bn254_field_clear(Field f);

//----------------------------------------------
//  ellptic curve : init clear
//----------------------------------------------
void ec_bn254_fp_group_new(EC_GROUP ec);
void ec_bn254_tw_group_new(EC_GROUP ec);

void ec_bn254_group_clear(EC_GROUP ec);

//----------------------------------------------
//  pairing : init clear
//----------------------------------------------
void ec_bn254_pairing_new(EC_PAIRING p);

void ec_bn254_pairing_clear(EC_PAIRING p);
