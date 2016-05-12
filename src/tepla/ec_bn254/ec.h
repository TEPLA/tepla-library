#include <tepla/ec.h>

//----------------------------------------------
//  finite field : init clear
//----------------------------------------------
void ec_bn254_fpa_new(Field f);
void ec_bn254_fp2a_new(Field f);
void ec_bn254_fp6a_new(Field f);
void ec_bn254_fp12a_new(Field f);

void ec_bn254_fpb_new(Field f);
void ec_bn254_fp2b_new(Field f);
void ec_bn254_fp6b_new(Field f);
void ec_bn254_fp12b_new(Field f);

void ec_bn254_field_clear(Field f);

//----------------------------------------------
//  ellptic curve : init clear
//----------------------------------------------
void ec_bn254_fpa_group_new(EC_GROUP ec);
void ec_bn254_twa_group_new(EC_GROUP ec);

void ec_bn254_fpb_group_new(EC_GROUP ec);
void ec_bn254_twb_group_new(EC_GROUP ec);

void ec_bn254_group_clear(EC_GROUP ec);

//----------------------------------------------
//  pairing : init clear
//----------------------------------------------
void ec_bn254_pairing_a_new(EC_PAIRING p);
void ec_bn254_pairing_b_new(EC_PAIRING p);

void ec_bn254_pairing_clear(EC_PAIRING p);
