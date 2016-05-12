#pragma once

#include <gmp.h>

#ifdef __cplusplus
extern "C" {
#endif

//------------------------------------------------------------------
//  The Function IHF1
//------------------------------------------------------------------
//   Input :
//      A string s in {0,1}*
//      A size of s
//      An integer n
//      A security parameter t in {80, 112, 128, 192, 256}
//      -- n <= 2^{2t}
//   Output :
//      An integer v in Zn
//------------------------------------------------------------------
void IHF1_SHA(mpz_t v, const unsigned char *s, size_t slen, const mpz_t n, int t);

//------------------------------------------------------------------
//  The Function mIHF
//------------------------------------------------------------------
//   Input :
//      A string s in {0,1}*
//      A size of s
//      A security parameter t in {80, 112, 128, 192, 256}
//   Output :
//      A digest d from s
//      A size of d
//	A size of d
//------------------------------------------------------------------
void mIHF_SHA(unsigned char *d, size_t *dlen, const char *s, size_t slen, int t);

#ifdef __cplusplus
}
#endif
