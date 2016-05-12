//-----------------------------------------------------------------
// Interface of Hash Function using TEPLA
//-----------------------------------------------------------------
#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>

#include <openssl/sha.h>

typedef unsigned char *(*Hash)(const unsigned char *s, size_t n, unsigned char *md);

//---------------------------------------------------------------------
//  Select Hash Function with security parameter t
//---------------------------------------------------------------------
void selectHash(Hash *H, int *size, int t)
{
	switch (t)
	{
	case 80:
		(*H) = SHA1;
		(*size) = SHA_DIGEST_LENGTH;
		break;
	case 112:
		(*H) = SHA224;
		(*size) = SHA224_DIGEST_LENGTH;
		break;
	case 128:
		(*H) = SHA256;
		(*size) = SHA256_DIGEST_LENGTH;
		break;
	case 192:
		(*H) = SHA384;
		(*size) = SHA384_DIGEST_LENGTH;
		break;
	case 256:
		(*H) = SHA512;
		(*size) = SHA512_DIGEST_LENGTH;
		break;
	default:
		fprintf(stderr, "We do not support the parameter t = %d\n", t);
		exit(100);
	}
}

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
void IHF1_SHA(mpz_t v, const unsigned char *s, size_t slen, const mpz_t n, int t)
{
	mpz_t v0, v1, v2, a1, a2;

	Hash H = NULL;

	int size;

	unsigned char *t1, *t2;
	unsigned char *h1, *h2;

	mpz_init(v0);
	mpz_init(v1);
	mpz_init(v2);
	mpz_init(a1);
	mpz_init(a2);

	selectHash(&H, &size, t);

	t1 = (unsigned char *)malloc(sizeof(unsigned char)*(size+(int)slen));
	t2 = (unsigned char *)malloc(sizeof(unsigned char)*(size+(int)slen));

	h1 = (unsigned char *)malloc(sizeof(unsigned char)*(size));
	h2 = (unsigned char *)malloc(sizeof(unsigned char)*(size));

	mpz_set_ui(v0, 0);

	//-----------------------------
	//  t1 = {0}^{2*t} || s
	//-----------------------------
	memset(t1, 0x00, size);
	memcpy(&(t1[size]), s, slen);

	//-----------------------------
	//  h1 = Hash(t1)
	//-----------------------------
	(*H)(t1, size+(int)slen, h1);

	//-----------------------------
	//  a1 = OS2IP(h1)
	//-----------------------------
	mpz_import(a1, size, 1, sizeof(*h1), 0, 0, h1);

	//-----------------------------
	//  t2 = h1 || s
	//-----------------------------
	memcpy(t2, h1, size);
	memcpy(&(t2[size]), s, slen);

	//-----------------------------
	//  h2 = Hash(t2)
	//-----------------------------
	(*H)(t2, size+(int)slen, h2);

	//-----------------------------
	//  a2 = OS2IP(h2)
	//-----------------------------
	mpz_import(a2, size, 1, sizeof(*h2), 0, 0, h2);

	//-----------------------------
	//  v2 = 2^{2t}*a1 + a2
	//-----------------------------
	mpz_mul_2exp(v1, a1, 2*t);
	mpz_add(v2, v1, a2);

	//-----------------------------
	//  Output: v = v2 mod n
	//-----------------------------
	mpz_mod(v, v2, n);

	mpz_clear(v0);
	mpz_clear(v1);
	mpz_clear(v2);
	mpz_clear(a1);
	mpz_clear(a2);

	free(t1);
	free(t2);
	free(h1);
	free(h2);
}

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
void mIHF_SHA(unsigned char *d, size_t *dlen, const char *s, size_t slen ,int t)
{
	Hash H = NULL;

	int len;

	selectHash(&H, &len, t);

	(*H)((const unsigned char *)s, slen, d);

	(*dlen) = (size_t)len;

	return;
}
