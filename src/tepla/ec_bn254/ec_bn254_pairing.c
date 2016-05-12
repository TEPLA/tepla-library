//=============================================================================
// pairing over BN curve (ec_bn254_pairing) implementation with GMP
//=============================================================================

#include "ec_bn254_lcl.h"

#define xcoord(p)   (p->x)
#define ycoord(p)   (p->y)
#define zcoord(p)   (p->z)

#define field(p)   (p->ec->field)
#define curve(p)   (p->ec)

//-------------------------------------------
//  precomputation for pairing
//-------------------------------------------
void ec_bn254_pairing_precomp(EC_PAIRING p)
{
	pairing_precomp_p precomp;

	int s[] = {0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,-1,-1,0,0,0,0,0,0,1,1};
	int t[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,1};

	int *sbuff, *tbuff;

	precomp = (pairing_precomp_p)malloc(sizeof(struct ec_pairing_precomp_st));

	sbuff = (int *)malloc(sizeof(s));
	memcpy(sbuff, s, sizeof(s));

	precomp->si = sbuff;
	precomp->slen = sizeof(s)/sizeof(*s);

	tbuff = (int *)malloc(sizeof(t));
	memcpy(tbuff, t, sizeof(t));

	precomp->ti = tbuff;
	precomp->tlen = sizeof(t)/sizeof(*t);

	p->precomp = (void*)precomp;
}

//-------------------------------------------
//  pairing
//-------------------------------------------
void ec_bn254_pairing_dob(EC_POINT T, Element l0, Element l3, Element l4, const EC_POINT P)
{
	Element *t = field(T)->tmp;

	bn254_fp2_sqr(t[4], zcoord(T)); //ZT_2 = ZT^2
	bn254_fp2_sqr(l0, xcoord(T));   //l0 = XT^2

	bn254_fp2_sqr(t[0], ycoord(T)); //tmp0 = YT^0
	bn254_fp2_sqr(t[1], t[0]);      //tmp1 = tmp0^2

	bn254_fp2_add(l3, xcoord(T), t[0]); //l3 = XT + tmp0
	bn254_fp2_sqr(l3, l3);       //l3 = l3^3
	bn254_fp2_sub(l3, l3, l0);   //l3 = l3 - l0
	bn254_fp2_sub(l3, l3, t[1]); //l3 = l3 - tmp1
	bn254_fp2_add(l3, l3, l3);   //l3 = 2*l3

	bn254_fp2_add(t[2], l0, l0);   //tmp2 = 3*l0
	bn254_fp2_add(t[2], t[2], l0); //tmp2 = 3*l0
	bn254_fp2_sqr(t[3], t[2]);     //tmp3 = tmp2^2

	bn254_fp2_add(l4, t[2], xcoord(T)); //l4 = tmp2 + XT
	bn254_fp2_sqr(l4, l4);       //l4 = l4^2
	bn254_fp2_sub(l4, l4, l0);   //l4 = l4 - l0
	bn254_fp2_sub(l4, l4, t[3]); //l4 = l4 - tmp3
	bn254_fp2_sub(xcoord(T), t[3], l3); //XT = tmp3 - l3
	bn254_fp2_sub(xcoord(T), xcoord(T), l3); //XT = XT - l3

	bn254_fp2_add(zcoord(T), ycoord(T), zcoord(T)); //ZT = YT + ZT
	bn254_fp2_sqr(zcoord(T), zcoord(T));       //ZT = ZT^2
	bn254_fp2_sub(zcoord(T), zcoord(T), t[0]); //ZT = ZT - tmp0
	bn254_fp2_sub(zcoord(T), zcoord(T), t[4]); //ZT = ZT - ZT_2

	bn254_fp2_sub(ycoord(T), l3, xcoord(T));   //YT = l3 - XT
	bn254_fp2_mul(ycoord(T), ycoord(T), t[2]); //YT = YT*tmp2
	bn254_fp2_add(t[1], t[1], t[1]);   //tmp1 = 2*tmp1
	bn254_fp2_add(t[1], t[1], t[1]);   //tmp1 = 2*tmp1
	bn254_fp2_add(t[1], t[1], t[1]);   //tmp1 = 2*tmp1
	bn254_fp2_sub(ycoord(T), ycoord(T), t[1]); //YT = YT - tmp1

	bn254_fp2_mul(l3, t[2], t[4]);      //l3 = tmp2*ZT
	bn254_fp2_mul(l0, zcoord(T), t[4]); //l0 = ZT*ZT_2

	bn254_fp2_neg(l3, l3);   //l3 = -l3
	bn254_fp2_mul_p(l3, l3, xcoord(P)); //l3 = l3*XP

	bn254_fp2_add(t[1], t[0], t[0]); //tmp1 = 2*tmp0
	bn254_fp2_add(t[1], t[1], t[1]); //tmp1 = 2*tmp1
	bn254_fp2_sub(l4, l4, t[1]);     //l4 = l4 - tmp1

	bn254_fp2_mul_p(l0, l0, ycoord(P)); //l0 = l0*YP

	bn254_fp2_add(l0, l0, l0);
	bn254_fp2_add(l3, l3, l3);
}

void ec_bn254_pairing_add(EC_POINT T, Element l0, Element l3, Element l4, const EC_POINT Q, const EC_POINT P)
{
	Element *t = field(T)->tmp;

	bn254_fp2_sqr(t[5], zcoord(T));   //ZT2 = ZT^2
	bn254_fp2_sqr(t[6], ycoord(Q));   //YQ2 = YQ^2

	bn254_fp2_mul(t[0], xcoord(Q), t[5]);    //t0 = XQ*ZT^2
	bn254_fp2_add(l3, ycoord(Q), zcoord(T)); //l3 = YQ*ZT
	bn254_fp2_sqr(l3, l3);         //l3 = l3^2
	bn254_fp2_sub(l3, l3, t[6]);   //l3 = l3 - YQ^2
	bn254_fp2_sub(l3, l3, t[5]);   //l3 = l3 - ZT^2
	bn254_fp2_mul(l3, l3, t[5]);   //l3 = l3*ZT^2

	bn254_fp2_sub(t[0], t[0], xcoord(T)); //tmp0 = tmp0 - XT
	bn254_fp2_sqr(t[1], t[0]);            //tmp1 = tmp0^2

	bn254_fp2_add(t[2], t[1], t[1]); //tmp2 = 2*tmp1
	bn254_fp2_add(t[2], t[2], t[2]); //tmp2 = 2*tmp2

	bn254_fp2_mul(t[3], t[0], t[2]); //tmp3 = tmp0*tmp2

	bn254_fp2_sub(t[4], l3, ycoord(T));   //tmp4 = l3 - YT
	bn254_fp2_sub(t[4], ycoord(T), t[4]); //tmp4 = YT - tmp4

	bn254_fp2_mul(l4, t[4], xcoord(Q));   //l4 = tmp4*XQ
	bn254_fp2_mul(t[2], t[2], xcoord(T)); //tmp2 = tmp2*XT
	bn254_fp2_sqr(xcoord(T), t[4]);       //XT = tmp4^2
	bn254_fp2_sub(xcoord(T), xcoord(T), t[3]); //XT = XT - t3
	bn254_fp2_sub(xcoord(T), xcoord(T), t[2]); //XT = XT - t2
	bn254_fp2_sub(xcoord(T), xcoord(T), t[2]); //XT = XT - t2

	bn254_fp2_add(zcoord(T), zcoord(T), t[0]); //ZT = ZT + tmp0
	bn254_fp2_sqr(zcoord(T), zcoord(T));       //ZT = ZT^2
	bn254_fp2_sub(zcoord(T), zcoord(T), t[5]); //ZT = ZT - ZT^2
	bn254_fp2_sub(zcoord(T), zcoord(T), t[1]); //ZT = ZT - tmp1

	bn254_fp2_sqr(t[5], zcoord(T));
	bn254_fp2_add(t[1], ycoord(Q), zcoord(T)); //tmp1 = YQ + ZT
	bn254_fp2_sub(t[2], xcoord(T), t[2]);      //tmp2 = XT - tmp2
	bn254_fp2_mul(t[2], t[2], t[4]);           //tmp2 = tmp2*tmp4

	bn254_fp2_mul(t[0], ycoord(T), t[3]); //tmp0 = YT*tmp3
	bn254_fp2_add(t[0], t[0], t[0]);      //tmp0 = 2*tmp0
	bn254_fp2_sub(ycoord(T), t[2], t[0]); //YT = tmp2 - tmp0

	bn254_fp2_sqr(t[1], t[1]);       //tmp1 = tmp1^2
	bn254_fp2_sub(t[1], t[1], t[6]); //tmp1 = tmp1 - YQ^2
	bn254_fp2_sub(t[1], t[5], t[1]); //tmp1 = ZT^2 - tmp1

	bn254_fp2_add(l4, l4, l4);   //l4 = 2*l4
	bn254_fp2_sub(l4, t[1], l4); //l4 = tmp1 - l4

	bn254_fp2_mul_p(l0, zcoord(T), ycoord(P)); // l0 = ZT * YP
	bn254_fp2_mul_p(l3, t[4], xcoord(P));      // l3 = tmp4 * XP

	bn254_fp2_add(l0, l0, l0);
	bn254_fp2_add(l3, l3, l3);
}

void ec_bn254_pairing_miller(Element z, const EC_POINT Q, const EC_POINT P, const EC_PAIRING p)
{
	int len, *s, i;

	EC_POINT T, R, S;

	Element xq, yq;
	Element f, l0, l3, l4;

	//--------------------------------
	//   init
	//--------------------------------
	element_init(f, z->field);
	element_init(xq, field(Q));
	element_init(yq, field(Q));

	point_init(R, curve(Q));
	point_init(T, curve(Q));
	point_init(S, curve(Q));

	element_init(l0, field(Q));
	element_init(l3, field(Q));
	element_init(l4, field(Q));

	len = ((pairing_precomp_p)(p->precomp))->slen;  // s = PAIRING->precomp->si
	s = ((pairing_precomp_p)(p->precomp))->si;

	bn254_fp12_set_one(f);        // f = 1
	ec_bn254_fp2_point_set(T, Q); // T = Q
	ec_bn254_fp2_neg(R, Q);       // R = -Q

	//-------------------------------
	//  Miller loop
	//-------------------------------
	for (i =len-2; i >= 0; i--)
	{
		ec_bn254_pairing_dob(T, l0, l3, l4, P);   //T = 2T, l = l(P)

		bn254_fp12_sqr(f, f);             // f = f^2*l
		bn254_fp12_mul_L(f, l0, l3, l4);  //

		if ( s[i] )
		{
			if ( s[i] < 0 )
			{
				ec_bn254_pairing_add(T, l0, l3, l4, R, P);   // T = T - Q, l = l(P)
			}
			else
			{
				ec_bn254_pairing_add(T, l0, l3, l4, Q, P);   // T = T + Q, l = l(P)
			}

			bn254_fp12_mul_L(f, l0, l3, l4);
		}
	}

	//--------------------------------
	//   addition part
	//--------------------------------
	ec_bn254_tw_frob(S, Q);
	ec_bn254_pairing_add(T, l0, l3, l4, S, P);   //addtion part 2
	bn254_fp12_mul_L(f, l0, l3, l4);

	ec_bn254_tw_frob2(S, Q);
	ec_bn254_fp2_neg(S, S);
	ec_bn254_pairing_add(T, l0, l3, l4, S, P);   //addtion part 2
	bn254_fp12_mul_L(f, l0, l3, l4);

	bn254_fp12_set(z, f);

	//--------------------------------
	//   relase
	//--------------------------------
	element_clear(f);
	element_clear(l0);
	element_clear(l3);
	element_clear(l4);
	element_clear(xq);
	element_clear(yq);

	point_clear(T);
	point_clear(R);
	point_clear(S);
}

void ec_bn254_pairing_finalexp(Element z, const Element x, const EC_PAIRING p)
{
	Element *t = z->field->tmp;

	int len, *u;

	len = ((pairing_precomp_p)(p->precomp))->tlen;
	u = ((pairing_precomp_p)(p->precomp))->ti;

	bn254_fp12_conj(t[0], x);      // f1 := conjugate of x
	bn254_fp12_inv(t[1], x);       // f2 := x^-1
	bn254_fp12_mul(z, t[0], t[1]); // z := f1*f2
	bn254_fp12_frob_p2(t[0], z);   // ftmp = z^(p^2)
	bn254_fp12_mul(z, z, t[0]);    // z := z*ftmp
	bn254_fp12_pow_forpairing(t[0], z, u, len);
	bn254_fp12_pow_forpairing(t[1], t[0], u, len);
	bn254_fp12_pow_forpairing(t[2], t[1], u, len);
	bn254_fp12_frob_p(t[3], z);       // fp1 = z^p
	bn254_fp12_frob_p2(t[4], z);      // fp2 = z^(p^2)
	bn254_fp12_frob_p3(t[5], z);      // fp3 = z^(p^3)
	bn254_fp12_mul(t[3], t[3], t[4]); // y0 = fp1*fp2
	bn254_fp12_mul(t[3], t[3], t[5]); // y0 = y0*fp3;
	bn254_fp12_conj(t[4], z);         // y1 = conjugate of z
	bn254_fp12_frob_p2(t[5], t[1]);   // y2 = ft2^(p^2)
	bn254_fp12_frob_p(t[6], t[0]);    // y3 = ft1^p
	bn254_fp12_conj(t[6], t[6]);      // y3 = conjugate of y3
	bn254_fp12_frob_p(t[7], t[1]);    // y4 = ft2^p
	bn254_fp12_mul(t[7], t[7], t[0]); // y4 = y4*ft1
	bn254_fp12_conj(t[7], t[7]);      // y4 = conjugate of y4
	bn254_fp12_conj(t[8], t[1]);      // y5 = conjugate of ft2
	bn254_fp12_frob_p(t[9], t[2]);    // y6 = ft3^p
	bn254_fp12_mul(t[9], t[9], t[2]); // y6 = y6*ft3
	bn254_fp12_conj(t[9], t[9]);      // y6 = conjugate of y6
	bn254_fp12_sqr(t[0], t[9]);       // t0 = y6^2
	bn254_fp12_mul(t[0], t[0], t[7]); // t0 = t0*y4
	bn254_fp12_mul(t[0], t[0], t[8]); // t0 = t0*y5
	bn254_fp12_mul(t[1], t[6], t[8]); // t1 = y3*y5
	bn254_fp12_mul(t[1], t[1], t[0]); // t1 = t1*t0
	bn254_fp12_mul(t[0], t[0], t[5]); // t0 = t0*y2
	bn254_fp12_sqr(t[1], t[1]);       // t1 = t1^2
	bn254_fp12_mul(t[1], t[1], t[0]); // t1 = t1*t0
	bn254_fp12_sqr(t[1], t[1]);       // t1 = t1^2
	bn254_fp12_mul(t[0], t[1], t[4]); // t0 = t1*y1
	bn254_fp12_mul(t[1], t[1], t[3]); // t1 = t1*y0
	bn254_fp12_sqr(t[0], t[0]);       // t0 = t0^2
	bn254_fp12_mul(z, t[1], t[0]);    // z = t1*t0
}

void ec_bn254_pairing(Element z, const EC_POINT Q, const EC_POINT P, const EC_PAIRING p)
{
	ec_bn254_pairing_miller(z, Q, P, p);
	ec_bn254_pairing_finalexp(z, z, p);
}
