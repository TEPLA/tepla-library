//============================================================================
//  elliptic curve over extension field ( bn254_fp2 ) implementation with GMP
//----------------------------------------------------------------------------
//  ec_bn254_twa tw_EC := Y^2 = X^3+5/i Beuchat et al.
//  ec_bn254_twb tw_EC := Y^2 = X^3+2/i Aranha et al.
//----------------------------------------------------------------------------
//  2015.10.31 created by kanbara
//============================================================================

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <tepla/hash.h>

#include "ec_bn254_lcl.h"

#define MAX(a, b) (((a)>(b))? (a): (b) )

#define mpz_rep(x) (*((mpz_t*)((x)->data)))
#define elt_rep(x, i) (mpz_rep(((Element*)((x)->data))[(i)]))

#define xcoord(p)  (p->x)
#define ycoord(p)  (p->y)
#define zcoord(p)  (p->z)

#define field(p)   (p->ec->field)
#define curve(p)   (p->ec)

//-------------------------------------------
//  initialization, clear, set
//-------------------------------------------
void ec_bn254_fp2_point_init(EC_POINT P)
{
    element_init(xcoord(P), field(P));
    element_init(ycoord(P), field(P));
    element_init(zcoord(P), field(P));

    point_set_infinity(P);
}

void ec_bn254_fp2_point_clear(EC_POINT P)
{
    element_clear(xcoord(P));
    element_clear(ycoord(P));
    element_clear(zcoord(P));

    P->ec = NULL;
}

void ec_bn254_fp2_point_set(EC_POINT z, const EC_POINT x)
{
    element_set(xcoord(z), xcoord(x));
    element_set(ycoord(z), ycoord(x));
    element_set(zcoord(z), zcoord(x));

    z->isinfinity = x->isinfinity;
}

void ec_bn254_fp2_point_set_str(EC_POINT P, const char* s)
{
    static const char infinity[] = "[0]";

    int i, hr;
    int len = strlen(s);

    char *tmp, *p1, *p2;

    if (strcmp(s, infinity) == 0) {
        point_set_infinity(P);
        return;
    }

    tmp = (char*)malloc(sizeof(char) * (len + 1));

    strcpy(tmp, s);

    hr = FALSE;
    p1 = NULL;
    p2 = NULL;

    for (i = 0; i < len; i++)
    {
        if (tmp[i] == '[') {
            p1 = &(tmp[i + 1]);
        }
        if (tmp[i] == ',') {
            p2 = &(tmp[i + 1]);
            tmp[i] = '\0';
        }
        if (tmp[i] == ']') {
            tmp[i] = '\0';
            hr = (p1 != NULL && p2 != NULL);
        }
    }

    if (!hr)
    {
        fprintf(stderr, "Please check input string : format = [x,y] or [0]\n");
        exit(200);
    }

    bn254_fp2_set_str(xcoord(P), p1);
    bn254_fp2_set_str(ycoord(P), p2);
    bn254_fp2_set_one(zcoord(P));

    P->isinfinity = FALSE;

    free(tmp);
}

void ec_bn254_fp2_point_get_str(char *s, const EC_POINT P)
{
    static const char infinity[] = "[0]";

    char sx[130], sy[130];

    EC_POINT T;

    if (point_is_infinity(P)) {
        strcpy(s, infinity);
        return;
    }

    point_init(T, curve(P));

    point_make_affine(T, P);

    bn254_fp2_get_str(sx, xcoord(T));
    bn254_fp2_get_str(sy, ycoord(T));

    sprintf(s, "[%s,%s]", sx, sy);

    point_clear(T);
}

void ec_bn254_fp2_point_set_xy(EC_POINT z, const Element x, const Element y)
{
    element_set(xcoord(z), x);
    element_set(ycoord(z), y);
    element_set_one(zcoord(z));

    z->isinfinity = FALSE;
}

void ec_bn254_fp2_point_set_infinity(EC_POINT P)
{
    element_set_one(xcoord(P));
    element_set_one(ycoord(P));
    element_set_zero(zcoord(P));

    P->isinfinity = TRUE;
}

//-------------------------------------------
//  arithmetic operation
//-------------------------------------------
void ec_bn254_fp2_add(EC_POINT R, const EC_POINT P, const EC_POINT Q)
{
    Element *t = field(R)->tmp;

    if (point_is_infinity(P)) {
        point_set(R, Q);
        return;
    }
    if (point_is_infinity(Q)) {
        point_set(R, P);
        return;
    }

    if (bn254_fp2_cmp(xcoord(P), xcoord(Q)) == 0)
    {
        if (bn254_fp2_cmp(ycoord(P), ycoord(Q)) == 0) {
            ec_bn254_fp2_dob(R, P);
            return;
        }
        point_set_infinity(R);
        return;
    }

    bn254_fp2_sub(t[0], xcoord(Q), xcoord(P));
    bn254_fp2_sub(t[1], ycoord(Q), ycoord(P));
    bn254_fp2_inv(t[0], t[0]);
    bn254_fp2_mul(t[2], t[1], t[0]);

    bn254_fp2_sqr(t[0], t[2]);
    bn254_fp2_sub(t[0], t[0], xcoord(P));
    bn254_fp2_sub(t[0], t[0], xcoord(Q));

    bn254_fp2_sub(t[1], xcoord(P), t[0]);
    bn254_fp2_mul(t[1], t[2], t[1]);
    bn254_fp2_sub(ycoord(R), t[1], ycoord(P));
    bn254_fp2_set(xcoord(R), t[0]);
    bn254_fp2_set_one(zcoord(R));

    R->isinfinity = FALSE;
}

void ec_bn254_fp2_dob(EC_POINT R, const EC_POINT P)
{
    Element *t = field(R)->tmp;

    if (point_is_infinity(P)) {
        point_set_infinity(R);
        return;
    }
    if (bn254_fp2_is_zero(ycoord(P))) {
        point_set_infinity(R);
        return;
    }

    bn254_fp2_add(t[0], ycoord(P), ycoord(P));
    bn254_fp2_inv(t[0], t[0]);
    bn254_fp2_sqr(t[1], xcoord(P));
    bn254_fp2_add(t[2], t[1], t[1]);
    bn254_fp2_add(t[1], t[1], t[2]);
    bn254_fp2_mul(t[2], t[1], t[0]);

    bn254_fp2_sqr(t[0], t[2]);
    bn254_fp2_sub(t[0], t[0], xcoord(P));
    bn254_fp2_sub(t[0], t[0], xcoord(P));

    bn254_fp2_sub(t[1], xcoord(P), t[0]);
    bn254_fp2_mul(t[1], t[2], t[1]);
    bn254_fp2_sub(ycoord(R), t[1], ycoord(P));
    bn254_fp2_set(xcoord(R), t[0]);
    bn254_fp2_set_one(zcoord(R));

    R->isinfinity = FALSE;
}

void ec_bn254_fp2_neg(EC_POINT Q, const EC_POINT P)
{
    if (point_is_infinity(P)) {
        point_set_infinity(Q);
        return;
    }

    bn254_fp2_set(xcoord(Q), xcoord(P));
    bn254_fp2_neg(ycoord(Q), ycoord(P));
    bn254_fp2_set(zcoord(Q), zcoord(P));

    Q->isinfinity = P->isinfinity;
}

void ec_bn254_fp2_sub(EC_POINT R, const EC_POINT P, const EC_POINT Q)
{
    EC_POINT T;

    point_init(T, curve(R));

    ec_bn254_fp2_neg(T, Q);
    ec_bn254_fp2_add(R, P, T);

    point_clear(T);
}

void ec_bn254_fp2_add_formul(EC_POINT R, const EC_POINT P, const EC_POINT Q)
{
    Element *t = field(R)->tmp;

    if (point_is_infinity(P)) {
        point_set(R, Q);
        return;
    }
    if (point_is_infinity(Q)) {
        point_set(R, P);
        return;
    }

    if (strcmp(P->ec->curve_name, "ec_bn254_twa") == 0)
    {
        bn254_fp2_sqr(t[0], zcoord(P));       // A = Pz^2
        bn254_fp2_mul(t[1], t[0], zcoord(P)); // B = Pz^3
        bn254_fp2_mul(t[2], t[0], xcoord(Q)); // C = Qx*A
        bn254_fp2_mul(t[1], t[1], ycoord(Q)); // D = Qy*B
        bn254_fp2_sub(t[2], t[2], xcoord(P)); // E = C-Px
        bn254_fp2_sub(t[1], t[1], ycoord(P)); // F = D-Py

        if (bn254_fp2_is_zero(t[2]))
        {
            if (bn254_fp2_is_zero(t[1])) {
                ec_bn254_fp2_dob_formul(R, P);
                return;
            }
            point_set_infinity(R);
            return;
        }

        bn254_fp2_mul(zcoord(R), zcoord(P), t[2]); // Rz = Pz*E

        bn254_fp2_sqr(t[0], t[2]);            // G = E^2
        bn254_fp2_mul(t[3], xcoord(P), t[0]); // I = Px*G
        bn254_fp2_mul(t[4], t[0], t[2]);      // H = G*E

        bn254_fp2_sqr(xcoord(R), t[1]);
        bn254_fp2_sub(xcoord(R), xcoord(R), t[4]);
        bn254_fp2_add(t[0], t[3], t[3]);
        bn254_fp2_sub(xcoord(R), xcoord(R), t[0]); // Rx = F^2-H-2*I-Rx

        bn254_fp2_sub(t[3], t[3], xcoord(R));
        bn254_fp2_mul(t[0], t[4], ycoord(P));
        bn254_fp2_mul(ycoord(R), t[1], t[3]);
        bn254_fp2_sub(ycoord(R), ycoord(R), t[0]); // Ry = F*(I-Rx)-Py*H
    }

    // Jacobian coordinate proposed by Aranha et al.
    if (strcmp(P->ec->curve_name, "ec_bn254_twb") == 0)
    {
        if (point_cmp(P, Q) == 0) {
            ec_bn254_fp2_dob_formul(R, P);
            return;
        }

        bn254_fp2_sqr(t[1], zcoord(P));        		// t1 = Pz^2
        bn254_fp2_mul(t[3], xcoord(Q), t[1]); 		// t3 = Qx*t1
        bn254_fp2_mul(t[1], t[1], zcoord(P));  		// t1 = t1*Pz
        bn254_fp2_sub(t[3], t[3], xcoord(P));  		// t3 = t3-Px
        bn254_fp2_mul(t[4], t[1], ycoord(Q));  		// t4 = t1*Qy
        bn254_fp2_mul(zcoord(R), zcoord(P), t[3]); 	// Rz = Pz*t3

        if (bn254_fp2_is_zero(zcoord(R)))
        {
            point_set_infinity(R);
            return;
        }

        bn254_fp2_sub(t[0], t[4], ycoord(P));  		// t0 = t4-Py
        bn254_fp2_sqr(t[1], t[3]);             		// t1 = t3^2
        bn254_fp2_mul(t[4], t[1], t[3]);	   		// t4 = t1*t3
        bn254_fp2_mul(t[1], t[1], xcoord(P));  		// t1 = t1*Px
        bn254_fp2_sqr(xcoord(R), t[0]);		   		// Rx = t0^2
        bn254_fp2_dob(t[3], t[1]);			   		// t3 = 2*t1
        bn254_fp2_sub(xcoord(R), xcoord(R), t[3]); 	// Rx = Rx-t3
        bn254_fp2_sub(xcoord(R), xcoord(R), t[4]); 	// Rx = Rx-t4
        bn254_fp2_sub(t[1], t[1], xcoord(R));  		// t1 = t1-Rx
        bn254_fp2_mul(t[2], t[0], t[1]); 			// t2 = t0*t1
        bn254_fp2_mul(t[3], t[4], ycoord(P)); 		// t3 = t4*Py
        bn254_fp2_sub(t[2], t[2], t[3]); 			// t2 = t2-t3
        bn254_fp2_mod(ycoord(R), t[2]);  			// Ry = t2 mod p
    }
    R->isinfinity = FALSE;
}

void ec_bn254_fp2_dob_formul(EC_POINT R, const EC_POINT P)
{
    Element *t = field(R)->tmp;

    if (point_is_infinity(P)) {
        point_set_infinity(R);
        return;
    }

    if (strcmp(P->ec->curve_name, "ec_bn254_twa") == 0)
    {
        bn254_fp2_sqr(t[0], ycoord(P));  // A = Py^2
        bn254_fp2_add(t[1], xcoord(P), xcoord(P));
        bn254_fp2_add(t[1], t[1], t[1]);
        bn254_fp2_mul(t[1], t[1], t[0]); // B = 4*Px*A
        bn254_fp2_sqr(t[0], t[0]);
        bn254_fp2_add(t[0], t[0], t[0]);
        bn254_fp2_add(t[0], t[0], t[0]);
        bn254_fp2_add(t[0], t[0], t[0]); // C = 8*A^2
        bn254_fp2_sqr(t[3], xcoord(P));
        bn254_fp2_add(t[2], t[3], t[3]);
        bn254_fp2_add(t[2], t[2], t[3]); // D = 3*Px^2

        bn254_fp2_sqr(xcoord(R), t[2]);
        bn254_fp2_add(t[3], t[1], t[1]);
        bn254_fp2_sub(xcoord(R), xcoord(R), t[3]); // Rx = D^2 - 2*B

        bn254_fp2_mul(zcoord(R), ycoord(P), zcoord(P));
        bn254_fp2_add(zcoord(R), zcoord(R), zcoord(P)); // Rz = 2*Py*Pz

        bn254_fp2_sub(ycoord(R), t[1], xcoord(R));
        bn254_fp2_mul(ycoord(R), ycoord(R), t[2]);
        bn254_fp2_sub(ycoord(R), ycoord(R), t[0]);  // Ry = D*(B-Rx)-C
    }

    // Jacobian coordinate proposed by Aranha et al.
    if (strcmp(P->ec->curve_name, "ec_bn254_twb") == 0)
    {
        bn254_fp2_sqr(t[0], xcoord(P)); 			// t0 = Px^2
        bn254_fp2_dob(t[1], t[0]);      			// t1 = 2*t0
        bn254_fp2_mul(zcoord(R), ycoord(P), zcoord(P)); // Rz = Py*Pz
        bn254_fp2_add(t[0], t[0], t[1]); 			// t0 = t0+t1
        bn254_fp2_sqr(t[3], ycoord(P));  			// t3 = Py^2
        bn254_fp2_div_2(t[0], t[0]);				// t0 = t0/2
        bn254_fp2_mul(t[1], t[3], xcoord(P)); 		// t1 = t3*Px
        bn254_fp2_dob(ycoord(R), t[1]); 	  		// Ry = 2*t1
        bn254_fp2_sqr(xcoord(R), t[0]);       		// Rx = t0^2
        bn254_fp2_sub(xcoord(R), xcoord(R), ycoord(R)); // Rx = Rx-Ry
        bn254_fp2_sub(t[1], t[1], xcoord(R)); 		// t1 = t1-Rx
        bn254_fp2_sqr(t[2], t[3]); 					// t2 = t3^2
        bn254_fp2_mul(t[1], t[0], t[1]);			// t1 = t0*t1
        bn254_fp2_sub(t[1], t[1], t[2]); 			// t1 = t1-t2
        bn254_fp2_mod(ycoord(R), t[1]); 			// Ry = t1 mod p
    }

    R->isinfinity = FALSE;
}

//--------------------------------------------------------------
//  Scalar Multiplication in Jacobian Coordinate
//--------------------------------------------------------------
void ec_bn254_fp2_mul(EC_POINT Q, const mpz_t s, const EC_POINT P)
{
    long t, i;
    EC_POINT R;

    point_init(R, curve(P));

    ec_bn254_fp2_point_set(R, P);

    t = mpz_sizeinbase(s, 2);

    for (i = t - 2; i >= 0; i--)
    {
        ec_bn254_fp2_dob_formul(R, R);
        if (mpz_tstbit(s, i))
        {
            ec_bn254_fp2_add_formul(R, R, P);
        }
    }

    point_make_affine(Q, R);

    point_clear(R);
}

//------------------------------------------------------
//  Scalar Multiplication with NAF
//------------------------------------------------------
void ec_bn254_fp2_mul_naf(EC_POINT Q, const mpz_t s, const EC_POINT P)
{
    long t, i;

    int *naf, nlen;

    EC_POINT R, mP;

    point_init(R, curve(P));
    point_init(mP, curve(P));

    ec_bn254_fp2_point_set(R, P);
    ec_bn254_fp2_neg(mP, P);

    t = mpz_sizeinbase(s, 2);

    naf = (int *)malloc(sizeof(int) * (t + 1));

    generate_naf(naf, &nlen, s);

    for (i = nlen - 2; i >= 0; i--)
    {
        ec_bn254_fp2_dob_formul(R, R);
        if (naf[i])
        {
            if (naf[i] < 0) {
                ec_bn254_fp2_add_formul(R, R, mP);
            }
            else {
                ec_bn254_fp2_add_formul(R, R, P);
            }
        }
    }

    point_make_affine(Q, R);

    point_clear(R);
    point_clear(mP);

    free(naf);
}

//----------------------------------------------------------------
//  Elliptic Curve Parameter Setting
//----------------------------------------------------------------
void ec_bn254_fp2_init_ec_data_aranha(EC_GROUP ec)
{
    ec_data_fp2 d = NULL;
    struct ec_field_st *fp2;

    mpz_t x;

    d = (ec_data_fp2)malloc(sizeof(struct ec_bn254_fp2_ec_data_st));

    mpz_init(x);
    mpz_set_str(x, "4080000000000001", 16); // x = -x

    mpz_init(d->_6x);
    mpz_init(d->_6x2);
    mpz_mul_ui(d->_6x, x, 6);
    mpz_mul(d->_6x2, d->_6x, x);

    fp2 = ec->field;

    element_init(d->vfrobx, fp2);
    element_init(d->vfroby, fp2);
    element_init(d->vfrobx2, fp2);
    element_init(d->vfroby2, fp2);
    element_init(d->vfrobx3, fp2);
    element_init(d->vfroby3, fp2);

    element_set_str(d->vfrobx, "0 25236482400000017080EB4000000006181800000000000CD98000000000000B");
    element_set_str(d->vfroby, "23DFC9D1A39F4DB8C69B87A8848AA075A7333A0E62D78CBF4B1B8EEAE58B81C5 23DFC9D1A39F4DB8C69B87A8848AA075A7333A0E62D78CBF4B1B8EEAE58B81C5");
    element_set_str(d->vfrobx2, "49B36240000000024909000000000006CD80000000000007 0");
    element_set_str(d->vfroby2, "2523648240000001BA344D80000000086121000000000013A700000000000012 0");
    element_set_str(d->vfrobx3, "0 1");
    element_set_str(d->vfroby3, "1439AB09C60B248F398C5D77B755F92B9EDC5F19D2873545BE471151A747E4E 1439AB09C60B248F398C5D77B755F92B9EDC5F19D2873545BE471151A747E4E");

    ec->ec_data = (void*)d;

    mpz_clear(x);
}

void ec_bn254_fp2_init_ec_data_beuchat(EC_GROUP ec)
{
    ec_data_fp2 d = NULL;
    struct ec_field_st *fp;

    mpz_t x;

    d = (ec_data_fp2)malloc(sizeof(struct ec_bn254_fp2_ec_data_st));

    mpz_init(x);
    mpz_set_str(x, "3FC0100000000000", 16);

    mpz_init(d->_6x);
    mpz_init(d->_6x2);

    mpz_mul_ui(d->_6x, x, 6);
    mpz_mul(d->_6x2, d->_6x, x);

    fp = ec->field->base;

    element_init(d->vfrobx, fp);
    element_init(d->vfroby, fp);
    element_init(d->vfrobx2, fp);
    element_init(d->vfroby2, fp);
    element_init(d->vfrobx3, fp);
    element_init(d->vfroby3, fp);

    element_set_str(d->vfrobx, "2370FB049D410FBE074D0D4C437281205F408FD005FFFFFF40BFD00000000000");
    element_set_str(d->vfroby, "A6A168991BF611A1A4DBA88BB5C4B1A96E9F8336A4488B464D5AA9B1E46A58E");
    element_set_str(d->vfrobx2, "2370FB049D410FBE074D0D4C437281205F408FD005FFFFFF40BFCFFFFFFFFFFF");
    element_set_str(d->vfroby2, "2370FB049D410FBE4E761A9886E502417D023F40180000017E80600000000000");
    element_set_str(d->vfrobx3, "2370FB049D410FBE4E761A9886E502417D023F40180000017E80600000000000");
    element_set_str(d->vfroby3, "1906E47B0B81AEA43428600FCB88B726E618470CADBB774D19AAB564E1B95A73");

    ec->ec_data = (void*)d;

    mpz_clear(x);
}

void ec_bn254_fp2_clear_ec_data(EC_GROUP ec)
{
    ec_data_fp2 d = (ec_data_fp2)(ec->ec_data);

    mpz_clear(d->_6x);
    mpz_clear(d->_6x2);

    element_clear(d->vfrobx);
    element_clear(d->vfroby);
    element_clear(d->vfrobx2);
    element_clear(d->vfroby2);
    element_clear(d->vfrobx3);
    element_clear(d->vfroby3);

    free(d);
    ec->ec_data = NULL;
}

//------------------------------------------------------
//  Decompose scalar
//    s = s0 + s1[6x] + s2[6x^2] + s3[36x^3]
// 	  x := 2^62 - 2^54 + 2^44  <- Beuchat's parameter
//    x := -(2^62+2^55+1)      <- Aranha's parameter
//------------------------------------------------------
void ec_bn254_fp2_decompose_scalar(mpz_t s0, mpz_t s1, mpz_t s2, mpz_t s3, const mpz_t s, ec_data_fp2 d)
{
    mpz_t A, B;

    mpz_init(A);
    mpz_init(B);

    mpz_fdiv_qr(B, A, s, d->_6x2);   // s = A + B*[6x^2]
    mpz_fdiv_qr(s1, s0, A, d->_6x);  // A = s0 + s1*[6x]
    mpz_fdiv_qr(s3, s2, B, d->_6x);  // B = s2 + s3*[6x]
    mpz_clear(A);
    mpz_clear(B);
}

//-------------------------------------------------------------
//  Frobenius map for point on Twist
//     frob(P) = frob(x,y)
//             = (v^{(p-1)/3}x^p, v^{(p-1)/2}*y^p)
//-------------------------------------------------------------
void ec_bn254_tw_frob(EC_POINT Q, const EC_POINT P)
{
    ec_data_fp2 d;

    if (point_is_infinity(P)) {
        point_set_infinity(Q);
        return;
    }

    d = (ec_data_fp2)(curve(P)->ec_data);

    bn254_fp2_frob_p(xcoord(Q), xcoord(P));
    bn254_fp2_frob_p(ycoord(Q), ycoord(P));

    if (strcmp(P->ec->curve_name, "ec_bn254_twa") == 0)
    {
        bn254_fp2_mul_p(xcoord(Q), xcoord(Q), d->vfrobx);
        bn254_fp2_mul_p(ycoord(Q), ycoord(Q), d->vfroby);
    }

    if (strcmp(P->ec->curve_name, "ec_bn254_twb") == 0)
    {
        bn254_fp2_mul(xcoord(Q), xcoord(Q), d->vfrobx);
        bn254_fp2_mul(ycoord(Q), ycoord(Q), d->vfroby);
    }

    element_set(zcoord(Q), zcoord(P));

    Q->isinfinity = P->isinfinity;
}

//-------------------------------------------------------------
//  Frobenius map for point on Twist
//     frob2(P) = frob2(x,y)
//              = (v^{(p^2-1)/3}x, v^{(p^2-1)/2}*y)
//-------------------------------------------------------------
void ec_bn254_tw_frob2(EC_POINT Q, const EC_POINT P)
{
    ec_data_fp2 d;

    if (point_is_infinity(P)) {
        point_set_infinity(Q);
        return;
    }

    d = (ec_data_fp2)(curve(P)->ec_data);

    if (strcmp(P->ec->curve_name, "ec_bn254_twa") == 0)
    {
        bn254_fp2_mul_p(xcoord(Q), xcoord(P), d->vfrobx2);
        bn254_fp2_neg(ycoord(Q), ycoord(P));
    }

    if (strcmp(P->ec->curve_name, "ec_bn254_twb") == 0)
    {
        bn254_fp2_mul(xcoord(Q), xcoord(P), d->vfrobx2);
        bn254_fp2_mul(ycoord(Q), ycoord(P), d->vfroby2);
    }

    element_set(zcoord(Q), zcoord(P));

    Q->isinfinity = P->isinfinity;
}

//-------------------------------------------------------------
//  Frobenius map for point on Twist
//     frob3(P) = frob3(x,y)
//              = (v^{(p^3-1)/3}x^p, v^{(p^3-1)/2}*y^p)
//-------------------------------------------------------------
void ec_bn254_tw_frob3(EC_POINT Q, const EC_POINT P)
{
    ec_data_fp2 d;

    if (point_is_infinity(P)) {
        point_set_infinity(Q);
        return;
    }

    d = (ec_data_fp2)(curve(P)->ec_data);

    bn254_fp2_frob_p(xcoord(Q), xcoord(P));
    bn254_fp2_frob_p(ycoord(Q), ycoord(P));

    if (strcmp(P->ec->curve_name, "ec_bn254_twa") == 0)
    {
        bn254_fp2_mul_p(xcoord(Q), xcoord(Q), d->vfrobx3);
        bn254_fp2_mul_p(ycoord(Q), ycoord(Q), d->vfroby3);
    }

    if (strcmp(P->ec->curve_name, "ec_bn254_twb") == 0)
    {
        bn254_fp2_mul(xcoord(Q), xcoord(Q), d->vfrobx3);
        bn254_fp2_mul(ycoord(Q), ycoord(Q), d->vfroby3);
    }

    element_set(zcoord(Q), zcoord(P));

    Q->isinfinity = P->isinfinity;
}

//------------------------------------------------------
//  concatinate function for bit
//------------------------------------------------------
int catbit(int b0, int b1, int b2, int b3)
{
    return b0 | (b1 << 1) | (b2 << 2) | (b3 << 3);
}

//---------------------------------------------------------
//  scalar multiplication with endomorphism
//---------------------------------------------------------
void ec_bn254_fp2_mul_end(EC_POINT Q, const mpz_t s, const EC_POINT P)
{
    long t, i, index;

    int t0, t1, t2, t3;

    mpz_t s0, s1, s2, s3;

    ec_data_fp2 d = (ec_data_fp2)(curve(P)->ec_data);

    EC_POINT R[16];

    mpz_init(s0);
    mpz_init(s1);
    mpz_init(s2);
    mpz_init(s3);

    if (ec_bn254_fp2_is_on_curve(P) == 0) {
        ec_bn254_fp2_mul_naf(Q, s, P);
        return;
    }

    //--------------------------------------------
    // s = s0 + s1[6x] + s2[6x^2] + s3[36x^3]
    //--------------------------------------------
    ec_bn254_fp2_decompose_scalar(s0, s1, s2, s3, s, d);

    //--------------------------------------------
    // create table R:
    //--------------------------------------------
    for (i = 0; i < 16; i++) {
        point_init(R[i], curve(P));
    }

    ec_bn254_tw_frob(R[4], P);
    point_add(R[2], P, R[4]);
    point_sub(R[8], P, R[4]);
    ec_bn254_tw_frob3(R[8], R[8]);

    if (strcmp(P->ec->curve_name, "ec_bn254_twa") == 0)
    {
        point_add(R[2], R[2], R[8]);
        point_neg(R[2], R[2]);
    }

    if (strcmp(P->ec->curve_name, "ec_bn254_twb") == 0)
    {
        point_add(R[2], R[2], R[8]);
    }

    ec_bn254_tw_frob(R[8], R[2]);

    point_set(R[1], P);
    point_add(R[3], R[1], R[2]);
    point_add(R[5], R[4], R[1]);
    point_add(R[6], R[4], R[2]);
    point_add(R[7], R[6], R[1]);
    point_add(R[9], R[8], R[1]);
    point_add(R[10], R[8], R[2]);
    point_add(R[11], R[8], R[3]);
    point_add(R[12], R[8], R[4]);
    point_add(R[13], R[8], R[5]);
    point_add(R[14], R[8], R[6]);
    point_add(R[15], R[8], R[7]);

    //--------------------------------------------
    //  init
    //--------------------------------------------
    point_set_infinity(Q);

    t0 = mpz_sizeinbase(s0, 2);
    t1 = mpz_sizeinbase(s1, 2);
    t2 = mpz_sizeinbase(s2, 2);
    t3 = mpz_sizeinbase(s3, 2);

    t = MAX(MAX(MAX(t0, t1), t2), t3);

    //--------------------------------------------
    //  [s]P = [s0, s1, s2, s3] P
    //--------------------------------------------
    for (i = t - 1; i >= 0; i--)
    {
        ec_bn254_fp2_dob_formul(Q, Q);

        index = catbit(mpz_tstbit(s0, i), mpz_tstbit(s1, i), mpz_tstbit(s2, i), mpz_tstbit(s3, i));

        if (index) {
            ec_bn254_fp2_add_formul(Q, Q, R[index]);
        }
    }

    point_make_affine(Q, Q);

    //--------------------------------------------
    //  release
    //--------------------------------------------
    for (i = 0; i < 16; i++) {
        point_clear(R[i]);
    }

    mpz_clear(s0);
    mpz_clear(s1);
    mpz_clear(s2);
    mpz_clear(s3);
}

//---------------------------------------------------------
//   Frobenius Map
//---------------------------------------------------------
void ec_bn254_fp2_frob_p(EC_POINT P, const EC_POINT Q)
{
    if (point_is_infinity(Q)) {
        point_set_infinity(P);
        return;
    }

    bn254_fp2_frob_p(xcoord(P), xcoord(Q));
    bn254_fp2_frob_p(ycoord(P), ycoord(Q));

    element_set_one(zcoord(P));

    P->isinfinity = FALSE;
}

//-------------------------------------------
//  comparison operation
//-------------------------------------------
int ec_bn254_fp2_cmp(const EC_POINT P, const EC_POINT Q)
{
    if (bn254_fp2_cmp(xcoord(P), xcoord(Q)) == 0)
    {
        if (bn254_fp2_cmp(ycoord(P), ycoord(Q)) == 0) {
            return 0;
        }
    }
    return 1;
}

int ec_bn254_fp2_is_on_curve(const EC_POINT P)
{
    int hr = FALSE;

    Element x, y;
    EC_POINT R;

    if (point_is_infinity(P)) { return TRUE; }

    element_init(x, field(P));
    element_init(y, field(P));
    point_init(R, curve(P));

    element_sqr(x, xcoord(P));
    element_mul(x, x, xcoord(P));
    element_add(x, x, curve(P)->b);
    element_sqr(y, ycoord(P));

    ec_bn254_fp2_mul_naf(R, curve(P)->order, P);

    hr = (element_cmp(x, y) == 0  && point_is_infinity(R));

    element_clear(x);
    element_clear(y);
    point_clear(R);

    return hr;
}

//-------------------------------------------
//  make affine, jacobian
//-------------------------------------------
void ec_bn254_fp2_make_affine(EC_POINT z, const EC_POINT x)
{
    if (point_is_infinity(x)) {
        point_set_infinity(z);
    }
    else
    {
        Element iz, iz2;

        element_init(iz, field(x));
        element_init(iz2, field(x));

        element_inv(iz, zcoord(x));
        element_sqr(iz2, iz);
        element_mul(xcoord(z), xcoord(x), iz2);
        element_mul(iz2, iz2, iz);
        element_mul(ycoord(z), ycoord(x), iz2);
        element_set_one(zcoord(z));

        z->isinfinity = FALSE;

        element_clear(iz);
        element_clear(iz2);
    }
}

//-------------------------------------------
//  random point
//-------------------------------------------
void ec_bn254_fp2_random(EC_POINT z)
{
    Element t0, t1, t2;

    const struct ec_field_st *f = field(z);

    element_init(t0, f);
    element_init(t1, f);
    element_init(t2, f);

    do {
        do {
            element_random(t0);      //t0 = random value in GF(p)

            element_sqr(t1, t0);     //t1 = t0^3 + b
            element_mul(t1, t1, t0); //
            element_add(t1, t1, curve(z)->b);

        } while (!element_sqrt(t2, t1));

        ec_bn254_fp2_point_set_xy(z, t0, t2);
        ec_bn254_fp2_mul_naf(z, curve(z)->cofactor, z);  // cofactor mult
    } while (point_is_infinity(z));

    element_clear(t0);
    element_clear(t1);
    element_clear(t2);
}

//===========================================
// map to point
//===========================================

//-------------------------------------------
// BS2FQE for fp2
//------------------------------------------
void bn254_fp2_BS2FQE(Element z, const unsigned char *os, const size_t oslen, int t)
{
    size_t tlen = oslen + 2;
    unsigned char *tmp = (unsigned char*)malloc(sizeof(unsigned char) * (tlen));

    const mpz_t *chr = field_get_char(z->field);

    memset(tmp, 0x00, 2);  // tmp = 0x0000 || os
    memcpy(&(tmp[2]), os, oslen);

    IHF1_SHA(elt_rep(z, 0), tmp, tlen, *chr, t);

    tmp[1] = 1; // tmp = 0x0001 || os

    IHF1_SHA(elt_rep(z, 1), tmp, tlen, *chr, t);

    free(tmp);
}

//------------------------------------------------------
//  Compare two value in integer level
//  input : x, y
//  output :
//     if x > y then return 1
//     if x = y then return 0
//     if x < y then return -1
//------------------------------------------------------
int bn254_fp2_compare(const Element x, const Element y)
{
    int ret = 0;

    ret = mpz_cmp(elt_rep(x, 1), elt_rep(y, 1));

    if (ret > 0) {
        return  1;
    }
    if (ret < 0) {
        return -1;
    }

    ret = mpz_cmp(elt_rep(x, 0), elt_rep(y, 0));

    if (ret > 0) {
        return  1;
    }
    if (ret < 0) {
        return -1;
    }

    return 0;
}

void ec_bn254_fp2_map_to_point(EC_POINT z, const char *s, size_t slen, int t)
{
    mpz_t i;

    unsigned char *d;   // d : For saving hash value of s (octet string)
    unsigned char *id;  // id : i||d (octet string)
    size_t dlen;        // length of d
    size_t idlen;       // length of id

    const struct ec_field_st *f;

    Element x0, y0, y1, y2, t0;

    d = (unsigned char*)malloc(sizeof(unsigned char) * (t / 4));
    id = (unsigned char*)malloc(sizeof(unsigned char) * (t / 4 + 2));

    mpz_init_set_ui(i, 0); // i = 0

    f = field(z);

    element_init(x0, f);
    element_init(y0, f);
    element_init(y1, f);
    element_init(y2, f);
    element_init(t0, f);

    mIHF_SHA(d, &dlen, s, slen, t);   //create digest for input ID

    do
    {
        cat_int_str(id, &idlen, i, d, dlen); // i||d (octet string)

        bn254_fp2_BS2FQE(x0, id, idlen, t); //create x0 by BS2FQE

        bn254_fp2_sqr(t0, x0);  //tmp = x0^3 + b
        bn254_fp2_mul(t0, t0, x0);
        bn254_fp2_add(t0, t0, curve(z)->b);

        if (element_is_zero(t0))
        {
            point_set_xy(z, x0, t0);   //z = (x0, 0)
            goto release;
        }

        mpz_add_ui(i, i, 1);   //i = i+1

    } while (!bn254_fp2_sqrt(y0, t0));

    bn254_fp2_set(y1, y0);   // y1 = y0
    bn254_fp2_neg(y2, y0);   // y2 = -y0

    (bn254_fp2_compare(y1, y2) < 0) ? bn254_fp2_set(y0, y1) : bn254_fp2_set(y0, y2);

    point_set_xy(z, x0, y0);

    ec_bn254_fp2_mul_naf(z, curve(z)->cofactor, z); // cofactor mult

release:
    element_clear(x0);
    element_clear(y0);
    element_clear(y1);
    element_clear(y2);
    element_clear(t0);

    mpz_clear(i);

    free(d);
    free(id);
}

//-------------------------------------------
//  i/o operation (octet string)
//-------------------------------------------
void ec_bn254_fp2_to_oct(unsigned char *os, size_t *size, const EC_POINT P)
{
    size_t sx, sy;

    unsigned char ox[64];
    unsigned char oy[64];

    if (point_is_infinity(P)) {
        os[0] = 0x00;
        (*size) = 1;
        return;
    }

    bn254_fp2_to_oct(ox, &sx, xcoord(P));
    bn254_fp2_to_oct(oy, &sy, ycoord(P));

    os[0] = 0x04;

    memset(&(os[1]), 0x00, 128);

    memcpy(&(os[1]), ox, sx);
    memcpy(&(os[65]), oy, sy);

    (*size) = 129;

    return;
}

void ec_bn254_fp2_from_oct(EC_POINT P, const unsigned char *os, size_t size)
{
    if (size != 1 && size != 129)
    {
        point_set_infinity(P);
        return;
    }

    switch (os[0])
    {
    case 0x00:
        point_set_infinity(P);
        break;
    case 0x04:
        bn254_fp2_from_oct(xcoord(P), &(os[1]), 64);
        bn254_fp2_from_oct(ycoord(P), &(os[65]), 64);
        bn254_fp2_set_one(zcoord(P));
        P->isinfinity = FALSE;
        break;
    }
}
