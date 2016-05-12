//============================================================================
//  elliptic curve over prime field implemented with GMP
//----------------------------------------------------------------------------
//  ec_bn254_fpa : EC_fp := Y^2 = X^3+5 Beuchat et al.
//  ec_bn254_fpb : EC_fp := Y^2 = X^3+2 Aranha et al.
//----------------------------------------------------------------------------
//  2015.10.31 created by kanbara
//============================================================================

#include "ec_bn254_lcl.h"

#define xcoord(p)  (p->x)
#define ycoord(p)  (p->y)
#define zcoord(p)  (p->z)

#define field(p)   (p->ec->field)
#define curve(p)   (p->ec)

#define mpz_rep(x) (*((mpz_t*)((x)->data)))

//------------------------------------------
//  initialization, clear, set
//-------------------------------------------
void ec_bn254_fp_point_init(EC_POINT P)
{
    element_init(xcoord(P), field(P));
    element_init(ycoord(P), field(P));
    element_init(zcoord(P), field(P));

    point_set_infinity(P);
}

void ec_bn254_fp_point_clear(EC_POINT P)
{
    element_clear(xcoord(P));
    element_clear(ycoord(P));
    element_clear(zcoord(P));

    P->ec = NULL;
}

void ec_bn254_fp_point_set(EC_POINT z, const EC_POINT x)
{
    element_set(xcoord(z), xcoord(x));
    element_set(ycoord(z), ycoord(x));
    element_set(zcoord(z), zcoord(x));

    z->isinfinity = x->isinfinity;
}

void ec_bn254_fp_point_set_str(EC_POINT P, const char* s)
{
    static const char infinity[] = "[0]";

    int i, hr;
    int len = strlen(s);

    char *tmp, *p1, *p2;

    if (strcmp(s, infinity) == 0) {
        point_set_infinity(P);
        return;
    }

    tmp = (char *)malloc(sizeof(char) * (len + 1));

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
        fprintf(stderr, "Please check input string : format = [x, y] or [0]\n");
        exit(200);
    }

    bn254_fp_set_str(xcoord(P), p1);
    bn254_fp_set_str(ycoord(P), p2);
    bn254_fp_set_one(zcoord(P));

    P->isinfinity = FALSE;

    free(tmp);
}

void ec_bn254_fp_point_get_str(char *s, const EC_POINT P)
{
    static const char infinity[] = "[0]";

    char sx[65], sy[65];

    EC_POINT T;

    if (point_is_infinity(P)) {
        strcpy(s, infinity);
        return;
    }

    point_init(T, curve(P));

    point_make_affine(T, P);

    bn254_fp_get_str(sx, xcoord(T));
    bn254_fp_get_str(sy, ycoord(T));

    sprintf(s, "[%s,%s]", sx, sy);

    point_clear(T);
}

void ec_bn254_fp_point_set_xy(EC_POINT z, const Element x, const Element y)
{
    element_set(xcoord(z), x);
    element_set(ycoord(z), y);
    element_set_one(zcoord(z));

    z->isinfinity = FALSE;
}

void ec_bn254_fp_point_set_infinity(EC_POINT P)
{
    element_set_one(xcoord(P));
    element_set_one(ycoord(P));
    element_set_zero(zcoord(P));

    P->isinfinity = TRUE;
}

//-------------------------------------------
//  arithmetic operation
//-------------------------------------------
void ec_bn254_fp_add(EC_POINT R, const EC_POINT P, const EC_POINT Q)
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

    if (element_cmp(xcoord(P), xcoord(Q)) == 0)
    {
        if (element_cmp(ycoord(P), ycoord(Q)) == 0)
        {
            ec_bn254_fp_dob(R, P);
            return;
        }
        point_set_infinity(R);
        return;
    }

    bn254_fp_sub(t[0], xcoord(Q), xcoord(P));
    bn254_fp_sub(t[1], ycoord(Q), ycoord(P));
    bn254_fp_inv(t[0], t[0]);
    bn254_fp_mul(t[2], t[1], t[0]);

    bn254_fp_sqr(t[0], t[2]);
    bn254_fp_sub(t[0], t[0], xcoord(P));
    bn254_fp_sub(t[0], t[0], xcoord(Q));

    bn254_fp_sub(t[1], xcoord(P), t[0]);
    bn254_fp_mul(t[1], t[2], t[1]);
    bn254_fp_sub(ycoord(R), t[1], ycoord(P));
    bn254_fp_set(xcoord(R), t[0]);
    bn254_fp_set_one(zcoord(R));

    R->isinfinity = FALSE;
}

void ec_bn254_fp_dob(EC_POINT R, const EC_POINT P)
{
    Element *t = field(R)->tmp;

    if (point_is_infinity(P)) {
        point_set_infinity(R);
        return;
    }
    if (element_is_zero(ycoord(P))) {
        point_set_infinity(R);
        return;
    }

    bn254_fp_add(t[0], ycoord(P), ycoord(P));
    bn254_fp_inv(t[0], t[0]);
    bn254_fp_sqr(t[1], xcoord(P));
    bn254_fp_add(t[2], t[1], t[1]);
    bn254_fp_add(t[1], t[1], t[2]);
    bn254_fp_mul(t[2], t[1], t[0]);

    bn254_fp_sqr(t[0], t[2]);
    bn254_fp_sub(t[0], t[0], xcoord(P));
    bn254_fp_sub(t[0], t[0], xcoord(P));

    bn254_fp_sub(t[1], xcoord(P), t[0]);
    bn254_fp_mul(t[1], t[2], t[1]);
    bn254_fp_sub(ycoord(R), t[1], ycoord(P));
    bn254_fp_set(xcoord(R), t[0]);
    bn254_fp_set_one(zcoord(R));

    R->isinfinity = FALSE;
}

void ec_bn254_fp_neg(EC_POINT z, const EC_POINT x)
{
    if (point_is_infinity(x)) {
        point_set_infinity(z);
        return;
    }

    bn254_fp_set(xcoord(z), xcoord(x));
    bn254_fp_neg(ycoord(z), ycoord(x));
    bn254_fp_set(zcoord(z), zcoord(x));

    z->isinfinity = x->isinfinity;
}

void ec_bn254_fp_sub(EC_POINT z, const EC_POINT x, const EC_POINT y)
{
    EC_POINT t0;

    point_init(t0, curve(z));

    point_neg(t0, y);
    point_add(z, x, t0);

    point_clear(t0);
}

void ec_bn254_fp_add_formul(EC_POINT R, const EC_POINT P, const EC_POINT Q)
{
    Element *t = field(R)->tmp;

    if (strcmp(P->ec->curve_name, "ec_bn254_fpa") == 0)
    {
        if (point_is_infinity(P)) {
            point_set(R, Q);
            return;
        }
        if (point_is_infinity(Q)) {
            point_set(R, P);
            return;
        }

        bn254_fp_sqr(t[0], zcoord(P));     // A = Pz^2
        bn254_fp_mul(t[1], t[0], zcoord(P)); // B = Pz^3
        bn254_fp_mul(t[2], t[0], xcoord(Q)); // C = Qx*A
        bn254_fp_mul(t[1], t[1], ycoord(Q)); // D = Qy*B
        bn254_fp_sub(t[2], t[2], xcoord(P)); // E = C-Px
        bn254_fp_sub(t[1], t[1], ycoord(P)); // F = D-Py

        if (element_is_zero(t[2]))
        {
            if (element_is_zero(t[1]))
            {
                ec_bn254_fp_dob_formul(R, P);
                return;
            }
            point_set_infinity(R);
            return;
        }

        bn254_fp_mul(zcoord(R), zcoord(P), t[2]); // Rz = Pz*E

        bn254_fp_sqr(t[0], t[2]);            // G = E^2
        bn254_fp_mul(t[3], xcoord(P), t[0]); // I = Px*G
        bn254_fp_mul(t[4], t[0], t[2]);      // H = G*E

        bn254_fp_sqr(xcoord(R), t[1]);
        bn254_fp_sub(xcoord(R), xcoord(R), t[4]);
        bn254_fp_add(t[0], t[3], t[3]);
        bn254_fp_sub(xcoord(R), xcoord(R), t[0]); // Rx = F^2-H-2*I

        bn254_fp_sub(t[3], t[3], xcoord(R));
        bn254_fp_mul(t[0], t[4], ycoord(P));
        bn254_fp_mul(ycoord(R), t[1], t[3]);
        bn254_fp_sub(ycoord(R), ycoord(R), t[0]); // Ry = F*(I-Rx)-Py*H
    }

    // Jacoban coordinate proposed by Aranha et al.
    if (strcmp(P->ec->curve_name, "ec_bn254_fpb") == 0)
    {
        if (point_is_infinity(P)) {
            point_set(R, Q);
            return;
        }
        if (point_is_infinity(Q)) {
            point_set(R, P);
            return;
        }

        if (point_cmp(P, Q) == 0) {
            ec_bn254_fp_dob_formul(R, P);
            return;
        }

        bn254_fp_sqr(t[1], zcoord(P));        		// t1 = Pz^2
        bn254_fp_mul(t[3], xcoord(Q), t[1]); 		// t3 = Qx*t1
        bn254_fp_mul(t[1], t[1], zcoord(P));  		// t1 = t1*Pz
        bn254_fp_sub(t[3], t[3], xcoord(P));  		// t3 = t3-Px
        bn254_fp_mul(t[4], t[1], ycoord(Q));  		// t4 = t1*Qy
        bn254_fp_mul(zcoord(R), zcoord(P), t[3]); 	// Rz = Pz*t3

        if (bn254_fp_is_zero(zcoord(R)))
        {
            point_set_infinity(R);
            return;
        }

        bn254_fp_sub(t[0], t[4], ycoord(P));  		// t0 = t4-Py
        bn254_fp_sqr(t[1], t[3]);             		// t1 = t3^2
        bn254_fp_mul(t[4], t[1], t[3]);	   			// t4 = t1*t3
        bn254_fp_mul(t[1], t[1], xcoord(P));  		// t1 = t1*Px
        bn254_fp_sqr(xcoord(R), t[0]);		   		// Rx = t0^2
        bn254_fp_dob(t[3], t[1]);			   		// t3 = 2*t1
        bn254_fp_sub(xcoord(R), xcoord(R), t[3]); 	// Rx = Rx-t3
        bn254_fp_sub(xcoord(R), xcoord(R), t[4]); 	// Rx = Rx-t4
        bn254_fp_sub(t[1], t[1], xcoord(R));  		// t1 = t1-Rx
        bn254_fp_mul(t[2], t[0], t[1]); 			// t2 = t0*t1
        bn254_fp_mul(t[3], t[4], ycoord(P)); 		// t3 = t4*Py
        bn254_fp_sub(t[2], t[2], t[3]); 			// t2 = t2-t3
        bn254_fp_set(ycoord(R), t[2]);  			// Ry = t2 mod p
    }
    R->isinfinity = FALSE;
}

void ec_bn254_fp_dob_formul(EC_POINT R, const EC_POINT P)
{
    Element *t = field(R)->tmp;

    if (point_is_infinity(P)) {
        point_set_infinity(R);
        return;
    }

    if (strcmp(P->ec->curve_name, "ec_bn254_fpa") == 0)
    {
        bn254_fp_sqr(t[0], ycoord(P));     // A = Py^2
        bn254_fp_add(t[1], xcoord(P), xcoord(P));
        bn254_fp_add(t[1], t[1], t[1]);
        bn254_fp_mul(t[1], t[1], t[0]);    // B = 4*Px*A
        bn254_fp_sqr(t[0], t[0]);
        bn254_fp_add(t[0], t[0], t[0]);
        bn254_fp_add(t[0], t[0], t[0]);
        bn254_fp_add(t[0], t[0], t[0]);    // C = 8*A^2
        bn254_fp_sqr(t[3], xcoord(P));
        bn254_fp_add(t[2], t[3], t[3]);
        bn254_fp_add(t[2], t[2], t[3]);    // D = 3*Px^2

        bn254_fp_sqr(xcoord(R), t[2]);
        bn254_fp_add(t[3], t[1], t[1]);
        bn254_fp_sub(xcoord(R), xcoord(R), t[3]); // Rx = D^2 - 2*B

        bn254_fp_mul(zcoord(R), ycoord(P), zcoord(P));
        bn254_fp_add(zcoord(R), zcoord(R), zcoord(P)); // Rz = 2*Py*Pz

        bn254_fp_sub(ycoord(R), t[1], xcoord(R));
        bn254_fp_mul(ycoord(R), ycoord(R), t[2]);
        bn254_fp_sub(ycoord(R), ycoord(R), t[0]);  // Ry = D*(B-Rx)-C
    }

    // Jacoban coordinate proposed by Aranha et al.
    if (strcmp(P->ec->curve_name, "ec_bn254_fpb") == 0)
    {
        bn254_fp_sqr(t[0], xcoord(P)); 				// t0 = Px^2
        bn254_fp_dob(t[1], t[0]);      				// t1 = 2*t0
        bn254_fp_mul(zcoord(R), ycoord(P), zcoord(P)); // Rz = Py*Pz
        bn254_fp_add(t[0], t[0], t[1]); 			// t0 = t0+t1
        bn254_fp_sqr(t[3], ycoord(P));  			// t3 = Py^2
        bn254_fp_div2(t[0], t[0]);			 		// t0 = t0/2
        bn254_fp_mul(t[1], t[3], xcoord(P)); 		// t1 = t3*Px
        bn254_fp_dob(ycoord(R), t[1]); 	  			// Ry = 2*t1
        bn254_fp_sqr(xcoord(R), t[0]);       		// Rx = t0^2
        bn254_fp_sub(xcoord(R), xcoord(R), ycoord(R)); // Rx = Rx-Ry
        bn254_fp_sub(t[1], t[1], xcoord(R)); 		// t1 = t1-Rx
        bn254_fp_sqr(t[2], t[3]); 					// t2 = t3^2
        bn254_fp_mul(t[1], t[0], t[1]);				// t1 = t0*t1
        bn254_fp_sub(t[1], t[1], t[2]); 			// t1 = t1-t2
        bn254_fp_set(ycoord(R), t[1]); 				// Ry = t1 mod p
    }

    R->isinfinity = FALSE;
}

//--------------------------------------------------------------
//  Scalar Multiplication in Affine Coordinate
//--------------------------------------------------------------
void ec_bn254_fp_mul_affine(EC_POINT Q, const mpz_t s, const EC_POINT P)
{
    long t, i;

    EC_POINT R;

    point_init(R, curve(P));
    point_set(R, P);

    t = mpz_sizeinbase(s, 2);

    for (i = t - 2; i >= 0; i--)
    {
        point_dob(R, R);
        if (mpz_tstbit(s, i))
        {
            point_add(R, P, R);
        }
    }

    point_set(Q, R);
    point_clear(R);
}

//--------------------------------------------------------------
//  Scalar Multiplication in Jacobian Coordinate
//--------------------------------------------------------------
void ec_bn254_fp_mul(EC_POINT Q, const mpz_t s, const EC_POINT P)
{
    long t, i;
    EC_POINT R;

    point_init(R, curve(P));

    ec_bn254_fp_point_set(R, P);

    t = mpz_sizeinbase(s, 2);

    for (i = t - 2; i >= 0; i--)
    {
        ec_bn254_fp_dob_formul(R, R);
        if (mpz_tstbit(s, i))
        {
            ec_bn254_fp_add_formul(R, R, P);
        }
    }

    point_make_affine(Q, R);

    point_clear(R);
}

//--------------------------------------------------------------
//  Generate NAF representation of s
//--------------------------------------------------------------
void generate_naf(int *naf, int *len, const mpz_t s)
{
    mpz_t k, r;

    int i = 0;

    mpz_init(k);
    mpz_init(r);

    mpz_abs(k, s);

    while (mpz_cmp_ui(k, 1) >= 0)
    {
        if (mpz_tstbit(k, 0))
        {
            int v = 2 - (int)mpz_mod_ui(r, k, 4);
            if (v > 0) {
                mpz_sub_ui(k, k, v);
            }
            if (v < 0) {
                mpz_add_ui(k, k, -v);
            }
            naf[i] = v;
        }
        else {
            naf[i] = 0;
        }
        mpz_fdiv_q_2exp(k, k, 1);
        i++;
    }

    (*len) = i;

    mpz_clear(k);
    mpz_clear(r);
}

//------------------------------------------------------
//  Scalar Multiplication with NAF
//------------------------------------------------------
void ec_bn254_fp_mul_naf(EC_POINT Q, const mpz_t s, const EC_POINT P)
{
    long t, i;

    int *naf, nlen;

    EC_POINT R, mP;

    point_init(R, curve(P));
    point_init(mP, curve(P));

    point_set(R, P);
    point_neg(mP, P);

    t = mpz_sizeinbase(s, 2);

    naf = (int *)malloc(sizeof(int) * (t + 1));

    generate_naf(naf, &nlen, s);

    for (i = nlen - 2; i >= 0; i--)
    {
        ec_bn254_fp_dob_formul(R, R);
        if (naf[i])
        {
            if (naf[i] < 0) {
                ec_bn254_fp_add_formul(R, R, mP);
            }
            else {
                ec_bn254_fp_add_formul(R, R, P);
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
//
//  beta : {beta}^3 = 1, beta in Fp
//  lambda : [lambda]P = (beta*x, y)
//
//----------------------------------------------------------------
void ec_bn254_fp_init_ec_data(EC_GROUP ec)
{
    ec_data_fp d;

    mpz_t lambda;

    d = (ec_data_fp)malloc(sizeof(struct ec_bn254_fp_ec_data_st));

    if (strcmp(ec->curve_name, "ec_bn254_fpa") == 0)
    {
        element_init(d->beta, ec->field);
        element_set_str(d->beta, "2370FB049D410FBE074D0D4C437281205F408FD005FFFFFF40BFCFFFFFFFFFFF");

        mpz_init(lambda);
        mpz_set_str(lambda, "2370FB049D410FBDC023FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF", 16);
    }

    if (strcmp(ec->curve_name, "ec_bn254_fpb") == 0)
    {
        element_init(d->beta, ec->field);
        element_set_str(d->beta, "25236482400000017080EB4000000006181800000000000CD98000000000000B");

        mpz_init(lambda);
        mpz_set_str(lambda, "252364824000000126CD8900000000024908FFFFFFFFFFFCF9FFFFFFFFFFFFF6", 16);
    }

    mpz_init(d->n);
    mpz_init(d->n2);
    mpz_init(d->a1);
    mpz_init(d->a2);
    mpz_init(d->b1);
    mpz_init(d->b2);

    mpz_set(d->n, *curve_get_order(ec));
    mpz_add(d->n2, d->n, d->n);

    ec_bn254_fp_decompose_scalar_init(d->a1, d->a2, d->b1, d->b2, d->n, lambda);

    ec->ec_data = (void*)d;

    mpz_clear(lambda);
}

void ec_bn254_fp_clear_ec_data(EC_GROUP ec)
{
    ec_data_fp d = (ec_data_fp)(ec->ec_data);

    element_clear(d->beta);

    mpz_clear(d->n);
    mpz_clear(d->n2);
    mpz_clear(d->a1);
    mpz_clear(d->a2);
    mpz_clear(d->b1);
    mpz_clear(d->b2);

    free(d);

    ec->ec_data = NULL;
}

//---------------------------------------------------------------
//  Endomorphism in Elliptic Curve E
//     E/Fp : Y^2 = X^3 + b ( p mod 3 = 1 )
//
//     P = (x,y) -> [lambda]P = (beta*x, y)
//
//       lambda^2 + lambda + 1 = 0 mod r
//---------------------------------------------------------------
void ec_bn254_fp_point_endomorphism(EC_POINT Q, const EC_POINT P)
{
    ec_data_fp d;

    if (point_is_infinity(P)) {
        point_set_infinity(Q);
        return;
    }

    d = (ec_data_fp)(curve(P)->ec_data);

    element_mul(xcoord(Q), xcoord(P), d->beta);
    element_set(ycoord(Q), ycoord(P));
    element_set(zcoord(Q), zcoord(P));

    Q->isinfinity = P->isinfinity;
}

//----------------------------------------------------------------------------------
//  Init function of Decompose scalar
//     k = ( k1 + k2*l ) mod n
//----------------------------------------------------------------------------------
void ec_bn254_fp_decompose_scalar_init(mpz_t a1, mpz_t a2, mpz_t b1, mpz_t b2, const mpz_t n, const mpz_t l)
{
    mpz_t q, _n;
    mpz_t r0, r1, r2;
    mpz_t s0, s1, s2;
    mpz_t t0, t1, t2;

    mpz_init(q);
    mpz_init(_n);
    mpz_init(s0);
    mpz_init(t0);
    mpz_init(r0);
    mpz_init(s1);
    mpz_init(t1);
    mpz_init(r1);
    mpz_init(r2);
    mpz_init(s2);
    mpz_init(t2);

    mpz_set_ui(s0, 1);
    mpz_set_ui(t0, 0);
    mpz_set(r0, n);
    mpz_set_ui(s1, 0);
    mpz_set_ui(t1, 1);
    mpz_set(r1, l);

    mpz_sqrt(_n, n);

    while (mpz_sgn(r1) != 0)
    {
        mpz_fdiv_qr(q, r2, r0, r1);
        mpz_mul(s2, q, s1);
        mpz_sub(s2, s0, s2);
        mpz_mul(t2, q, t1);
        mpz_sub(t2, t0, t2);
        if (mpz_cmp(r1, _n) < 0) {
            break;
        }
        mpz_set(r0, r1);
        mpz_set(r1, r2);
        mpz_set(s0, s1);
        mpz_set(s1, s2);
        mpz_set(t0, t1);
        mpz_set(t1, t2);
    }

    mpz_set(a1, r1);
    mpz_neg(b1, t1);

    mpz_mul(s0, r0, r0);
    mpz_mul(s1, t0, t0);
    mpz_add(s0, s0, s1);
    mpz_mul(s1, r2, r2);
    mpz_mul(s2, t2, t2);
    mpz_add(s1, s1, s2);

    if (mpz_cmp(s0, s1) <= 0) {
        mpz_set(a2, r0);
        mpz_neg(b2, t0);
    }
    else {
        mpz_set(a2, r2);
        mpz_neg(b2, t2);
    }

    if (mpz_sgn(b2) < 0) {
        mpz_swap(a1, a2);
        mpz_swap(b1, b2);
    }

    mpz_clear(q);
    mpz_clear(_n);
    mpz_clear(s0);
    mpz_clear(t0);
    mpz_clear(r0);
    mpz_clear(s1);
    mpz_clear(t1);
    mpz_clear(r1);
    mpz_clear(r2);
    mpz_clear(s2);
    mpz_clear(t2);
}

//----------------------------------------------------------------------------------
//  Decompose scalar
//     k = ( k1 + k2*l ) mod n
//----------------------------------------------------------------------------------
void ec_bn254_fp_decompose_scalar(mpz_t k1, mpz_t k2, const mpz_t k, ec_data_fp d)
{
    mpz_t c1, c2, t1, t2;

    mpz_init(c1);
    mpz_init(c2);
    mpz_init(t1);
    mpz_init(t2);

    mpz_mul(t1, d->b2, k);
    mpz_add(t1, t1, t1);
    mpz_mul(t2, d->b1, k);
    mpz_add(t2, t2, t2);
    mpz_neg(t2, t2);

    mpz_abs(c1, t1);
    mpz_add(c1, c1, d->n);
    mpz_fdiv_q(c1, c1, d->n2);
    mpz_mul_si(c1, c1, mpz_sgn(t1));
    mpz_abs(c2, t2);
    mpz_add(c2, c2, d->n);
    mpz_fdiv_q(c2, c2, d->n2);
    mpz_mul_si(c2, c2, mpz_sgn(t2));

    mpz_mul(t1, c1, d->a1);
    mpz_mul(t2, c2, d->a2);
    mpz_sub(k1, k, t1);
    mpz_sub(k1, k1, t2);
    mpz_mul(t2, c1, d->b1);
    mpz_neg(k2, t2);
    mpz_mul(t2, c2, d->b2);
    mpz_sub(k2, k2, t2);

    mpz_clear(c1);
    mpz_clear(c2);
    mpz_clear(t1);
    mpz_clear(t2);
}

//-----------------------------------------------------
//  Scalar Multiplication with Endomorphism
//-----------------------------------------------------
void ec_bn254_fp_mul_end(EC_POINT Q, const mpz_t s, const EC_POINT P)
{
    mpz_t s1, s2;

    int sl1, sl2, *sn1, *sn2;
    int i, t;

    EC_POINT P1, P2, mP1, mP2;

    ec_data_fp d;

    mpz_init(s1);
    mpz_init(s2);

    d = (ec_data_fp)(curve(P)->ec_data);

    ec_bn254_fp_decompose_scalar(s1, s2, s, d);

    sl1 = mpz_sizeinbase(s1, 2);
    sl2 = mpz_sizeinbase(s2, 2);

    t = (sl1 > sl2) ? sl1 : sl2;

    sn1 = (int*)malloc(sizeof(int) * (t + 1));
    sn2 = (int*)malloc(sizeof(int) * (t + 1));

    generate_naf(sn1, &sl1, s1);
    generate_naf(sn2, &sl2, s2);

    t = (sl1 > sl2) ? sl1 : sl2;

    for (i = sl1; i < t; i++) {
        sn1[i] = 0;
    }
    for (i = sl2; i < t; i++) {
        sn2[i] = 0;
    }

    if (mpz_sgn(s1) < 0) {
        for (i = 0; i < sl1; i++) {
            sn1[i] *= -1;
        }
    }
    if (mpz_sgn(s2) < 0) {
        for (i = 0; i < sl2; i++) {
            sn2[i] *= -1;
        }
    }

    point_init(P1, curve(P));
    point_init(P2, curve(P));
    point_init(mP1, curve(P));
    point_init(mP2, curve(P));

    ec_bn254_fp_point_set(P1, P);
    ec_bn254_fp_point_endomorphism(P2, P1);

    ec_bn254_fp_neg(mP1, P1);
    ec_bn254_fp_neg(mP2, P2);

    ec_bn254_fp_point_set_infinity(Q);

    for (i = t - 1; i >= 0; i--)
    {
        ec_bn254_fp_dob_formul(Q, Q);

        if (sn1[i])
        {
            if (sn1[i] < 0) {
                ec_bn254_fp_add_formul(Q, Q, mP1);
            }
            else {
                ec_bn254_fp_add_formul(Q, Q, P1);
            }
        }

        if (sn2[i])
        {
            if (sn2[i] < 0) {
                ec_bn254_fp_add_formul(Q, Q, mP2);
            }
            else {
                ec_bn254_fp_add_formul(Q, Q, P2);
            }
        }
    }

    point_make_affine(Q, Q);

    mpz_clear(s1);
    mpz_clear(s2);

    point_clear(P1);
    point_clear(P2);
    point_clear(mP1);
    point_clear(mP2);

    free(sn1);
    free(sn2);
}

//-------------------------------------------
//  comparison operation
//-------------------------------------------
int ec_bn254_fp_is_infinity(const EC_POINT P)
{
    return (P->isinfinity == TRUE);
}

int ec_bn254_fp_is_on_curve(const EC_POINT P)
{
    int hr = FALSE;

    Element x, y;

    if (point_is_infinity(P)) {
        return TRUE;
    }

    element_init(x, field(P));
    element_init(y, field(P));

    element_sqr(x, xcoord(P));
    element_mul(x, x, xcoord(P));
    element_add(x, x, curve(P)->b);
    element_sqr(y, ycoord(P));

    hr = (element_cmp(x, y) == 0);

    element_clear(x);
    element_clear(y);

    return hr;
}

int ec_bn254_fp_cmp(const EC_POINT P, const EC_POINT Q)
{
    if (element_cmp(xcoord(P), xcoord(Q)) == 0)
    {
        if (element_cmp(ycoord(P), ycoord(Q)) == 0) {
            return 0;
        }
    }
    return 1;
}

//-------------------------------------------
//  make affine, jacobian
//-------------------------------------------
void ec_bn254_fp_make_affine(EC_POINT z, const EC_POINT x)
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
//  random and map to point
//-------------------------------------------
void ec_bn254_fp_random(EC_POINT z)
{
    Element t0, t1, t2;

    const struct ec_field_st *f = field(z);

    element_init(t0, f);
    element_init(t1, f);
    element_init(t2, f);

    do {
        element_random(t0);      //t0 = random value in GF(p)

        element_sqr(t1, t0);     //t1 = t0^3 + b
        element_mul(t1, t1, t0); //
        element_add(t1, t1, curve(z)->b);

    } while (!element_sqrt(t2, t1));

    point_set_xy(z, t0, t2);

    element_clear(t0);
    element_clear(t1);
    element_clear(t2);
}

//===========================================
//  map to point
//===========================================

//-------------------------------------------
// BS2FQE for fp
//------------------------------------------
void bn254_fp_BS2FQE(Element z, const unsigned char *os, const size_t oslen, int t)
{
    size_t tlen = oslen + 2;
    unsigned char *tmp = (unsigned char *)malloc(sizeof(unsigned char) * (tlen));

    memset(tmp, 0x00, 2);          // os0 = 0 || os
    memcpy(&(tmp[2]), os, oslen);  //

    IHF1_SHA(mpz_rep(z), tmp, tlen, *field_get_char(z->field), t);

    free(tmp);
}

//---------------------------------------------
// this function calculate i||str (octet string)
//----------------------------------------------
void cat_int_str(unsigned char *os, size_t *oslen, const mpz_t i, const unsigned char *s, const size_t slen)
{
    unsigned char ios[2];
    size_t ilen;

    memset(ios, 0x00, 2);
    mpz_export(ios, &ilen, -1, sizeof(*ios), 1, 0, i);     //change i mpz_t to octet string

    (*oslen) = slen + 2;      //oslen = slen + ilen

    os[0] = ios[1];
    os[1] = ios[0];

    memcpy(&(os[2]), s, slen);
}

void ec_bn254_fp_map_to_point(EC_POINT z, const char *s, size_t slen, int t)
{
    mpz_t i;              // counter i

    unsigned char *d;     // d : For saving hash value of s (octet string)
    unsigned char *id;    // id : i||d (octet string)
    size_t dlen;          // length of d
    size_t idlen;         // length of id

    const struct ec_field_st *f;

    Element x0, y0, y1, y2, t0;

    d = (unsigned char *)malloc(sizeof(unsigned char) * (t / 4));
    id = (unsigned char *)malloc(sizeof(unsigned char) * (t / 4 + 2));

    mpz_init_set_ui(i, 0);      // i = 0

    f = field(z);

    element_init(x0, f);
    element_init(y0, f);
    element_init(y1, f);
    element_init(y2, f);
    element_init(t0, f);

    mIHF_SHA(d, &dlen, s, slen, t); //create digest for input ID

    do
    {
        cat_int_str(id, &idlen, i, d, dlen); // i||d (octet string)

        bn254_fp_BS2FQE(x0, id, idlen, t); //create x0 by BS2FQE

        bn254_fp_sqr(t0, x0);    // t0 = x0^3 + b
        bn254_fp_mul(t0, t0, x0);
        bn254_fp_add(t0, t0, curve(z)->b);

        if (bn254_fp_is_zero(t0))
        {
            ec_bn254_fp_point_set_xy(z, x0, t0);   //z = (x0, 0)
            goto release;
        }

        mpz_add_ui(i, i, 1);   //i = i+1

    } while (!bn254_fp_sqrt(y0, t0));

    bn254_fp_set(y1, y0);   //y1 = y0
    bn254_fp_neg(y2, y0);   //y2 = -y0

    (mpz_cmp(mpz_rep(y1), mpz_rep(y2)) > 0) ? bn254_fp_set(y0, y2) : bn254_fp_set(y0, y1);

    point_set_xy(z, x0, y0);

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
void ec_bn254_fp_to_oct(unsigned char *os, size_t *size, const EC_POINT P)
{
    size_t sx, sy;

    unsigned char ox[32];
    unsigned char oy[32];

    if (point_is_infinity(P)) {
        os[0] = 0x00;
        (*size) = 1;
        return;
    }

    bn254_fp_to_oct(ox, &sx, xcoord(P));
    bn254_fp_to_oct(oy, &sy, ycoord(P));

    os[0] = 0x04;

    memset(&(os[1]), 0x00, 64);

    memcpy(&(os[1]), ox, sx);
    memcpy(&(os[33]), oy, sy);

    (*size) = 65;

    return;
}

void ec_bn254_fp_from_oct(EC_POINT z, const unsigned char *os, size_t size)
{
    if (size != 1 && size != 65)
    {
        point_set_infinity(z);
        return;
    }

    switch (os[0])
    {
    case 0x00:
        point_set_infinity(z);
        break;
    case 0x04:
        bn254_fp_from_oct(z->x, &(os[1]), 32);
        bn254_fp_from_oct(z->y, &(os[33]), 32);
        bn254_fp_set_one(z->z);
        z->isinfinity = FALSE;
        break;
    }
}
