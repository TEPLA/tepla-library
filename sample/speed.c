#include <stdio.h>
#include <time.h>

#include <tepla/ec.h>

#define N 1000

int main(void)
{
    EC_PAIRING p;

    clock_t t1, t2, t_tmp;

    EC_POINT P1, P2;
    EC_POINT Q1, Q2;

    Element a, b, g;

    mpz_t exp, order;
    gmp_randstate_t s;

    int i;

    pairing_init(p, "ECBN254b");

    point_init(P1, p->g1);
    point_init(P2, p->g1);
    point_init(Q1, p->g2);
    point_init(Q2, p->g2);

    element_init(a, p->g3);
    element_init(b, p->g3);
    element_init(g, p->g3);

    mpz_init(exp);
    mpz_init(order);
    gmp_randinit_default(s);
    mpz_set(order, *pairing_get_order(p));

    t_tmp = 0;
    for (i = 0; i < N; i++)
    {
        point_random(P1);
        mpz_urandomm(exp, s, order);
        t1 = clock();
        point_mul(P2, exp, P1);
        t2 = clock();
        t_tmp += t2 - t1;
    }
    printf("Scalar Multiplication on G1: %2.6lf [msec]\n", ((double)t_tmp * 1000) / N / CLOCKS_PER_SEC);

    t_tmp = 0;
    for (i = 0; i < N; i++)
    {
        point_random(Q1);
        mpz_urandomm(exp, s, order);
        t1 = clock();
        point_mul(Q2, exp, Q1);
        t2 = clock();
        t_tmp += t2 - t1;
    }
    printf("Scalar Multiplication on G2: %2.6lf [msec]\n", ((double)t_tmp * 1000) / N / CLOCKS_PER_SEC);

    t_tmp = 0;
    for (i = 0; i < N; i++)
    {
        point_random(P1);
        point_random(Q1);
        t1 = clock();
        pairing_map(g, P1, Q1, p);
        t2 = clock();
        t_tmp += t2 - t1;
    }
    printf("Pairing G2*G1->Gt: %2.6lf [msec]\n", ((double)t_tmp * 1000) / N / CLOCKS_PER_SEC);

    t_tmp = 0;
    for (i = 0; i < N; i++)
    {
        element_random(a);
        mpz_urandomm(exp, s, order);
        t1 = clock();
        element_pow(b, a, exp);
        t2 = clock();
        t_tmp += t2 - t1;
    }
    printf("Power on Gt: %2.6lf [msec]\n", ((double)t_tmp * 1000) / N / CLOCKS_PER_SEC);

    t_tmp = 0;
    for (i = 0; i < N; i++)
    {
        char msg[] = "abc";
        t1 = clock();
        point_map_to_point(P1, msg, sizeof(msg), 128);
        t2 = clock();
        t_tmp += t2 - t1;
    }
    printf("Map to Point on G1: %2.6lf [msec]\n", ((double)t_tmp * 1000) / N / CLOCKS_PER_SEC);

    t_tmp = 0;
    for (i = 0; i < N; i++)
    {
        char msg[] = "abc";
        t1 = clock();
        point_map_to_point(Q1, msg, sizeof(msg), 128);
        t2 = clock();
        t_tmp += t2 - t1;
    }
    printf("Map to Point on G2: %2.6lf [msec]\n", ((double)t_tmp * 1000) / N / CLOCKS_PER_SEC);

    element_clear(a);
    element_clear(b);
    element_clear(g);

    point_clear(P1);
    point_clear(P2);
    point_clear(Q1);
    point_clear(Q2);

    mpz_clear(exp);
    mpz_clear(order);
    gmp_randclear(s);

    pairing_clear(p);

    return 0;
}
