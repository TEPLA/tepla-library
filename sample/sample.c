#include <stdio.h>
#include <stdlib.h>

#include <tepla/ec.h>

int main(void)
{
    EC_PAIRING p;
    EC_POINT a, b, c;
    Element d;

    pairing_init(p, "ECBN254a");

    point_init(a, p->g1);
    point_init(b, p->g2);
    point_init(c, p->g1);

    element_init(d, p->g3);

    point_set_str(a,
                  "["
                  "0000000000000000000000000000000000000000000000000000000000000001,"
                  "0D45589B158FAAF6AB0E4AD38D998E9982E7FF63964EE1460342A592677CCCB0"
                  "]"
                 );

    point_set_str(b,
                  "["
                  "19850140BC38957238BDEB56EC7B97FE30A6A65D15C4BA07CEF54DB5026C7210 "
                  "1DEB7F4B6C1AEFAEBD0EB750B841BD8ABF916EB750FDF7291F99DFD290C28CE0,"
                  "14C164D6D18CBC7F64559076E00789C75FF001D1BE0968D210C19FB0D3AD649A "
                  "059A2ABA101B7A3C1FA3CAF4DF6B38F2CB4976287488E33F526FA7E8C5441B4B"
                  "]"
                 );

    pairing_map(d, a, b, p);

    point_print(a);
    point_print(b);

    element_print(d);

    char msg[] = "abc";

    point_map_to_point(c, msg, sizeof(msg), 80);

    point_print(c);

    point_clear(a);
    point_clear(b);
    point_clear(c);

    element_clear(d);
    pairing_clear(p);

    return 0;
}
