=============================================================
 TEPLA Elliptic Curve and Pairing Library
=============================================================

TEPLA (TEPLA Elliptic Curve and Pairing Library) is a
free C library providing functions required to develop pairing based
cryptographis systems. The aim of TEPLA is to provide a library which
can be used in various platform.

TEPLA use GNU MP library and OpenSSL.

TEPLA provide following function

- Finite Field Arithmetic with 254-bit prime number
- Elliptic Curve Arithmetic on Barreto-Neahrig Curve.
- Pairing Arithmetic using Optimal Ate Pairing over BN Curve.

=============
 Components
=============

File composition is as follows

src : source directory
  > tepla
    - ec_lib.c : interfaces of operation function.
    - hash.c : interfaces of hash function.
    >> ec_bn254
      - bn254_fp.c  : prime field.
      - bn254_fp2.c : quadratic extension field.
      - bn254_fp6.c : sextic extension field.
      - bn254_fp12.c : twelvetic extension field.
      - ec_bn254_fp.c : elliptic curve over prime field.
      - ec_bn254_fp2.c : twisted elliptic curve over quadratic extension field.
      - ec_bn254_pairing.c : pairing function.
      - ec_bn254_lib.c : functions to create a instance of pairing structure.
      >>> test
        - test programs (e.g. test_bn254_fp.c, test_bn254_fp12.c)

include : include directory
  > tepla
    - ec.h : declaration of curve group, point and pairing structures.
    - hash.h : declaration hash funtions.
    >> ec
      - field.h : declaration of field, element.

sample : sample source
  - sample.c : basic sample code.
  - speed.c : for evaluation of running time.
