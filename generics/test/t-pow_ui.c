/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include "generics.h"

void
test_pow_ui(flint_rand_t state, const ring_t ring, const long * size, long iters)
{
    long iter;

    for (iter = 0; iter < iters; iter++)
    {
        elem_ptr A, B, C, D, E, F;
        ulong e, f;

        A = elem_new(ring);
        B = elem_new(ring);
        C = elem_new(ring);
        D = elem_new(ring);
        E = elem_new(ring);
        F = elem_new(ring);

        e = n_randint(state, 8);

        elem_randtest(A, state, size, ring);
        elem_randtest(B, state, size, ring);
        elem_pow_ui(B, A, e, ring);
        elem_pow_ui(A, A, e, ring);
        if (!elem_equal(B, A, ring))
        {
            printf("FAIL: aliasing B, A\n");
            ring_print(ring); printf("\n\n");
            elem_print(A, ring); printf("\n\n");
            elem_print(B, ring); printf("\n\n");
            abort();
        }

        elem_randtest(A, state, size, ring);
        elem_randtest(B, state, size, ring);
        elem_mul(C, A, B, ring);
        elem_pow_ui(D, C, e, ring);
        elem_pow_ui(E, A, e, ring);
        elem_pow_ui(F, B, e, ring);
        elem_mul(F, E, F, ring);
        if (!elem_equal(D, F, ring))
        {
            printf("FAIL: (A * B)^e = A^e * B^e \n");
            ring_print(ring); printf("\n\n");
            elem_print(A, ring); printf("\n\n");
            elem_print(B, ring); printf("\n\n");
            elem_print(C, ring); printf("\n\n");
            elem_print(D, ring); printf("\n\n");
            elem_print(E, ring); printf("\n\n");
            elem_print(F, ring); printf("\n\n");
            abort();
        }

        e = n_randint(state, 8);
        f = n_randint(state, 8);
        elem_randtest(A, state, size, ring);
        elem_pow_ui(B, A, e, ring);
        elem_pow_ui(C, A, f, ring);
        elem_mul(D, B, C, ring);
        elem_pow_ui(E, A, e + f, ring);
        if (!elem_equal(D, E, ring))
        {
            printf("FAIL: A^e * A^f = A^(e+f) \n");
            ring_print(ring); printf("\n\n");
            elem_print(A, ring); printf("\n\n");
            elem_print(B, ring); printf("\n\n");
            elem_print(C, ring); printf("\n\n");
            elem_print(D, ring); printf("\n\n");
            elem_print(E, ring); printf("\n\n");
            abort();
        }

        elem_del(A, ring);
        elem_del(B, ring);
        elem_del(C, ring);
        elem_del(D, ring);
        elem_del(E, ring);
        elem_del(F, ring);
    }
}

int main()
{
    flint_rand_t state;
    long i;

    printf("pow_ui....");
    fflush(stdout);

    flint_randinit(state);

    /* polynomials over (fmpz) integers */
    {
        ring_t Z, Zx, Zxy, Zxyz;
        long size[4] = {6, 6, 6, 6};

        ring_init_fmpz(Z);
        ring_init_poly(Zx, Z);
        ring_init_poly(Zxy, Zx);
        ring_init_poly(Zxyz, Zxy);

        test_pow_ui(state, Z, size, 1000);
        test_pow_ui(state, Zx, size, 1000);
        test_pow_ui(state, Zxy, size, 10);
        /* slow test_pow_ui(state, Zxyz, size, 1000); */

        ring_clear(Zxyz);
        ring_clear(Zxy);
        ring_clear(Zx);
        ring_clear(Z);
    }

    /* polynomials over (fmpz) integers mod n */
    for (i = 0; i < 20; i++)
    {
        ring_t Z, Zn, Znx, Znxy, Znxyz;
        fmpz_t mod;
        long size[4] = {6, 6, 6, 6};

        ring_init_fmpz(Z);

        fmpz_init(mod);
        fmpz_set_ui(mod, n_randtest_prime(state, 0));
        ring_init_mod(Zn, Z, mod);

        ring_init_poly(Znx, Zn);
        ring_init_poly(Znxy, Znx);
        ring_init_poly(Znxyz, Znxy);

        test_pow_ui(state, Zn, size, 10);
        test_pow_ui(state, Znx, size, 10);
        test_pow_ui(state, Znxy, size, 10);
        /* slow test_pow_ui(state, Znxyz, size, 10); */

        ring_clear(Znxyz);
        ring_clear(Znxy);
        ring_clear(Znx);
        ring_clear(Zn);
        fmpz_clear(mod);
        ring_clear(Z);
    }

    /* polynomials over (nmod) integers mod n */
    for (i = 0; i < 20; i++)
    {
        ring_t Z, Zn, Znx, Znxy, Znxyz;
        mp_limb_t mod;
        long size[4] = {6, 6, 6, 6};

        ring_init_limb(Z);
        mod = n_randtest_prime(state, 0);
        ring_init_mod(Zn, Z, &mod);

        ring_init_poly(Znx, Zn);
        ring_init_poly(Znxy, Znx);
        ring_init_poly(Znxyz, Znxy);

        test_pow_ui(state, Zn, size, 10);
        test_pow_ui(state, Znx, size, 10);
        test_pow_ui(state, Znxy, size, 10);
        /* slow test_pow_ui(state, Znxyz, size, 10); */

        ring_clear(Znxyz);
        ring_clear(Znxy);
        ring_clear(Znx);
        ring_clear(Zn);
        ring_clear(Z);
    }

    /* polynomials over (fmpz) integer fractions */
    {
        ring_t Z, Zx, Zxy, Zq, Zqx, Zxq, Zqxy, Zxqy, Zxyq;
        long size[4] = {6, 6, 6, 6};

        ring_init_fmpz(Z);
        ring_init_poly(Zx, Z);
        ring_init_poly(Zxy, Zx);
        ring_init_frac(Zq, Z, Z);
        ring_init_poly(Zqx, Zq);
        ring_init_frac(Zxq, Zx, Z);
        ring_init_poly(Zqxy, Zqx);
        ring_init_poly(Zxqy, Zxy);
        ring_init_frac(Zxyq, Zxy, Z);

        test_pow_ui(state, Zq, size, 100);
        test_pow_ui(state, Zqx, size, 100);
        test_pow_ui(state, Zxq, size, 100);
        test_pow_ui(state, Zqxy, size, 10);
        test_pow_ui(state, Zxqy, size, 10);
        test_pow_ui(state, Zxyq, size, 10);

        ring_clear(Zxyq);
        ring_clear(Zxqy);
        ring_clear(Zqxy);
        ring_clear(Zxq);
        ring_clear(Zqx);
        ring_clear(Zq);
        ring_clear(Zxy);
        ring_clear(Zx);
        ring_clear(Z);
    }

    /* Complex numbers */
    {
        ring_t Z, Zx, Zi, Zxi, Zix, Zq, Zqi;
        long size[4] = {6, 6, 6, 6};

        ring_init_fmpz(Z);
        ring_init_poly(Zx, Z);
        ring_init_complex(Zi, Z);
        ring_init_complex(Zxi, Zx);
        ring_init_poly(Zix, Zi);
        ring_init_frac(Zq, Z, Z);
        ring_init_complex(Zqi, Zq);

        test_pow_ui(state, Zi, size, 100);
        test_pow_ui(state, Zxi, size, 100);
        test_pow_ui(state, Zix, size, 100);
        test_pow_ui(state, Zqi, size, 100);

        ring_clear(Zqi);
        ring_clear(Zq);
        ring_clear(Zix);
        ring_clear(Zxi);
        ring_clear(Zi);
        ring_clear(Zx);
        ring_clear(Z);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

