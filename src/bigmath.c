#include <ctype.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>

#include "../include/bignum.h"

/*

void bn_jacobi(bignum* r, bignum* a, bignum* m) {

        if (bn_is_zero(m) || bn_is_even(m)) {
                bn_set_u64(r, 0);
                return;
        }

        bignum a_rem, ua, m_1, t, a_1, m_2, z_0;
        bn_init(&a_rem);
        bn_init(&ua);
        bn_init(&m_1);
        bn_init(&a_1);
        bn_init(&m_2);
        bn_init(&z_0);
        bn_mod(&a_rem, a, m);
        if (a_rem.is_neg) {
                bn_add(&a_rem, &a_rem, m);
        }

        bn_init(&t);
        bn_set_u64(&t, 1);
        bn_copy(&m_1, m);
        bn_rshift1(&m_1);
        bn_and_u64(&m_1);

        while (!bn_is_zero(&ua)) {
                u64 z = bn_cnt_led_zero(&ua);
                bn_rshift(&ua, &ua, z);

                bn_copy(&a_1, &ua);
                bn_rshift1(&a_1);
                bn_and_u64(&a_1, 1);

                bn_copy(&m_2, m);
                bn_rshift(&m_2, &m_2, 2);
                bn_and_u64(m_2, 1);

                bn_copy(&z_0, &z);
                bn_add_u64(&z, &z , 1);

                // TODO :
                // if ((z_0 & (m_1 ^ m_2)) ^ (a_1 & m_1)) {
                //	t = -t;
                // }

                bignum tmp;
                bn_init(&tmp);
                bn_copy(&tmp, &ua);
                bn_copy(&ua, m);
                bn_copy(m, &tmp);
                bn_copy(&m_1, &a_1);
                bn_mod(&ua, &ua, m);
        }

        bignum one;
        bn_init(&one);
        bn_set_u64(&one, 1);

        if (bn_cmp(m, &one) != 0) {
                bn_set_u64(r, 0);
                return;
        } else {
                bn_copy(r, &t);
                return;
        }
}


        */
