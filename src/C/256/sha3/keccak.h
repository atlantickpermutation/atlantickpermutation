#ifndef KECCAK_FIPS202_H
#define KECCAK_FIPS202_H
#define __STDC_WANT_LIB_EXT1__ 1
#include <stdint.h>
#include <stdlib.h>

#define decshake(bits) \
    int shake##bits(uint8_t *, size_t, const uint8_t *, size_t);

#define decsha3(bits) \
    int sha3_##bits(uint8_t *, size_t, const uint8_t *, size_t);

decsha3(256)
#endif