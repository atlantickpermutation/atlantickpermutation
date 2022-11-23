#ifndef KECCAK_FIPS202_H
#define KECCAK_FIPS202_H
#define __STDC_WANT_LIB_EXT1__ 1
#include <stdint.h>
#include <stdlib.h>

#define decshake(bits) \
    int shake##bits(uint8_t *, size_t, const uint8_t *, size_t);

#define decsha3(bits) \
    int sha3_##bits(uint8_t *, size_t, const uint8_t *, size_t);

decshake(128)
    decshake(192)
        decshake(256)
            decsha3(224)
                decsha3(256)
                    decshake(320)
                        decsha3(384)
                            decsha3(512)
                                decsha3(768)
                                    decsha3(896)
                                        decshake(1024)
                                            decshake(2048)
                                                decshake(4096)
                                                    decshake(8192)
                                                        decshake(16384)
#endif