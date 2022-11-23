#include <stdio.h>
#include <stdint.h>

// Initial Constants from MD6 constants: Constante MSB_32(Const 166 MD6)
// #define delta 0x1c2cf22c

// Shift Constants from Atlantick
// ------- Contants a ----------
int a[] = {1, 9, 12, 5, 11, 15, 13, 3, 16, 6, 4, 8, 14, 7, 10, 2};

// ------- Contants b ----------
int b[] = {8, 12, 10, 6, 15, 5, 2, 3, 16, 1, 14, 4, 13, 9, 7, 11};

// ------- Contants e ----------
int e[] = {6, 5, 10, 2, 12, 15, 3, 11, 8, 16, 4, 7, 1, 9, 14, 13};

// ------- Contants h ----------
int h[] = {10, 12, 15, 5, 11, 2, 3, 13, 6, 16, 1, 8, 14, 4, 7, 9};

// ------- Contants v ----------
int v[] = {7, 10, 14, 16, 15, 2, 13, 11, 8, 6, 1, 4, 3, 5, 9, 12};

// ------- Contants w ----------
int w[] = {6, 9, 12, 16, 11, 4, 2, 14, 5, 1, 15, 8, 3, 7, 10, 13};

/*************************************************
* Names:        Mapd, MapDelta, MapD
*
* Description: Elementary Diffusion Functions.
*
* Arguments:   - uint32_t x: a 32 bits word
*              - int v: the right shift parameter
*              - int w: the left shift parameter
*
* Returns a 32 bits word (success)
**************************************************/
uint32_t Mapd(uint32_t x, int v, int w)
{
    return (~x >> v) ^ (x << w) ^ x;
}

///////////////////////////////
// Round constants Rc
//////////////////////////////
uint32_t Constants_Rc(uint32_t Constants_Rc[], int z)
{
    int delta[8];
    uint32_t delta_;
    for (int r = 0; r <= z - 1; r++)
    {

        delta[0] = (v[r % 16] + a[r % 15]) % 16;
        delta[1] = (v[(r + 1) % 16] + b[r % 15]) % 16;
        delta[2] = (v[(r + 2) % 16] + e[r % 15]) % 16;
        delta[3] = (v[(r + 3) % 16] + h[r % 15]) % 16;
        delta[4] = (w[r % 16] + a[(r + 1) % 15]) % 16;
        delta[5] = (w[(r + 1) % 16] + b[(r + 1) % 15]) % 16;
        delta[6] = (w[(r + 2) % 16] + e[(r + 1) % 15]) % 16;
        delta[7] = (w[(r + 3) % 16] + h[(r + 1) % 15]) % 16;
        delta_ = delta[7] + delta[6] * (1 << 4) + delta[5] * (1 << 8) + delta[4] * (1 << 12) + delta[3] * (1 << 16) + delta[2] * (1 << 20) + delta[1] * (1 << 24) + delta[0] * (1 << 28);
        Constants_Rc[r] = Mapd((delta_ + ((uint32_t)r)), v[r % 16], w[r % 15]);
    }
}

void main()
{
    uint32_t roundconstants_Rc[26];

    printf("// ------- Contants Rc ----------\n");
    Constants_Rc(roundconstants_Rc, 26);
    for (int r = 0; r <= 25; r++)
    {
        printf("0x%08x, ", roundconstants_Rc[r]);
    }
    //----------------------------------------
}