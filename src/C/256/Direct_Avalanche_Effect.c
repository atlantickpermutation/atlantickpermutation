/** Avalanche Effect Study of the Direct Atlantick Permutation
 * Implementors: Atlantick Team
 * License: CC0, attribution kindly requested. Blame taken too,
 * but not liability.
 */

////////////////////////////////////////////////////////////////////////////////////////////
//  This implementation is for q=1 (One SuperBlock), 32-bits words and a block size of 256 //
////////////////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include "sha3/sha3.c"

// Round Constantes (for 26 rounds at most)
uint32_t Rc[] = {0x100cd809, 0x8cc76981, 0xbbfdf360, 0xef5210ad, 0x14deb5c6, 0x035e2f35, 0x1b5e8771, 0x527b79bf, 0x55b5b64c, 0x463314df, 0x3aee3d1a, 0x05221a65, 0x9a550a2d, 0xee3ed3bb, 0x75e0a6f4, 0x1d12d04e, 0x7f895ac6, 0x86a82676, 0x49855833, 0xf96a8950, 0x9b15d051, 0x004381ad, 0xc58fe08b, 0xa40763e1, 0x29254a93, 0xe4cf79bb};

// Left and Right Rotation constants

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

///////////////////////////////// Permutation /////////////////

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

uint32_t MapDelta(uint32_t x, int v, int w)
{
    return (x >> v) ^ (x << w) ^ x;
}

uint32_t MapD(uint32_t x, int v, int w)
{
    return (x >> v) ^ (x << w);
}

/*************************************************
* Names:        AtlantickLinearLayer
*
* Description: Atlantick Linear.
*
* Arguments:   - uint32_t State: a block (point to an array) 
*                       of 2*p 32 bits words
*              - int p: a parameter that defines the
*                       block size              
*
* Returns a pointer to the output block of length 2*p 32 bits
*         words(success)
**************************************************/
uint32_t *AtlantickLinearLayer(uint32_t *State, int p)
{
    uint32_t C, l = 2, temp[l + 1];

    ///////////////////////////////
    State[2 * p - 1] = State[2 * p - 1] ^ Mapd(State[0], v[0], w[0]);
    ///////////////////////////////

    ////////////// MIXA in dimension 2p-1
    temp[0] = Mapd(State[2 * p - 1], v[1], w[1]);
    temp[1] = Mapd(State[1], v[2], w[2]);
    temp[2] = State[2];
    for (int i = 3; i <= 2 * p - 2; i++)
    {
        temp[2] ^= State[i];
    }

    temp[0] = temp[0] ^ temp[1] ^ temp[2];
    temp[0] = Mapd(temp[0], v[3], w[3]);

    State[0] = State[0] ^ temp[0];
    ///////////////////////////

    ////////////// MIXA in dimension 2
    temp[0] = Mapd(State[0], v[4], w[4]);
    temp[1] = Mapd(State[2 * p - 1], v[5], w[5]);

    C = temp[0] ^ temp[1];
    for (int j = 0; j < l; j++)
    {
        temp[j] = MapDelta(C, v[j % 16], w[j % 16]);
    }

    // Sequential implementation (Compute in parallel if wanted)
    int j, u;
    for (int k = 1; k <= 2 * p - 2; k++)
    {
        j = k % l;
        u = (k - j) / l;

        C = MapDelta(temp[j], e[u % 16], a[u % 15]);
        State[k] = State[k] ^ C;
    }
    ///////////////////////////

    ////////////// MIXA in dimension 2p-1
    temp[0] = Mapd(State[2 * p - 2], v[6], w[6]);
    temp[1] = Mapd(State[0], v[7], w[7]);
    temp[2] = State[2];
    for (int i = 3; i <= 2 * p - 2; i++)
    {
        temp[2] ^= State[i];
    }

    temp[0] = temp[0] ^ temp[1] ^ temp[2];
    temp[0] = Mapd(temp[0], v[8], w[8]);

    State[2 * p - 1] = State[2 * p - 1] ^ temp[0];
    ///////////////////////////

    ///////////////////////////
    State[0] = State[0] ^ Mapd(State[2 * p - 1], v[9], w[9]);
    ///////////////////////////

    return State;
}

/*************************************************
* Names:        NonLinear
*
* Description: Elementary non linear map.
*
* Arguments:   - uint32_t tempB: a block (point to an array) 
*                       of p 32 bits words
*              - int p: a parameter that defines the
*                       block size              
*
* Returns 0 (success)
**************************************************/
int NonLinear(uint32_t *tempB, int p)
{
    uint32_t tmp1 = tempB[0], tmp2 = tempB[1];
    for (int i = 0; i <= p - 3; i++)
    {
        tempB[i] = (~tempB[i + 1]) & tempB[i + 2];
    }
    tempB[p - 2] = (~tempB[p - 1]) & tmp1;
    tempB[p - 1] = (~tmp1) & tmp2;
    return 0;
}

/*************************************************
* Names:        MIXASTAR
*
* Description: MIXA* map.
*
* Arguments:   - uint32_t tempR: a block (point to an array) 
*                       of p 32 bits words
*              - int p: a parameter that defines the
*                       block size              
*
* Returns a pointer to the output block of length p 32 bits
*         words(success)
**************************************************/
uint32_t *MIXASTAR(uint32_t *tempR, int p)
{
    uint32_t C, temp[3], tempB[p], *tempo;
    ////////////// MIXA STAR
    temp[0] = Mapd(tempR[p - 1], v[10], w[10]);
    temp[1] = Mapd(tempR[0], v[11], w[11]);

    temp[2] = tempR[1];

    C = temp[0] ^ temp[1] ^ temp[2];
    temp[0] = MapD(C, v[0], w[0]);

    // Sequential implementation (Parallel implementation possible)
    for (int k = 0; k <= p - 1; k++)
    {
        tempB[k] = tempR[k] ^ MapD(temp[0], b[k % 16], h[k % 15]);
    }

    tempo = &tempB[0];
    return tempo;
}

/*************************************************
* Names:        MIXB
*
* Description: MIXB map.
*
* Arguments:   - uint32_t tempR: a block (point to an array) 
*                       of p 32 bits words
*              - int p: a parameter that defines the
*                       block size              
*
* Returns a pointer to the output block of length p 32 bits
*         words(success)
**************************************************/
uint32_t *MIXB(uint32_t *tempR, int p)
{
    uint32_t *tempB;

    //////////////////// MIXA*
    tempB = MIXASTAR(tempR, p);
    //////////////////////////

    ////////////// NonLinear
    NonLinear(tempB, p);
    //////////////////////////

    return tempB;
}

/*************************************************
* Names:        SumBlock
*
* Description: Sum of two 32 bits words Block of size p.
*
* Arguments:   - uint32_t input1: a block (point to an array) 
*                       of p 32 bits words
*              - uint32_t input2: a block (point to an array) 
*                       of p 32 bits words
*              - int p: a parameter that defines the
*                       block size              
*
* Returns 0 (success)
**************************************************/
int SumBlock(uint32_t input1[], uint32_t input2[], int p)
{
    for (int i = 0; i < p; i++)
    {
        input1[i] = input1[i] ^ input2[i];
    }
    return 0;
}

/*************************************************
* Names:        RotBlock
*
* Description: RotBlock map.
*
* Arguments:   - uint32_t tempB: a block (point to an array) 
*                       of p 32 bits words
*              - int p: a parameter that defines the
*                       block size             
*
* Returns a pointer to the output block of length p 32 bits
*         words(success)
**************************************************/
uint32_t *RotBlock(uint32_t *tempB, int p)
{
    uint32_t tmp[p], *output;
    tmp[p - 1] = tempB[0];
    for (int i = 0; i <= p - 2; i++)
        tmp[i] = tempB[i + 1];

    output = &tmp[0];
    return output;
}

/*************************************************
* Names:        AtlantickSPbox
*
* Description: Atlantick SPbox.
*
* Arguments:   - uint32_t State: a block (point to an array) 
*                       of 2*p 32 bits words
*              - int p: a parameter that defines the
*                       block size              
*
* Returns a pointer to the output block of length 2*p 32 bits
*         words(success)
**************************************************/
uint32_t *AtlantickSPbox(uint32_t *State, int p)
{
    uint32_t tempL[p], tempR[p];
    // Sequential computation (Compute in parallel if wanted)

    for (int i = 0; i <= p - 1; i++)
    {
        tempL[i] = State[i];
        tempR[i] = State[p + i];
    }

    SumBlock(tempR, RotBlock(tempL, p), p);
    SumBlock(tempL, MIXB(tempR, p), p);
    SumBlock(tempR, MIXB(tempL, p), p);
    SumBlock(tempL, RotBlock(tempR, p), p);

    for (int i = 0; i <= p - 1; i++)
    {
        State[i] = tempL[i];
        State[p + i] = tempR[i];
    }
    return State;
}

/*************************************************
* Names:        AtlantickRoundFunction
*
* Description: Atlantick Round Permutation.
*
* Arguments:   - uint32_t M: a block (point to an array) 
*                       of 2*p 32 bits words
*              - int p: a parameter that defines the
*                       block size  
*              - int r: the current round
*
* Returns a pointer to the output block of length 2*p 32 bits
*         words(success)
**************************************************/
uint32_t *AtlantickRoundFunction(uint32_t *M, int p, int r)
{
    M = AtlantickLinearLayer(M, p);
    M[0] = M[0] ^ Rc[2 * r];
    M[p] = M[p] ^ (uint32_t)p ^ Rc[2 * r + 1];
    M = AtlantickSPbox(M, p);
    return M;
}

/*************************************************
* Names:        AtlantickPermutation
*
* Description: Atlantick Permutation.
*
* Arguments:   - uint32_t M: a block (point to an array) 
*                       of 2*p 32 bits words
*              - int p: a parameter that defines the
*                       block size
*              - int s: the number of rounds    
*
* Returns a pointer to the output block of length 2*p 32 bits
*         words(success)
**************************************************/
uint32_t *AtlantickPermutation(uint32_t *M, int p, int s)
{
    for (int r = 0; r <= s - 1; r++)
    {
        M = AtlantickRoundFunction(M, p, r);
    }
    return M;
}

///////////////////////////////////////////////////////////////

/*************************************************
* Names:        _Sj
*
* Description: Generate a 2*p*32 bits block with all
*               bits to 0 except bit at position j.
*
* Arguments:   - int j: an integer between 0 and 2*p*32-1
*                       block size
*              - uint32_t M: a block (point to an array) 
*                       of 2*p 32 bits words
*              - int p: a parameter that defines the
*                       block size 
*
* Returns 0 (success)
**************************************************/
int _Sj(int j, uint32_t M[], int p)
{
    for (int i = 0; i <= 2 * p - 1; i++)
    {
        if ((j >= (i * 32)) && (j < ((i + 1) * 32)))
        {
            M[i] = (1 << ((j % 32)));
        }
        else
        {
            M[i] = 0;
        }
    }
    return 0;
}

/*************************************************
* Names:        SumState
*
* Description: Sum of two 32 bits words Block.
*
* Arguments:   - uint32_t input1: a block (point to an array) 
*                       of "size" 32 bits words
*              - uint32_t input2: a block (point to an array) 
*                       of "size" 32 bits words
*              - int size: a parameter that defines the
*                       block size              
*
* Returns 0 (success)
**************************************************/
int SumState(uint32_t input1[], uint32_t input2[], int size)
{
    for (int i = 0; i < size; i++)
    {
        input1[i] = input1[i] ^ input2[i];
    }
    return 0;
}

/*************************************************
* Names:        hwt
*
* Description: Computes the Hamming weight of a 32 bits word.
*
* Arguments:   - uint32_t n: a 32 bits word          
*
* Returns the Hamming weight of n (success)
**************************************************/
int hwt(uint32_t n)
{
    return (n & 1) + (n == 0 ? 0 : hwt(n >> 1));
}

/*************************************************
* Names:        block_hwt
*
* Description: Computes the Hamming weight of a 32 bits word 
*              block of size "size".
*
* Arguments:   - uint32_t block: a block (point to an array) 
*                       of "size" 32 bits words
*              - int size: a parameter that defines the
*                       block size      
*
* Returns the Hamming weight of block (success)
**************************************************/
int block_hwt(uint32_t block[], int size)
{
    int out = 0;
    for (int i = 0; i <= size - 1; i++)
    {
        out += hwt(block[i]);
    }
    return out;
}

/*************************************************
* Names:        SumBlockP
*
* Description: Sum of two 32 bits words Block.
*
* Arguments:   - uint32_t input1: a block (point to an array) 
*                       of "size" 32 bits words
*              - uint32_t input2: a block (point to an array) 
*                       of "size" 32 bits words
*              - uint32_t output: a block (point to an array) 
*                       of "size" 32 bits words
*              - int size: a parameter that defines the
*                       block size              
*
* Returns 0 (success)
**************************************************/
int SumBlockP(uint32_t input1[], uint32_t input2[], uint32_t output[], int size)
{
    for (int i = 0; i < size; i++)
    {
        output[i] = input1[i] ^ input2[i];
    }
    return 0;
}

/*************************************************
* Names:        convert_32_to_8_bit
*
* Description: Convert a 32 bits word to 8 bits  words
*
* Arguments:   - uint32_t x: the 32 bits word to convert
*              - int out: a block (point to an array) 
*                       8 bits words to output     
*
* Returns 0 (success)
**************************************************/
int convert_32_to_8_bit(uint32_t x, uint8_t out[])
{
    for (int i = 0; i < 3; i++)
    {
        out[i] = (uint8_t)(x / (0x10 << (4 * (5 - 2 * i))));
        x = x % (0x10 << (4 * (5 - 2 * i)));
    }
    out[3] = x;
    return 0;
}

/*************************************************
* Names:        convert_array_32_to_8_bit
*
* Description: Convert a 32 bits words block to 8 bits 
*               words block
*
* Arguments:   - uint32_t a: a block (point to an array) 
*                       of size_a 32 bits words
*              - int size_a: a parameter that defines the
*                       block size 
*              - uint32_t out: a block (point to an array) 
*                       of 8 bits words
*
* Returns a pointer to the output block of length 8 bits
*         words (success)
**************************************************/
uint8_t *convert_array_32_to_8_bit(uint32_t *a, int size_a, uint8_t out[])
{

    for (int i = 0; i <= size_a - 1; i++)
    {
        uint8_t tmp[4];
        convert_32_to_8_bit(a[i], tmp);
        for (int j = 0; j <= 3; j++)
        {
            out[i * 4 + j] = tmp[j];
        }
    }

    return out;
}

/*************************************************
* Names:        convert_array_8_to_32_bit
*
* Description: Convert an 8 bits words block to 32 bits 
*               words block
*
* Arguments:   - uint32_t a: a block (point to an array) 
*                       of size_a 8 bits words
*              - int size_a: a parameter that defines the
*                       block size 
*              - uint32_t out: a block (point to an array) 
*                       of 32 bits words
*
* Returns 0 (success)
**************************************************/
int convert_array_8_to_32_bit(uint8_t *a, int size_a, uint32_t out[])
{

    uint32_t tmp = 0;
    int out_size = (size_a + 3) / 4;
    for (int i = 0; i <= size_a - 1; i++)
    {
        if ((i % 4) == 0)
        {
            tmp += a[size_a - i - 1];
            if (i == (size_a - 1))
            {
                out[0] = tmp;
            }
        }
        else if ((i % 4) == 1)
        {
            tmp += a[size_a - i - 1] * (0x10 << 4);
            if (i == (size_a - 1))
            {
                out[out_size - i / 4 - 1] = tmp;
                tmp = 0;
            }
        }
        else if ((i % 4) == 3)
        {
            tmp += a[size_a - i - 1] * (0x10 << 20);
            out[out_size - i / 4 - 1] = tmp;
            tmp = 0;
        }
        else
        {
            tmp += a[size_a - i - 1] * (0x10 << 12);
            if (i == (size_a - 1))
            {
                out[0] = tmp;
            }
        }
    }
    return 0;
}

/*************************************************
* Names:        add_t
*
* Description: Add a 32 bits word to a block
*
* Arguments:   - uint32_t Z: a block (point to an array) 
*                       of size_z 32 bits words
*              - int size_z: a parameter that defines the
*                       block size 
*              - uint32_t t: a 32 bits word
*              - uint32_t Zn: a block (point to an array) 
*                       of size_z+1 32 bits words
*
* Returns 0 (success)
**************************************************/
int add_t(uint32_t Z[], int size_z, uint32_t t, uint32_t Zn[])
{
    for (int i = 0; i <= size_z - 1; i++)
    {
        Zn[i] = Z[i];
    }
    Zn[size_z] = t;
    return 0;
}

void main()
{
    int round = 7;
    double result[round][5];

    int p = 4, size = 2 * p * 32 / 8;
    uint32_t Z[8] = {0}, Sj[2 * p], S[2 * p], tempo[2 * p];
    uint8_t tmp[size], out[size], temp[size + 4];
    int X[2 * p * 32][1000];

    AtlantickRoundFunction(Z, p, 0);
    sha3_256(out, size, convert_array_32_to_8_bit(Z, 2 * p, tmp), 2 * p * 32 / 8);
    convert_array_8_to_32_bit(out, size, Z);

    for (int s = 1; s <= round; s++)
    {
        result[s - 1][3] = 2 * p * 32;
        result[s - 1][4] = 0;

        double arith_mean = 0, quad_mean = 0, std = 0;
        for (int j = 0; j <= 2 * p * 32 - 1; j++)
        {
            _Sj(j, Sj, p);
            for (int t = 0; t <= 999; t++)
            {
                uint32_t Zn[2 * p + 1];
                add_t(Z, 2 * p, (uint32_t)t, Zn);
                sha3_256(out, size, convert_array_32_to_8_bit(Zn, 2 * p + 1, temp), (2 * p + 1) * 32 / 8);
                convert_array_8_to_32_bit(out, size, S);

                SumBlockP(S, Sj, tempo, 2 * p);

                AtlantickPermutation(S, p, s);

                AtlantickPermutation(tempo, p, s);

                SumBlockP(S, tempo, S, 2 * p);

                X[j][t] = block_hwt(S, 2 * p);

                arith_mean += (double)(X[j][t]) / ((double)2 * p * 32 * 1000);
                quad_mean += pow((double)X[j][t], 2) / ((double)2 * p * 32 * 1000);
                if (result[s - 1][3] > X[j][t])
                {
                    result[s - 1][3] = X[j][t];
                }
                if (result[s - 1][4] < X[j][t])
                {
                    result[s - 1][4] = X[j][t];
                }
            }
        }
        quad_mean = sqrt((double)quad_mean);
        for (int j = 0; j <= 2 * p * 32 - 1; j++)
        {
            for (int t = 0; t <= 999; t++)
            {
                std += pow((double)(X[j][t] - arith_mean), 2) / ((double)2 * p * 32 * 1000);
            }
        }
        std = sqrt((double)std);

        result[s - 1][0] = arith_mean;
        result[s - 1][1] = quad_mean;
        result[s - 1][2] = std;
        printf("\nRounds\tarith_mean\tquad_mean\tstd\tmin\tmax\n");
        printf("%d\t%.2lf\t%14.2lf\t%12.2lf\t%.0lf\t%.0lf\n", s, result[s - 1][0], result[s - 1][1], result[s - 1][2], result[s - 1][3], result[s - 1][4]);
    }
    printf("\n\nRounds\tarith_mean\tquad_mean\tstd\tmin\tmax\n");
    for (int s = 1; s <= round; s++)
    {
        printf("%d\t%.2lf\t%14.2lf\t%12.2lf\t%.0lf\t%.0lf\n", s, result[s - 1][0], result[s - 1][1], result[s - 1][2], result[s - 1][3], result[s - 1][4]);
    }
}