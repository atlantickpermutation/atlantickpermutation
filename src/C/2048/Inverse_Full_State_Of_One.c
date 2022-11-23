/** Full State of One Study of the Inverse Atlantick Permutation
 * Implementors: Atlantick Team
 * License: CC0, attribution kindly requested. Blame taken too,
 * but not liability.
 */

//////////////////////////////////////////////////////////////////////////////////////////////
//  This implementation is for q=1 (One SuperBlock), 32-bits words and a block size of 2048  //
//////////////////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdint.h>
#include "sha3/sha3.c"
#include "Permutation.c"

/*************************************************
* Names:        S
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
uint32_t S(int j, uint32_t M[], int p)
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

///////////////////////////////// Permutation /////////////////

/*************************************************
* Names:        MapdOR, MapDeltaOR, MapDOR
*
* Description: Elementary Diffusion Functions.
*
* Arguments:   - uint32_t x: a 32 bits word
*              - int v: the right shift parameter
*              - int w: the left shift parameter
*
* Returns a 32 bits word (success)
**************************************************/
uint32_t MapdOR(uint32_t x, int v, int w)
{
    return (~x >> v) | (x << w) | x;
}
uint32_t MapDeltaOR(uint32_t x, int v, int w)
{
    return (x >> v) | (x << w) | x;
}
uint32_t MapDOR(uint32_t x, int v, int w)
{
    return (x >> v) | (x << w);
}

/*************************************************
* Names:        InvAtlantickLinearLayerOR
*
* Description: Inverse Atlantick Linear replacing ^ by |.
*
* Arguments:   - uint32_t State: a block (point to an array) 
*                       of 2*p 32 bits words
*              - int p: a parameter that defines the
*                       block size              
*
* Returns a pointer to the output block of length 2*p 32 bits
*         words(success)
**************************************************/
uint32_t *InvAtlantickLinearLayerOR(uint32_t *State, int p)
{
    uint32_t C, l = 2, temp[l + 1];

    ///////////////////////////
    State[0] = State[0] | MapdOR(State[2 * p - 1], v[9], w[9]);
    ///////////////////////////

    ////////////// MIXA in dimension 2p-1
    temp[0] = MapdOR(State[2 * p - 2], v[6], w[6]);
    temp[1] = MapdOR(State[0], v[7], w[7]);
    temp[2] = 0;
    for (int i = 1; i <= 2 * p - 3; i++)
    {
        temp[2] |= State[i];
    }

    temp[0] = temp[0] | temp[1] | temp[2];
    temp[0] = MapdOR(temp[0], v[8], w[8]);

    State[2 * p - 1] = State[2 * p - 1] | temp[0];
    ///////////////////////////

    ////////////// MIXA in dimension 2
    temp[0] = MapdOR(State[0], v[4], w[4]);
    temp[1] = MapdOR(State[2 * p - 1], v[5], w[5]);

    C = temp[0] | temp[1];
    for (int j = 0; j < l; j++)
    {
        temp[j] = MapDeltaOR(C, v[j % 16], w[j % 16]);
    }

    // Sequential implementation (Compute in parallel if wanted)
    int j, u;
    for (int k = 2 * p - 2; k >= 1; k--)
    {
        j = k % l;
        u = (k - j) / l;

        C = MapDeltaOR(temp[j], e[u % 16], a[u % 15]);
        State[k] = State[k] | C;
    }
    ///////////////////////////

    ////////////// MIXA in dimension 2p-1
    temp[0] = MapdOR(State[2 * p - 1], v[1], w[1]);
    temp[1] = MapdOR(State[1], v[2], w[2]);
    temp[2] = State[2];
    for (int i = 3; i <= 2 * p - 2; i++)
    {
        temp[2] |= State[i];
    }

    temp[0] = temp[0] | temp[1] | temp[2];
    temp[0] = MapdOR(temp[0], v[3], w[3]);

    State[0] = State[0] | temp[0];
    ///////////////////////////

    ///////////////////////////////
    State[2 * p - 1] = State[2 * p - 1] | MapdOR(State[0], v[0], w[0]);
    ///////////////////////////////

    return State;
}

/*************************************************
* Names:        NonLinearOR
*
* Description: Elementary non linear map replacing ^ by |.
*
* Arguments:   - uint32_t tempB: a block (point to an array) 
*                       of p 32 bits words
*              - int p: a parameter that defines the
*                       block size              
*
* Returns 0 (success)
**************************************************/
int NonLinearOR(uint32_t *tempB, int p)
{
    uint32_t tmp1 = tempB[0], tmp2 = tempB[1];
    for (int i = 0; i <= p - 3; i++)
    {
        tempB[i] = (~tempB[i + 1]) | tempB[i + 2];
    }
    tempB[p - 2] = (~tempB[p - 1]) | tmp1;
    tempB[p - 1] = (~tmp1) | tmp2;
    return 0;
}

/*************************************************
* Names:        MIXASTAROR
*
* Description: MIXA* map replacing ^ by |.
*
* Arguments:   - uint32_t tempR: a block (point to an array) 
*                       of p 32 bits words
*              - int p: a parameter that defines the
*                       block size              
*
* Returns a pointer to the output block of length p 32 bits
*         words(success)
**************************************************/
uint32_t *MIXASTAROR(uint32_t *tempR, int p)
{
    uint32_t C, temp[3], tempB[p], *tempo;
    ////////////// MIXA STAR
    temp[0] = MapdOR(tempR[p - 1], v[10], w[10]);
    temp[1] = MapdOR(tempR[0], v[11], w[11]);

    C = temp[0] | temp[1];
    temp[0] = MapDOR(C, v[0], w[0]);

    // Sequential implementation (Parallel implementation possible)
    for (int k = 0; k <= p - 1; k++)
    {
        tempB[k] = tempR[k] | MapDOR(temp[0], b[k % 16], h[k % 15]);
    }

    tempo = &tempB[0];
    return tempo;
}

/*************************************************
* Names:        MIXBOR
*
* Description: MIXB map replacing ^ by |.
*
* Arguments:   - uint32_t tempR: a block (point to an array) 
*                       of p 32 bits words
*              - int p: a parameter that defines the
*                       block size              
*
* Returns a pointer to the output block of length p 32 bits
*         words(success)
**************************************************/
uint32_t *MIXBOR(uint32_t *tempR, int p)
{
    uint32_t *tempB;

    //////////////////// MIXA*
    tempB = MIXASTAROR(tempR, p);
    //////////////////////////

    ////////////// NonLinearOR
    NonLinearOR(tempB, p);

    return tempB;
}

/*************************************************
* Names:        SumBlockOR
*
* Description: Sum of two 32 bits words Block of size p, 
*              replacing ^ by |.
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
int SumBlockOR(uint32_t input1[], uint32_t input2[], int p)
{
    for (int i = 0; i < p; i++)
    {
        input1[i] = input1[i] | input2[i];
    }
    return 0;
}

/*************************************************
* Names:        InvAtlantickSPboxOR
*
* Description: Inverse Atlantick SPbox replacing ^ by |.
*
* Arguments:   - uint32_t State: a block (point to an array) 
*                       of 2*p 32 bits words
*              - int p: a parameter that defines the
*                       block size              
*
* Returns a pointer to the output block of length 2*p 32 bits
*         words(success)
**************************************************/
uint32_t *InvAtlantickSPboxOR(uint32_t *State, int p)
{
    uint32_t tempL[p], tempR[p];
    // Sequential computation (Compute in parallel if wanted)
    for (int i = 0; i <= p - 1; i++)
    {
        tempL[i] = State[i];
        tempR[i] = State[p + i];
    }

    SumBlockOR(tempL, RotBlock(tempR, p), p);
    SumBlockOR(tempR, MIXBOR(tempL, p), p);
    SumBlockOR(tempL, MIXBOR(tempR, p), p);
    SumBlockOR(tempR, RotBlock(tempL, p), p);

    for (int i = 0; i <= p - 1; i++)
    {
        State[i] = tempL[i];
        State[p + i] = tempR[i];
    }
    return State;
}

/*************************************************
* Names:        InvAtlantickRoundFunctionOR
*
* Description: Inverse Atlantick Permutation replacing 
*               ^ by |.
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
uint32_t *InvAtlantickRoundFunctionOR(uint32_t *M, int p, int r)
{
    M = InvAtlantickSPboxOR(M, p);
    M[p] = M[p] | (uint32_t)p | Rc[2 * r + 1];
    M[0] = M[0] | Rc[2 * r];
    M = InvAtlantickLinearLayerOR(M, p);
    return M;
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

///////////////////////////////////////////////////////////////

void main()
{
    int p = 32, round = 11;
    uint32_t M[2 * p], Z[2 * p], F[2 * p], N[2 * p], Sj[2 * p], temp[2 * p];
    int Lr[2 * p * 32][100], Lrm[2 * p * 32];
    int size = 2 * p * 32 / 8;
    uint8_t tmp[size], out[size], tempo[size + 4];

    AtlantickRoundFunction(F, p, 0);
    sha3_2048(out, size, convert_array_32_to_8_bit(F, 2 * p, tmp), 2 * p * 32 / 8);
    convert_array_8_to_32_bit(out, size, F);

    //Initializing Z to 0^(2pn)
    for (int i = 0; i <= 2 * p - 1; i++)
    {
        Z[i] = 0;
    }

    //Setting of Z to 1^(2pn)
    for (int j = 0; j <= 2 * p * 32 - 1; j++)
    {
        S(j, M, p);
        SumState(Z, M, 2 * p);
    }

    for (int j = 0; j <= 2 * p * 32 - 1; j++)
    {
        S(j, Sj, p);
        SumState(Sj, F, 2 * p);
        AtlantickRoundFunction(Sj, p, 0);
        uint32_t Zn[2 * p + 1];
        Lrm[j] = 0;
        for (int t = 0; t <= 99; t++)
        {
            add_t(Sj, 2 * p, (uint32_t)t, Zn);
            sha3_2048(out, size, convert_array_32_to_8_bit(Zn, 2 * p + 1, tempo), (2 * p + 1) * 32 / 8);
            convert_array_8_to_32_bit(out, size, M);

            int r = 0;
            InvAtlantickRoundFunctionOR(M, p, r);
            SumBlockP(M, Z, temp, 2 * p);
            int hw = block_hwt(temp, 2 * p);

            while ((hw != 0) && (r < round))
            {
                r += 1;
                InvAtlantickRoundFunctionOR(M, p, r);
                SumBlockP(M, Z, temp, 2 * p);
                int hw = block_hwt(temp, 2 * p);
            }
            Lr[j][t] = r + 1;
            if (r + 1 > Lrm[j])
            {
                Lrm[j] = r + 1;
            }
        }
    }

    printf("N° Active Bit\t|\tNbr of Rounds\n");
    int max = 0;
    for (int i = 0; i <= 2 * p * 32 - 1; i++)
    {
        printf("%10d\t\t|\t %d \n", i + 1, Lrm[i]);

        if (Lrm[i] > max)
        {
            max = Lrm[i];
        }
    }
    printf("Maximum Number of Rounds ==> %d\n", max);
}